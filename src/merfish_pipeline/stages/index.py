"""``index`` stage — scan raw data and produce manifest + normalized positions.

This stage creates two canonical files consumed by all downstream stages:

1. ``manifest.csv`` — one row per raw image, with columns:
   ``round, fov, z_slice, channel, wavelength, abs_path, file_size, image_shape``

2. ``positions.normalized.csv`` — microscope-agnostic position data, with columns:
   ``round, tile_number, stage_pos_x, stage_pos_y, z_position_0, ..., z_position_N``
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from merfish_pipeline.microscopes import get_adapter
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage


@register_stage("index")
class IndexStage(PipelineStage):
    """Scan raw data → ``manifest.csv`` + ``positions.normalized.csv``."""

    description = "Index raw data and generate manifest"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        errors: list[str] = []
        raw_dir = Path(self.config.paths.raw_data_dir)
        if not raw_dir.exists():
            errors.append(f"Raw data directory does not exist: {raw_dir}")
        elif not raw_dir.is_dir():
            errors.append(f"Raw data path is not a directory: {raw_dir}")
        return errors

    def check_outputs_exist(self) -> bool:
        out = self.get_output_dir()
        manifest = out / "manifest.csv"
        positions = out / "positions.normalized.csv"
        return manifest.exists() and positions.exists()

    def run(self, dry_run: bool = False) -> StageResult:
        output_dir = self.get_output_dir()
        output_dir.mkdir(parents=True, exist_ok=True)

        raw_dir = Path(self.config.paths.raw_data_dir)
        adapter = get_adapter(self.config)

        # --- Build manifest ---
        self.logger.info("Scanning raw data in %s ...", raw_dir)
        file_records = adapter.discover_raw_files(raw_dir)

        if not file_records:
            return StageResult(
                status="failed",
                error=f"No raw files discovered in {raw_dir}",
            )

        manifest_df = pd.DataFrame(file_records)

        # Ensure canonical column order
        canonical_cols = [
            "round", "fov", "z_slice", "channel", "wavelength",
            "abs_path", "file_size",
        ]
        # Keep extra columns (like image_shape) if present
        extra_cols = [c for c in manifest_df.columns if c not in canonical_cols]
        manifest_df = manifest_df[
            [c for c in canonical_cols if c in manifest_df.columns] + extra_cols
        ]

        manifest_path = output_dir / "manifest.csv"
        manifest_df.to_csv(manifest_path, index=False)
        self.logger.info("Wrote manifest: %s (%d rows)", manifest_path, len(manifest_df))

        # --- Build normalized positions ---
        self.logger.info("Reading position data...")
        positions_df = adapter.read_positions(raw_dir)

        positions_path = output_dir / "positions.normalized.csv"
        if not positions_df.empty:
            positions_df.to_csv(positions_path, index=False)
            self.logger.info(
                "Wrote positions: %s (%d rows)", positions_path, len(positions_df)
            )
        else:
            # Write an empty file with header so downstream stages don't fail
            positions_df = pd.DataFrame(
                columns=["round", "tile_number", "stage_pos_x", "stage_pos_y"]
            )
            positions_df.to_csv(positions_path, index=False)
            self.logger.warning("No position data found — wrote empty positions file.")

        # --- Summary ---
        n_rounds = manifest_df["round"].nunique()
        n_fovs = manifest_df["fov"].nunique()
        n_z = manifest_df["z_slice"].nunique()
        n_channels = manifest_df["channel"].nunique()
        total_size_mb = manifest_df["file_size"].sum() / (1024 * 1024)

        self.logger.info(
            "Index summary: %d rounds, %d FOVs, %d z-slices, %d channels, %.1f MB total",
            n_rounds, n_fovs, n_z, n_channels, total_size_mb,
        )

        # --- Validate completeness ---
        completeness_issues = self._check_completeness(manifest_df)
        if completeness_issues:
            for issue in completeness_issues[:10]:  # Limit output
                self.logger.warning("Completeness: %s", issue)
            if len(completeness_issues) > 10:
                self.logger.warning(
                    "... and %d more completeness issues.", len(completeness_issues) - 10
                )

        return StageResult(
            status="completed",
            output_files=[str(manifest_path), str(positions_path)],
            metadata={
                "n_rounds": n_rounds,
                "n_fovs": n_fovs,
                "n_z_slices": n_z,
                "n_channels": n_channels,
                "total_files": len(manifest_df),
                "total_size_mb": round(total_size_mb, 2),
                "completeness_issues": len(completeness_issues),
            },
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _check_completeness(df: pd.DataFrame) -> list[str]:
        """Flag missing (round, fov, z) combinations."""
        issues: list[str] = []

        if df.empty:
            return issues

        rounds = sorted(df["round"].unique())
        fovs = sorted(df["fov"].unique())
        z_slices = sorted(df["z_slice"].unique())
        channels = sorted(df["channel"].unique())

        # Build expected set of (round, fov, z, channel) tuples
        existing = set(
            zip(df["round"], df["fov"], df["z_slice"], df["channel"])
        )

        for r in rounds:
            for f in fovs:
                for z in z_slices:
                    for ch in channels:
                        if (r, f, z, ch) not in existing:
                            issues.append(
                                f"Missing: round={r}, fov={f}, z={z}, channel={ch}"
                            )

        return issues
