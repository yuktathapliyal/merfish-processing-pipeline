"""``convert`` stage -- merge per-channel TIFF planes into stacked MERlin format.

Used for ONI and NIKON microscopes where raw acquisition produces one TIFF
file per (round, fov, z-slice, wavelength) combination.  This stage
reassembles those individual planes into a single multi-frame TIFF stack per
(round, fov), with frames interleaved as MERlin expects::

    z0_wv0, z0_wv1, z0_wv2, z1_wv0, z1_wv1, z1_wv2, ...

Algorithm
---------
1. Read ``manifest.csv`` produced by the *index* stage.
2. Group rows by ``(round, fov)`` -- each group becomes one merged stack.
3. For every group, sort by ``z_slice`` then by ``wavelength`` (numerically)
   and build an ordered list of source plane paths.
4. Call :func:`~merfish_pipeline.io.tiff_io.merge_stack` to read all planes
   in parallel and write the merged TIFF.
5. Validate the output with
   :func:`~merfish_pipeline.io.tiff_io.validate_tiff_stack`.

Output filename convention::

    merFISH_merged_{round:02d}_{fov:03d}.tiff
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd

from merfish_pipeline.io.sheet_io import read_sheet
from merfish_pipeline.io.tiff_io import merge_stack, validate_tiff_stack
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_source_dir(config) -> Path:
    """Return the directory that contains the source plane images.

    If reregistration was enabled **and** ``remapped_data_dir`` exists on
    disk, that directory is used (planes have been spatially remapped).
    Otherwise the original ``raw_data_dir`` is returned.
    """
    if config.reregistration.enabled:
        remapped = Path(config.paths.remapped_data_dir)
        if remapped.exists():
            return remapped
    return Path(config.paths.raw_data_dir)


def _build_plane_order(
    group_df: pd.DataFrame,
    wavelengths: list[int],
) -> list[Path]:
    """Return an ordered list of plane paths for a single (round, fov) group.

    The interleave order expected by MERlin is::

        for z in sorted(z_slices):
            for wv in wavelengths:   # numerically sorted
                yield path

    Parameters
    ----------
    group_df:
        Subset of the manifest for one ``(round, fov)`` combination.
        Must contain at least ``z_slice``, ``wavelength``, and ``abs_path``
        columns.
    wavelengths:
        Globally sorted list of unique wavelength values (ints) so that the
        channel order is consistent across all groups.

    Returns
    -------
    list[Path]
        Ordered plane paths ready for :func:`merge_stack`.
    """
    # Index the sub-frame by (z_slice, wavelength) for fast lookup.
    lookup: dict[tuple[int, int], Path] = {}
    for _, row in group_df.iterrows():
        key = (int(row["z_slice"]), int(row["wavelength"]))
        lookup[key] = Path(row["abs_path"])

    z_slices = sorted(group_df["z_slice"].unique())

    plane_paths: list[Path] = []
    for z in z_slices:
        for wv in wavelengths:
            key = (int(z), int(wv))
            if key not in lookup:
                raise KeyError(
                    f"Missing plane for z_slice={z}, wavelength={wv} "
                    f"in (round={group_df['round'].iloc[0]}, "
                    f"fov={group_df['fov'].iloc[0]})"
                )
            plane_paths.append(lookup[key])

    return plane_paths


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("convert")
class ConvertStage(PipelineStage):
    """Merge per-channel TIFF planes into stacked MERlin format."""

    description = "Merge individual TIFF planes into MERlin-compatible stacks"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        errors: list[str] = []

        manifest_path = self._manifest_path()
        if not manifest_path.exists():
            errors.append(
                f"Manifest file not found: {manifest_path}. "
                "Has the 'index' stage been run?"
            )

        source_dir = _resolve_source_dir(self.config)
        if not source_dir.exists():
            errors.append(f"Source data directory does not exist: {source_dir}")

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if the run_metadata.json already exists.

        A more thorough check would verify every expected merged file, but
        the metadata sentinel is sufficient for the skip-if-done logic;
        individual missing files are caught by the force/re-run path.
        """
        metadata_path = self.get_output_dir() / "run_metadata.json"
        return metadata_path.exists()

    def run(self, dry_run: bool = False) -> StageResult:
        start_time = datetime.now()
        stage_output_dir = self.get_output_dir()
        stage_output_dir.mkdir(parents=True, exist_ok=True)

        merlin_data_dir = Path(self.config.paths.merlin_data_dir)
        merlin_data_dir.mkdir(parents=True, exist_ok=True)

        workers = self.config.execution.max_workers

        # ----------------------------------------------------------
        # Step 1: Read manifest
        # ----------------------------------------------------------
        manifest_path = self._manifest_path()
        self.logger.info("Reading manifest: %s", manifest_path)
        manifest_df = read_sheet(manifest_path)

        if manifest_df.empty:
            return StageResult(
                status="failed",
                error=f"Manifest is empty: {manifest_path}",
            )

        # ----------------------------------------------------------
        # Step 1b: Remap source paths if reregistration was used
        # ----------------------------------------------------------
        source_dir = _resolve_source_dir(self.config)
        raw_dir = Path(self.config.paths.raw_data_dir)

        if source_dir != raw_dir:
            self.logger.info(
                "Using remapped data from %s (reregistration enabled)",
                source_dir,
            )
            manifest_df["abs_path"] = manifest_df["abs_path"].apply(
                lambda p: str(Path(str(p).replace(str(raw_dir), str(source_dir))))
            )

        # ----------------------------------------------------------
        # Step 2: Determine global wavelength order
        # ----------------------------------------------------------
        wavelengths = sorted(manifest_df["wavelength"].unique().astype(int))
        self.logger.info(
            "Wavelengths (sorted): %s (%d channels)",
            wavelengths,
            len(wavelengths),
        )

        # ----------------------------------------------------------
        # Step 3: Group by (round, fov) and merge
        # ----------------------------------------------------------
        grouped = manifest_df.groupby(["round", "fov"])
        n_groups = len(grouped)
        self.logger.info("Found %d (round, fov) groups to merge.", n_groups)

        if dry_run:
            self.logger.info("[DRY RUN] Would merge %d stacks.", n_groups)
            return StageResult(
                status="skipped",
                metadata={"dry_run": True, "n_groups": n_groups},
            )

        output_files: list[str] = []
        skipped = 0
        failed: list[str] = []
        force = self.config.pipeline.force

        for idx, ((rnd, fov), group_df) in enumerate(sorted(grouped), start=1):
            rnd_int = int(rnd)
            fov_int = int(fov)
            output_name = f"merFISH_merged_{rnd_int:02d}_{fov_int:03d}.tiff"
            output_path = merlin_data_dir / output_name

            # Skip existing files unless --force
            if output_path.exists() and not force:
                self.logger.debug(
                    "Skipping existing file: %s", output_path
                )
                output_files.append(str(output_path))
                skipped += 1
                continue

            # Build interleaved plane order
            try:
                plane_paths = _build_plane_order(group_df, wavelengths)
            except KeyError as exc:
                self.logger.error("Failed to build plane order: %s", exc)
                failed.append(str(exc))
                continue

            n_planes = len(plane_paths)
            self.logger.info(
                "Merging FOV %d, round %d: %d planes -> %s",
                fov_int,
                rnd_int,
                n_planes,
                output_path,
            )

            # Merge planes into a single stack
            try:
                merge_stack(
                    plane_paths=plane_paths,
                    output_path=output_path,
                    workers=workers,
                )
            except Exception as exc:
                self.logger.error(
                    "Failed to merge FOV %d, round %d: %s",
                    fov_int,
                    rnd_int,
                    exc,
                )
                failed.append(f"round={rnd_int}, fov={fov_int}: {exc}")
                continue

            # Validate the merged stack
            if not validate_tiff_stack(output_path, expected_frames=n_planes):
                msg = (
                    f"Validation failed for {output_path}: "
                    f"expected {n_planes} frames"
                )
                self.logger.error(msg)
                failed.append(msg)
                continue

            output_files.append(str(output_path))

            if idx % 20 == 0:
                self.logger.info(
                    "  Progress: %d / %d groups merged ...", idx, n_groups
                )

        # ----------------------------------------------------------
        # Step 4: Summarise and write metadata
        # ----------------------------------------------------------
        n_merged = len(output_files) - skipped
        self.logger.info(
            "Convert complete: %d merged, %d skipped, %d failed (of %d total).",
            n_merged,
            skipped,
            len(failed),
            n_groups,
        )

        if failed:
            status = "completed" if output_files else "failed"
            error_summary = f"{len(failed)} group(s) failed: {'; '.join(failed[:5])}"
            if len(failed) > 5:
                error_summary += f" ... and {len(failed) - 5} more"
        else:
            status = "completed"
            error_summary = ""

        result = StageResult(
            status=status,
            output_files=output_files,
            metadata={
                "n_groups": n_groups,
                "n_merged": n_merged,
                "n_skipped": skipped,
                "n_failed": len(failed),
                "n_wavelengths": len(wavelengths),
                "wavelengths": wavelengths,
                "source_dir": str(source_dir),
                "merlin_data_dir": str(merlin_data_dir),
            },
            error=error_summary,
        )

        self.write_run_metadata(
            result,
            start_time,
            parameters={
                "max_workers": workers,
                "source_dir": str(source_dir),
                "force": force,
            },
        )

        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _manifest_path(self) -> Path:
        """Return the path to the index-stage manifest."""
        return Path(self.config.paths.output_dir) / "index" / "manifest.csv"
