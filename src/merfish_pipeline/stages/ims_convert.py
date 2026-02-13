"""``ims_convert`` stage -- convert ANDOR IMS (HDF5) files to merged TIFF stacks for MERlin.

This stage is specific to the ANDOR microscope workflow.  ANDOR raw data
arrives as Imaris ``.ims`` (HDF5) files organised in round-specific
directories.  MERlin expects a flat directory of merged TIFF stacks where
each file contains all z-slices interleaved with all channels for a single
FOV, plus a ``stagePos_Round#<n>.csv`` positions file for each round.

Algorithm
---------
1. Scan round folders in ``raw_data_dir`` using the same regex patterns as
   the ANDOR microscope adapter (``1st round``, ``round 1``, ``R1``, etc.).
2. For each round folder, find all ``.ims`` files and extract FOV numbers.
3. For each IMS file:

   a. Read metadata: n_channels, n_z_slices, image_shape.
   b. Apply channel-order mapping: IMS channels are reordered according to
      ``config.raw_data.andor.channel_order`` (default ``[0, 2, 1]``).
   c. Build an output array of shape ``(n_z * n_channels, height, width)``
      with dtype ``uint16``.
   d. For each z-slice, for each mapped channel, read the 2-D image and
      place it at frame index ``z * n_channels + ch_idx``.
   e. Write the merged TIFF as ``merFISH_merged_{round:02d}_{fov:03d}.tiff``.

4. Generate per-round stage-positions CSV files:
   ``stagePos_Round#{round}.csv`` with columns ``tile_number``,
   ``stage_pos_x``, ``stage_pos_y``, ``z_position_0``, ..., ``z_position_N``.
"""

from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from merfish_pipeline.io.hdf5_io import (
    extract_stage_positions,
    read_ims_channel,
    read_ims_metadata,
)
from merfish_pipeline.io.tiff_io import write_tiff
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

# ---------------------------------------------------------------------------
# Round-folder and FOV regex patterns (mirrored from andor.py)
# ---------------------------------------------------------------------------

_ROUND_PATTERNS = [
    re.compile(r"(\d+)(?:st|nd|rd|th)\s*round", re.IGNORECASE),
    re.compile(r"round\s*(\d+)", re.IGNORECASE),
    re.compile(r"R(\d+)", re.IGNORECASE),
]

_FOV_RE = re.compile(r"_F(\d+)\.ims$", re.IGNORECASE)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _find_round_folders(raw_dir: Path) -> dict[int, Path]:
    """Detect round directories under *raw_dir*.

    Returns a mapping of ``{round_number: directory_path}`` for every
    subdirectory whose name matches one of the ANDOR round-naming
    conventions.
    """
    result: dict[int, Path] = {}
    if not raw_dir.is_dir():
        return result
    for entry in sorted(raw_dir.iterdir()):
        if not entry.is_dir():
            continue
        round_num = _get_round_number(entry.name)
        if round_num is not None:
            result[round_num] = entry
    return result


def _get_round_number(folder_name: str) -> int | None:
    """Extract the round number from a folder name.

    Supports patterns like ``1st round``, ``round 1``, and ``R1``.
    """
    for pattern in _ROUND_PATTERNS:
        m = pattern.search(folder_name)
        if m:
            return int(m.group(1))
    return None


def _get_fov_number(ims_path: Path) -> int | None:
    """Extract the FOV number from an IMS filename (``_F<digits>.ims``)."""
    m = _FOV_RE.search(ims_path.name)
    if m:
        return int(m.group(1))
    return None


def _convert_single_ims(
    ims_path: Path,
    output_path: Path,
    channel_order: list[int],
    n_z_override: int | None = None,
) -> dict[str, Any]:
    """Convert one IMS file to a merged TIFF stack.

    Parameters
    ----------
    ims_path:
        Path to the source ``.ims`` file.
    output_path:
        Destination path for the merged TIFF.
    channel_order:
        List mapping MERlin channel index to IMS channel index.
        For example ``[0, 2, 1]`` means MERlin channel 0 reads IMS
        channel 0, MERlin channel 1 reads IMS channel 2, and MERlin
        channel 2 reads IMS channel 1.
    n_z_override:
        If provided, override the number of z-slices read from metadata
        (useful for truncating stacks during testing).

    Returns
    -------
    dict
        Conversion summary with keys ``n_channels``, ``n_z_slices``,
        ``image_shape``, and ``output_path``.
    """
    meta = read_ims_metadata(ims_path)
    n_channels = meta["n_channels"]
    n_z_slices = n_z_override if n_z_override is not None else meta["n_z_slices"]
    height, width = meta["image_shape"]

    # Validate channel_order against available channels
    effective_order = channel_order[:n_channels]
    n_mapped = len(effective_order)

    # Build the merged output array: frames = n_z * n_mapped_channels
    total_frames = n_z_slices * n_mapped
    merged = np.zeros((total_frames, height, width), dtype=np.uint16)

    for z in range(n_z_slices):
        for ch_idx, ims_channel in enumerate(effective_order):
            frame_index = z * n_mapped + ch_idx
            plane = read_ims_channel(ims_path, ims_channel, z)
            merged[frame_index] = plane.astype(np.uint16)

    write_tiff(merged, output_path, bigtiff=True)

    return {
        "n_channels": n_mapped,
        "n_z_slices": n_z_slices,
        "image_shape": (height, width),
        "output_path": str(output_path),
    }


def _generate_positions_csv(
    ims_files: list[tuple[int, Path]],
    output_path: Path,
) -> None:
    """Create a MERlin-compatible stage-positions CSV from IMS metadata.

    Parameters
    ----------
    ims_files:
        List of ``(fov_number, ims_path)`` tuples, sorted by FOV.
    output_path:
        Destination path for the CSV file.
    """
    records: list[dict[str, Any]] = []

    for fov, ims_path in ims_files:
        positions = extract_stage_positions(ims_path)
        record: dict[str, Any] = {
            "tile_number": fov,
            "stage_pos_x": positions["stage_pos_x"],
            "stage_pos_y": positions["stage_pos_y"],
        }
        for i, z_pos in enumerate(positions.get("z_positions", [])):
            record[f"z_position_{i}"] = z_pos
        records.append(record)

    df = pd.DataFrame(records)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


# ---------------------------------------------------------------------------
# Pipeline stage
# ---------------------------------------------------------------------------


@register_stage("ims_convert")
class IMSConvertStage(PipelineStage):
    """Convert ANDOR IMS (HDF5) files to merged TIFF stacks for MERlin consumption."""

    description = "Convert ANDOR IMS files to merged TIFF stacks and generate stage-position CSVs"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        """Check that the raw data directory exists and contains round folders with IMS files."""
        errors: list[str] = []

        raw_dir = Path(self.config.paths.raw_data_dir)
        if not raw_dir.exists():
            errors.append(f"Raw data directory does not exist: {raw_dir}")
            return errors
        if not raw_dir.is_dir():
            errors.append(f"Raw data path is not a directory: {raw_dir}")
            return errors

        round_folders = _find_round_folders(raw_dir)
        if not round_folders:
            errors.append(
                f"No round folders detected in raw data directory: {raw_dir}. "
                f"Expected folder names like '1st round', 'round 1', or 'R1'."
            )
            return errors

        # Check that at least one round folder contains IMS files.
        has_ims = False
        for round_dir in round_folders.values():
            ims_files = [
                p for p in round_dir.iterdir()
                if p.suffix.lower() == ".ims"
            ]
            if ims_files:
                has_ims = True
                break

        if not has_ims:
            errors.append(
                f"No .ims files found in any round folder under: {raw_dir}"
            )

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if merged TIFFs already exist in merlin_data_dir.

        A lightweight check: we verify that the stage metadata file and
        at least one merged TIFF exist.
        """
        metadata_path = self.get_output_dir() / "run_metadata.json"
        if not metadata_path.exists():
            return False

        merlin_dir = Path(self.config.paths.merlin_data_dir)
        if not merlin_dir.exists():
            return False

        merged_tiffs = list(merlin_dir.glob("merFISH_merged_*.tiff"))
        return len(merged_tiffs) > 0

    def run(self, dry_run: bool = False) -> StageResult:
        """Execute IMS-to-TIFF conversion for all rounds and FOVs."""
        start_time = datetime.now()

        raw_dir = Path(self.config.paths.raw_data_dir)
        merlin_dir = Path(self.config.paths.merlin_data_dir)
        max_workers = self.config.execution.max_workers

        # Channel order from config (fall back to [0, 2, 1] if andor config is absent)
        andor_cfg = getattr(self.config.raw_data, "andor", None)
        if andor_cfg is not None:
            channel_order = list(andor_cfg.channel_order)
        else:
            channel_order = [0, 2, 1]

        # Discover round folders
        round_folders = _find_round_folders(raw_dir)

        if dry_run:
            n_ims = sum(
                len([p for p in rd.iterdir() if p.suffix.lower() == ".ims"])
                for rd in round_folders.values()
            )
            self.logger.info(
                "[DRY RUN] Would convert %d IMS files across %d rounds "
                "(channel_order=%s, output=%s)",
                n_ims,
                len(round_folders),
                channel_order,
                merlin_dir,
            )
            return StageResult(
                status="skipped",
                metadata={"dry_run": True, "n_rounds": len(round_folders), "n_ims": n_ims},
            )

        # Ensure output directory exists
        merlin_dir.mkdir(parents=True, exist_ok=True)

        output_files: list[str] = []
        total_converted = 0
        total_skipped = 0
        conversion_errors: list[str] = []

        for round_num, round_dir in sorted(round_folders.items()):
            self.logger.info(
                "Processing round %d: %s", round_num, round_dir.name
            )

            # Find all IMS files and extract FOV numbers
            ims_entries: list[tuple[int, Path]] = []
            for ims_path in sorted(round_dir.iterdir()):
                if ims_path.suffix.lower() != ".ims":
                    continue
                fov = _get_fov_number(ims_path)
                if fov is None:
                    self.logger.debug("Skipping non-FOV file: %s", ims_path)
                    total_skipped += 1
                    continue
                ims_entries.append((fov, ims_path))

            if not ims_entries:
                self.logger.warning(
                    "No IMS files with FOV numbers found in round %d (%s)",
                    round_num,
                    round_dir,
                )
                continue

            ims_entries.sort(key=lambda x: x[0])

            self.logger.info(
                "  Found %d IMS files for round %d", len(ims_entries), round_num
            )

            # ----------------------------------------------------------
            # Convert IMS files (parallel within each round)
            # ----------------------------------------------------------
            def _convert_one(
                fov: int, ims_path: Path, rnd: int = round_num
            ) -> tuple[int, str | None, str | None]:
                """Worker function for parallel conversion.

                Returns ``(fov, output_path_str, error_message)``.
                """
                out_name = f"merFISH_merged_{rnd:02d}_{fov:03d}.tiff"
                out_path = merlin_dir / out_name
                try:
                    _convert_single_ims(ims_path, out_path, channel_order)
                    return fov, str(out_path), None
                except Exception as exc:
                    return fov, None, f"FOV {fov} round {rnd}: {exc}"

            effective_workers = min(max_workers, len(ims_entries))

            with ThreadPoolExecutor(max_workers=effective_workers) as pool:
                futures = {
                    pool.submit(_convert_one, fov, ims_path): fov
                    for fov, ims_path in ims_entries
                }
                for future in as_completed(futures):
                    fov_num, out_path_str, error = future.result()
                    if error:
                        self.logger.error("  Conversion failed: %s", error)
                        conversion_errors.append(error)
                    else:
                        output_files.append(out_path_str)
                        total_converted += 1
                        if total_converted % 10 == 0:
                            self.logger.info(
                                "  Converted %d files so far ...",
                                total_converted,
                            )

            # ----------------------------------------------------------
            # Generate stage positions CSV for this round
            # ----------------------------------------------------------
            positions_name = f"stagePos_Round#{round_num}.csv"
            positions_path = merlin_dir / positions_name

            try:
                _generate_positions_csv(ims_entries, positions_path)
                output_files.append(str(positions_path))
                self.logger.info(
                    "  Wrote positions file: %s (%d FOVs)",
                    positions_path,
                    len(ims_entries),
                )
            except Exception as exc:
                error_msg = (
                    f"Failed to generate positions CSV for round {round_num}: {exc}"
                )
                self.logger.error("  %s", error_msg)
                conversion_errors.append(error_msg)

        # ----------------------------------------------------------
        # Build result
        # ----------------------------------------------------------
        if conversion_errors and total_converted == 0:
            status = "failed"
            error_str = "; ".join(conversion_errors)
        elif conversion_errors:
            status = "completed"
            error_str = f"{len(conversion_errors)} error(s): " + "; ".join(
                conversion_errors[:5]
            )
        else:
            status = "completed"
            error_str = ""

        result = StageResult(
            status=status,
            output_files=output_files,
            metadata={
                "n_rounds": len(round_folders),
                "total_converted": total_converted,
                "total_skipped": total_skipped,
                "conversion_errors": len(conversion_errors),
                "channel_order": channel_order,
                "merlin_data_dir": str(merlin_dir),
            },
            error=error_str,
        )

        self.logger.info(
            "IMS conversion complete: %d files converted across %d rounds "
            "(%d errors, %d skipped)",
            total_converted,
            len(round_folders),
            len(conversion_errors),
            total_skipped,
        )

        self.write_run_metadata(
            result,
            start_time,
            parameters={
                "channel_order": channel_order,
                "max_workers": max_workers,
                "raw_data_dir": str(raw_dir),
                "merlin_data_dir": str(merlin_dir),
            },
        )

        return result
