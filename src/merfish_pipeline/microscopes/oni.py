"""ONI microscope adapter.

ONI raw data is organized as TIFF files in wavelength-specific directories::

    <raw_dir>/
        488nm, Raw/
            merFISH_01_000_00.TIFF
            merFISH_01_001_00.TIFF
            ...
        561nm, Raw/
            ...
        640nm, Raw/
            ...
        stagePos_Round#1.csv
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

from merfish_pipeline.io.path_utils import find_files_matching
from merfish_pipeline.io.sheet_io import read_sheet
from merfish_pipeline.microscopes.base import MicroscopeAdapter

# Regex that extracts wavelength, imaging round, FOV, and z-slice from the
# full path.  The wavelength is picked up from the directory name (e.g.
# ``488nm, Raw``), and the three numeric fields from the filename.
_PATH_RE = re.compile(
    r"(\d+)(?=nm)"           # wavelength (digits before "nm")
    r".*_(\d+)\D(\d+)\D(\d+)"  # ir, fov, z from filename
)


class ONIAdapter(MicroscopeAdapter):
    """Adapter for ONI Nanoimager raw data."""

    name = "oni"

    # ------------------------------------------------------------------
    # File discovery
    # ------------------------------------------------------------------

    def discover_raw_files(self, raw_dir: Path) -> list[dict[str, Any]]:
        raw_dir = Path(raw_dir)
        file_pattern = self.config.microscope.file_pattern  # e.g. "merFISH_{ir}_{fov}_{z}.TIFF"
        glob_pattern = file_pattern.replace("{ir}", "*").replace("{fov}", "*").replace("{z}", "*")

        records: list[dict[str, Any]] = []

        # Walk wavelength subdirectories
        for wv_dir in sorted(raw_dir.iterdir()):
            if not wv_dir.is_dir():
                continue
            # Must contain "nm" to be a wavelength folder
            if "nm" not in wv_dir.name.lower():
                continue

            matches = find_files_matching(wv_dir, glob_pattern)
            for fpath in matches:
                info = self.parse_filename(fpath)
                if info is None:
                    continue
                info["abs_path"] = str(fpath)
                info["file_size"] = fpath.stat().st_size
                records.append(info)

        self.logger.info(
            "%s discover: found %d files in %s", self.name.upper(), len(records), raw_dir
        )
        return records

    # ------------------------------------------------------------------
    # Filename parsing
    # ------------------------------------------------------------------

    def parse_filename(self, path: Path) -> dict[str, Any] | None:
        path = Path(path)
        full = str(path)
        m = _PATH_RE.search(full)
        if m is None:
            return None
        wavelength = int(m.group(1))
        ir = int(m.group(2))
        fov = int(m.group(3))
        z_slice = int(m.group(4))
        channel = f"{wavelength}nm, Raw"
        return {
            "round": ir,
            "fov": fov,
            "z_slice": z_slice,
            "channel": channel,
            "wavelength": wavelength,
        }

    # ------------------------------------------------------------------
    # Position reading
    # ------------------------------------------------------------------

    def read_positions(self, raw_dir: Path) -> pd.DataFrame:
        raw_dir = Path(raw_dir)
        mic = self.config.microscope
        pos_pattern = mic.position_file_pattern  # e.g. "stagePos_Round#{round}.csv"

        # Find all position files matching the pattern
        all_records: list[dict[str, Any]] = []
        glob_pat = pos_pattern.replace("{round}", "*")
        pos_files = find_files_matching(raw_dir, glob_pat)

        if not pos_files:
            self.logger.warning("No position files found in %s matching %s", raw_dir, glob_pat)
            return pd.DataFrame()

        for pf in sorted(pos_files):
            # Extract round number from filename
            round_num = self._extract_round_from_filename(pf.name, pos_pattern)
            df = read_sheet(pf)
            x_col = self.config.raw_data.stage_x_heading or mic.stage_x_heading
            y_col = self.config.raw_data.stage_y_heading or mic.stage_y_heading

            for idx, row in df.iterrows():
                record: dict[str, Any] = {
                    "round": round_num,
                    "tile_number": int(row["tile_number"]) if "tile_number" in df.columns else int(idx),
                    "stage_pos_x": float(row[x_col]),
                    "stage_pos_y": float(row[y_col]),
                }
                # Preserve grid coordinates if present in the position file
                for grid_col in ("grid_pos_x", "grid_pos_y"):
                    if grid_col in df.columns:
                        record[grid_col] = int(row[grid_col])
                # Include z-positions if present
                z_cols = [c for c in df.columns if c.lower().startswith("z")]
                for i, zc in enumerate(z_cols):
                    record[f"z_position_{i}"] = float(row[zc])
                all_records.append(record)

        result = pd.DataFrame(all_records)
        self.logger.info(
            "%s positions: %d rows from %d files", self.name.upper(), len(result), len(pos_files)
        )
        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_round_from_filename(filename: str, pattern: str) -> int:
        """Extract the round number from a position filename."""
        # Turn the pattern into a regex: stagePos_Round#{round}.csv -> stagePos_Round#(\d+).csv
        regex = re.escape(pattern).replace(r"\{round\}", r"(\d+)")
        # Also handle # character (common in ONI filenames)
        regex = regex.replace(r"\#", "#")
        m = re.search(regex, filename)
        if m:
            return int(m.group(1))
        # Fallback: find any number in the filename
        nums = re.findall(r"(\d+)", filename)
        return int(nums[-1]) if nums else 0
