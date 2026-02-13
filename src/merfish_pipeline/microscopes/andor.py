"""ANDOR (Imaris IMS) microscope adapter.

ANDOR raw data is organized as HDF5 ``.ims`` files in round-specific
directories::

    <raw_dir>/
        1st round/
            test_2025-10-23_test1_F000.ims
            test_2025-10-23_test1_F001.ims
            ...
        2nd round/
            ...
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd

from merfish_pipeline.io.hdf5_io import extract_stage_positions, read_ims_metadata
from merfish_pipeline.microscopes.base import MicroscopeAdapter

# Regex patterns for detecting round numbers from folder names.
_ROUND_PATTERNS = [
    re.compile(r"(\d+)(?:st|nd|rd|th)\s*round", re.IGNORECASE),
    re.compile(r"round\s*(\d+)", re.IGNORECASE),
    re.compile(r"R(\d+)", re.IGNORECASE),
]

# Regex for extracting FOV number from IMS filenames.
_FOV_RE = re.compile(r"_F(\d+)\.ims$", re.IGNORECASE)


class ANDORAdapter(MicroscopeAdapter):
    """Adapter for ANDOR Imaris (IMS/HDF5) raw data."""

    name = "andor"

    # ------------------------------------------------------------------
    # File discovery
    # ------------------------------------------------------------------

    def discover_raw_files(self, raw_dir: Path) -> list[dict[str, Any]]:
        raw_dir = Path(raw_dir)
        round_folders = self._find_round_folders(raw_dir)

        if not round_folders:
            self.logger.warning("No round folders found in %s", raw_dir)
            return []

        records: list[dict[str, Any]] = []

        for round_num, round_dir in sorted(round_folders.items()):
            ims_files = sorted(
                p for p in round_dir.iterdir()
                if p.suffix.lower() == ".ims"
            )
            for ims_path in ims_files:
                fov = self._get_fov_number(ims_path)
                if fov is None:
                    self.logger.debug("Skipping non-FOV file: %s", ims_path)
                    continue

                # Read IMS metadata for channel and z-slice info
                try:
                    meta = read_ims_metadata(ims_path)
                except Exception as exc:
                    self.logger.warning("Failed to read metadata from %s: %s", ims_path, exc)
                    continue

                n_channels = meta["n_channels"]
                n_z = meta["n_z_slices"]
                channel_names = meta.get("channel_names", [])

                for ch in range(n_channels):
                    ch_name = channel_names[ch] if ch < len(channel_names) else f"Channel {ch}"
                    for z in range(n_z):
                        records.append({
                            "round": round_num,
                            "fov": fov,
                            "z_slice": z,
                            "channel": ch_name,
                            "wavelength": ch,  # Channel index for IMS
                            "abs_path": str(ims_path),
                            "file_size": ims_path.stat().st_size,
                            "image_shape": meta["image_shape"],
                        })

        self.logger.info(
            "ANDOR discover: found %d file entries across %d rounds in %s",
            len(records), len(round_folders), raw_dir,
        )
        return records

    # ------------------------------------------------------------------
    # Filename parsing
    # ------------------------------------------------------------------

    def parse_filename(self, path: Path) -> dict[str, Any] | None:
        path = Path(path)
        if path.suffix.lower() != ".ims":
            return None

        fov = self._get_fov_number(path)
        if fov is None:
            return None

        # Try to get round number from parent directory
        round_num = self._get_round_number(path.parent.name)
        if round_num is None:
            round_num = 0

        return {
            "round": round_num,
            "fov": fov,
            "z_slice": 0,  # IMS files contain all z-slices
            "channel": "all",
            "wavelength": 0,
        }

    # ------------------------------------------------------------------
    # Position reading
    # ------------------------------------------------------------------

    def read_positions(self, raw_dir: Path) -> pd.DataFrame:
        raw_dir = Path(raw_dir)
        round_folders = self._find_round_folders(raw_dir)

        if not round_folders:
            self.logger.warning("No round folders found in %s", raw_dir)
            return pd.DataFrame()

        all_records: list[dict[str, Any]] = []

        for round_num, round_dir in sorted(round_folders.items()):
            ims_files = sorted(
                p for p in round_dir.iterdir()
                if p.suffix.lower() == ".ims"
            )
            for ims_path in ims_files:
                fov = self._get_fov_number(ims_path)
                if fov is None:
                    continue

                try:
                    positions = extract_stage_positions(ims_path)
                except Exception as exc:
                    self.logger.warning(
                        "Failed to extract positions from %s: %s", ims_path, exc
                    )
                    continue

                record: dict[str, Any] = {
                    "round": round_num,
                    "tile_number": fov,
                    "stage_pos_x": positions["stage_pos_x"],
                    "stage_pos_y": positions["stage_pos_y"],
                }
                # Add per-z positions
                for i, zpos in enumerate(positions.get("z_positions", [])):
                    record[f"z_position_{i}"] = zpos

                all_records.append(record)

        result = pd.DataFrame(all_records)
        self.logger.info(
            "ANDOR positions: %d rows from %d rounds", len(result), len(round_folders)
        )
        return result

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _find_round_folders(self, raw_dir: Path) -> dict[int, Path]:
        """Return ``{round_number: Path}`` for all detected round folders."""
        result: dict[int, Path] = {}
        for entry in sorted(raw_dir.iterdir()):
            if not entry.is_dir():
                continue
            round_num = self._get_round_number(entry.name)
            if round_num is not None:
                result[round_num] = entry
        return result

    @staticmethod
    def _get_round_number(folder_name: str) -> int | None:
        """Extract the round number from a folder name."""
        for pattern in _ROUND_PATTERNS:
            m = pattern.search(folder_name)
            if m:
                return int(m.group(1))
        return None

    @staticmethod
    def _get_fov_number(path: Path) -> int | None:
        """Extract the FOV number from an IMS filename."""
        m = _FOV_RE.search(path.name)
        if m:
            return int(m.group(1))
        return None

    def get_bead_dir(self, raw_dir: Path) -> Path:
        """ANDOR does not have a separate bead directory.

        Returns the first round folder as a fallback.
        """
        round_folders = self._find_round_folders(raw_dir)
        if round_folders:
            return next(iter(round_folders.values()))
        return raw_dir
