"""Abstract base class for microscope adapters.

Each microscope type (ONI, NIKON, ANDOR) subclasses ``MicroscopeAdapter`` and
provides concrete implementations for file discovery, filename parsing, and
position reading.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


class MicroscopeAdapter(ABC):
    """Abstract interface that every microscope adapter must implement.

    Parameters
    ----------
    config:
        Fully resolved ``PipelineConfig`` instance.
    """

    name: str = ""

    def __init__(self, config: Any):
        self.config = config
        self.logger = logging.getLogger(f"merfish_pipeline.microscopes.{self.name}")

    # ------------------------------------------------------------------
    # File discovery
    # ------------------------------------------------------------------

    @abstractmethod
    def discover_raw_files(self, raw_dir: Path) -> list[dict[str, Any]]:
        """Scan *raw_dir* and return a list of file-info dicts.

        Each dict must have at least:
        ``round``, ``fov``, ``z_slice``, ``channel``, ``wavelength``,
        ``abs_path``, ``file_size``.
        """

    # ------------------------------------------------------------------
    # Filename parsing
    # ------------------------------------------------------------------

    @abstractmethod
    def parse_filename(self, path: Path) -> dict[str, Any] | None:
        """Parse a raw-data filename and return extracted metadata.

        Returns ``None`` if the filename does not match the expected pattern.
        """

    # ------------------------------------------------------------------
    # Position reading
    # ------------------------------------------------------------------

    @abstractmethod
    def read_positions(self, raw_dir: Path) -> pd.DataFrame:
        """Read and normalize position data into a standard DataFrame.

        The returned DataFrame must have columns:
        ``round``, ``tile_number``, ``stage_pos_x``, ``stage_pos_y``,
        and optionally ``z_position_0 .. z_position_N``.
        """

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def get_bead_dir(self, raw_dir: Path) -> Path:
        """Return the bead-channel image directory."""
        bead_folder = self.config.raw_data.bead_channel_folder
        return raw_dir / bead_folder

    def get_image_shape(self, path: Path) -> tuple[int, ...]:
        """Read the shape of a single image without loading pixel data."""
        from merfish_pipeline.io.tiff_io import read_tiff_shape

        return read_tiff_shape(path)
