"""I/O utilities for the merFISH processing pipeline.

Submodules
----------
sheet_io
    Unified CSV / Excel / TSV reader and writer.
path_utils
    Cross-platform (Windows / WSL) path helpers.
tiff_io
    Optimized TIFF I/O with parallel read support.
hdf5_io
    HDF5 / IMS reader for ANDOR microscope files.
codebook
    Codebook loading and blank-barcode helpers shared across stages.
columns
    Shared column-name candidate lists and detection helper for barcode
    tables.
"""

from merfish_pipeline.io.codebook import BLANK_RE, is_blank, load_codebook
from merfish_pipeline.io.columns import (
    FOV_CANDIDATES,
    GLOBAL_X_CANDIDATES,
    GLOBAL_Y_CANDIDATES,
    LOCAL_X_CANDIDATES,
    LOCAL_Y_CANDIDATES,
    Z_CANDIDATES,
    detect_column,
)
from merfish_pipeline.io.hdf5_io import (
    extract_stage_positions,
    list_ims_contents,
    read_ims_channel,
    read_ims_metadata,
    read_ims_stack,
)
from merfish_pipeline.io.path_utils import (
    ensure_dir,
    find_files_matching,
    normalize_path,
    windows_to_wsl_path,
    wsl_to_windows_path,
)
from merfish_pipeline.io.sheet_io import read_sheet, write_sheet
from merfish_pipeline.io.tiff_io import (
    merge_stack,
    read_planes_parallel,
    read_tiff,
    read_tiff_shape,
    validate_tiff_stack,
    write_tiff,
)

__all__ = [
    # sheet_io
    "read_sheet",
    "write_sheet",
    # path_utils
    "windows_to_wsl_path",
    "wsl_to_windows_path",
    "normalize_path",
    "ensure_dir",
    "find_files_matching",
    # tiff_io
    "read_tiff",
    "read_tiff_shape",
    "write_tiff",
    "read_planes_parallel",
    "merge_stack",
    "validate_tiff_stack",
    # hdf5_io
    "read_ims_metadata",
    "read_ims_channel",
    "read_ims_stack",
    "extract_stage_positions",
    "list_ims_contents",
    # codebook
    "BLANK_RE",
    "is_blank",
    "load_codebook",
    # columns
    "FOV_CANDIDATES",
    "LOCAL_X_CANDIDATES",
    "LOCAL_Y_CANDIDATES",
    "GLOBAL_X_CANDIDATES",
    "GLOBAL_Y_CANDIDATES",
    "Z_CANDIDATES",
    "detect_column",
]
