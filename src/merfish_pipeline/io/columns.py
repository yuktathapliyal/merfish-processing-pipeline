"""Shared column-name candidate lists and detection helper for barcode tables.

MERlin and the various pipeline stages disagree on the exact spelling of
common columns (``fov`` vs ``FOV``, ``z`` vs ``zIndex``, ``x`` vs
``global_x``).  Several stages used to maintain identical lists of
fallback names plus a private ``_detect_column`` helper; those have been
consolidated here.

Two distinct X / Y candidate lists exist on purpose:

- ``LOCAL_X_CANDIDATES`` / ``LOCAL_Y_CANDIDATES`` -- per-FOV pixel
  coordinates used by ``cell_assignment``, ``barcode_qc``, etc.
- ``GLOBAL_X_CANDIDATES`` / ``GLOBAL_Y_CANDIDATES`` -- stitched
  global-stage coordinates preferred by ``spatial_visualization``;
  ``global_x`` ranks first so that whole-experiment plots use stitched
  coordinates when available.
"""

from __future__ import annotations

import pandas as pd

#: Names commonly used for the FOV / tile-number column.
FOV_CANDIDATES: list[str] = ["fov", "FOV", "Fov"]

#: Per-FOV pixel coordinate column names (preferred by stages that work
#: against single-FOV mask / image data).
LOCAL_X_CANDIDATES: list[str] = ["x", "X"]
LOCAL_Y_CANDIDATES: list[str] = ["y", "Y"]

#: Global / stitched coordinate column names (preferred by visualisation
#: stages).  ``global_x`` is first so the stitched value wins when present.
GLOBAL_X_CANDIDATES: list[str] = ["global_x", "x", "X"]
GLOBAL_Y_CANDIDATES: list[str] = ["global_y", "y", "Y"]

#: Z-slice column names.
Z_CANDIDATES: list[str] = ["z", "Z", "zIndex", "z_index", "zPos", "zpos"]


def detect_column(
    df: pd.DataFrame, candidates: list[str], label: str
) -> str:
    """Return the first matching column name from *candidates*.

    Parameters
    ----------
    df:
        DataFrame to inspect.
    candidates:
        Ordered list of column names to try.
    label:
        Human-readable name used in the error message when no match is
        found (e.g. ``"FOV"`` or ``"x"``).

    Raises
    ------
    ValueError
        When none of the *candidates* are present in ``df.columns``.
    """
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(
        f"Cannot auto-detect {label} column. "
        f"Tried {candidates}; available: {list(df.columns)}"
    )
