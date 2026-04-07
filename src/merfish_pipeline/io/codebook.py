"""Codebook loading and blank-barcode helpers shared across stages.

Several post-MERlin stages (``barcode_qc``, ``anndata_export``,
``spatial_visualization``, ``correlation``) need to:

1. Load and normalise a MERlin codebook CSV (column renaming, default IDs).
2. Identify "blank" / control barcodes by name.

Both operations were duplicated as private helpers in three stage files
before being consolidated here.  Each consumer should now import from this
module instead of redefining the logic locally.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from merfish_pipeline.io.sheet_io import read_sheet

#: Regex matching blank / control barcode names such as ``Blank_01``,
#: ``blank-3``, ``Blank02``, etc.  Used to flag rows in the codebook /
#: barcodes table that should be excluded from gene-level analyses.
BLANK_RE = re.compile(r"^[Bb]lank[-_]?\d+$")


def is_blank(gene_symbol: str) -> bool:
    """Return ``True`` when *gene_symbol* matches the blank-barcode pattern."""
    return bool(BLANK_RE.match(str(gene_symbol)))


def load_codebook(codebook_path: Path) -> pd.DataFrame:
    """Load a MERlin codebook and normalise its column names.

    Guarantees the returned DataFrame has at least the columns
    ``barcode_id`` and ``gene_symbol``:

    - ``barcode_id`` defaults to the row index when missing.
    - ``gene_symbol`` is renamed from ``name`` or ``gene_name`` when those
      legacy column names are present.

    Other columns (e.g. per-bit codeword columns) are passed through
    unchanged.
    """
    cb = read_sheet(codebook_path)
    if "barcode_id" not in cb.columns:
        cb["barcode_id"] = cb.index
    if "gene_symbol" not in cb.columns:
        if "name" in cb.columns:
            cb = cb.rename(columns={"name": "gene_symbol"})
        elif "gene_name" in cb.columns:
            cb = cb.rename(columns={"gene_name": "gene_symbol"})
    return cb
