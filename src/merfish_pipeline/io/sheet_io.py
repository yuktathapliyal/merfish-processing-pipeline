"""Unified CSV / Excel / TSV reader and writer.

Replaces the duplicated ``read_sheet()`` helper that was copy-pasted across
multiple legacy scripts.  Supports *.csv*, *.tsv*, *.xlsx*, *.xls*, and
*.csv.gz* files with automatic delimiter detection for CSV variants.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_EXCEL_EXTENSIONS = {".xlsx", ".xls"}
_CSV_EXTENSIONS = {".csv", ".tsv", ".gz"}  # .gz handled via suffixes check


def _suffixes_key(path: Path) -> str:
    """Return a normalised extension key.

    For compound extensions like ``.csv.gz`` the key is ``".csv.gz"``;
    otherwise it is the single suffix (lower-cased).
    """
    suffixes = [s.lower() for s in path.suffixes]
    if len(suffixes) >= 2 and suffixes[-1] == ".gz":
        return "".join(suffixes[-2:])
    return suffixes[-1] if suffixes else ""


def _sniff_delimiter(path: Path, encoding: str = "utf-8") -> str:
    """Peek at the first few lines to guess the CSV delimiter."""
    import gzip

    opener = gzip.open if _suffixes_key(path) == ".csv.gz" else open

    try:
        with opener(path, "rt", encoding=encoding, errors="replace") as fh:  # type: ignore[arg-type]
            sample = fh.read(8192)
    except Exception:
        return ","

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        return str(dialect.delimiter)
    except csv.Error:
        return ","


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def read_sheet(path: str | Path, **kwargs: object) -> pd.DataFrame:
    """Read a tabular file into a :class:`~pandas.DataFrame`.

    The file format is inferred from the extension:

    * **.csv** / **.csv.gz** -- comma-separated (delimiter auto-detected)
    * **.tsv** -- tab-separated
    * **.xlsx** / **.xls** -- Excel workbook (first sheet by default; pass
      *sheet_name* to select a different one)

    All extra *kwargs* are forwarded to the underlying pandas reader
    (``read_csv`` or ``read_excel``).

    Parameters
    ----------
    path:
        Filesystem path to the input file.
    **kwargs:
        Passed through to ``pd.read_csv`` or ``pd.read_excel``.

    Returns
    -------
    pd.DataFrame

    Raises
    ------
    ValueError
        If the extension is not recognised.
    FileNotFoundError
        If *path* does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    ext = _suffixes_key(path)

    if ext in _EXCEL_EXTENSIONS:
        logger.debug("Reading Excel file: %s", path)
        return pd.read_excel(path, **kwargs)  # type: ignore[arg-type]

    if ext in {".csv", ".csv.gz", ".tsv"}:
        # Let the caller override the separator; otherwise auto-detect.
        if "sep" not in kwargs and "delimiter" not in kwargs:
            if ext == ".tsv":
                kwargs["sep"] = "\t"
            else:
                kwargs["sep"] = _sniff_delimiter(path)
        logger.debug("Reading delimited file (%s): %s", kwargs.get("sep", "?"), path)
        return pd.read_csv(path, **kwargs)  # type: ignore[arg-type]

    raise ValueError(
        f"Unsupported file extension '{ext}' for: {path}. "
        "Expected one of: .csv, .tsv, .xlsx, .xls, .csv.gz"
    )


def write_sheet(df: pd.DataFrame, path: str | Path, **kwargs: object) -> None:
    """Write a :class:`~pandas.DataFrame` to a CSV or Excel file.

    The output format is inferred from the extension of *path*.

    * **.csv** / **.tsv** -- uses ``DataFrame.to_csv``
    * **.xlsx** -- uses ``DataFrame.to_excel`` (requires *openpyxl*)

    Extra *kwargs* are forwarded to the underlying writer.

    Parameters
    ----------
    df:
        The DataFrame to persist.
    path:
        Destination path.
    **kwargs:
        Forwarded to ``to_csv`` / ``to_excel``.

    Raises
    ------
    ValueError
        If the extension is not recognised.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    ext = _suffixes_key(path)

    # Provide sensible defaults but allow caller overrides.
    if ext == ".csv" or ext == ".csv.gz":
        kwargs.setdefault("index", False)  # type: ignore[arg-type]
        logger.debug("Writing CSV: %s", path)
        df.to_csv(path, **kwargs)  # type: ignore[arg-type]
    elif ext == ".tsv":
        kwargs.setdefault("index", False)  # type: ignore[arg-type]
        kwargs.setdefault("sep", "\t")  # type: ignore[arg-type]
        logger.debug("Writing TSV: %s", path)
        df.to_csv(path, **kwargs)  # type: ignore[arg-type]
    elif ext in _EXCEL_EXTENSIONS:
        kwargs.setdefault("index", False)  # type: ignore[arg-type]
        logger.debug("Writing Excel: %s", path)
        df.to_excel(path, **kwargs)  # type: ignore[arg-type]
    else:
        raise ValueError(
            f"Unsupported file extension '{ext}' for: {path}. "
            "Expected one of: .csv, .tsv, .xlsx"
        )
