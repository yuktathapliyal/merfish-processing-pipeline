"""Optimized TIFF I/O with parallel read support.

Uses :mod:`tifffile` for all low-level access.  Large files are written
with the BigTIFF format by default.  Multiple planes can be loaded in
parallel using a :class:`~concurrent.futures.ThreadPoolExecutor` to speed
up Z-stack assembly.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import tifffile

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Single-file operations
# ---------------------------------------------------------------------------


def read_tiff(path: str | Path) -> np.ndarray:
    """Read a TIFF file and return its contents as a NumPy array.

    Parameters
    ----------
    path:
        Path to the TIFF file.

    Returns
    -------
    np.ndarray
        Image data.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"TIFF file not found: {path}")
    logger.debug("Reading TIFF: %s", path)
    return tifffile.imread(str(path))


def read_tiff_shape(path: str | Path) -> tuple[int, ...]:
    """Read only the shape metadata of a TIFF file -- no pixel data is loaded.

    This is much faster than :func:`read_tiff` when you only need
    dimensions (e.g. to validate frame counts).

    Parameters
    ----------
    path:
        Path to the TIFF file.

    Returns
    -------
    tuple[int, ...]
        The array shape that ``tifffile.imread`` would return.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"TIFF file not found: {path}")
    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series
        if series:
            return tuple(series[0].shape)
        # Fallback: read page shapes and count pages.
        pages = tif.pages
        if not pages:
            return (0,)
        first_shape = pages[0].shape
        if len(pages) == 1:
            return first_shape
        return (len(pages), *first_shape)


def write_tiff(
    data: np.ndarray,
    path: str | Path,
    bigtiff: bool = True,
) -> None:
    """Write a NumPy array as a TIFF file.

    Parameters
    ----------
    data:
        Image data to write.
    path:
        Destination path.
    bigtiff:
        If *True* (default), the BigTIFF format is used, which supports
        files larger than 4 GiB.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    logger.debug("Writing TIFF (%s, bigtiff=%s): %s", data.shape, bigtiff, path)
    tifffile.imwrite(str(path), data, bigtiff=bigtiff)


# ---------------------------------------------------------------------------
# Parallel / multi-file operations
# ---------------------------------------------------------------------------


def read_planes_parallel(
    paths: list[Path | str],
    workers: int = 4,
) -> list[np.ndarray]:
    """Read multiple TIFF files in parallel.

    Parameters
    ----------
    paths:
        Ordered list of TIFF file paths.
    workers:
        Number of threads for the :class:`~concurrent.futures.ThreadPoolExecutor`.

    Returns
    -------
    list[np.ndarray]
        Images in the **same order** as *paths*.
    """
    n = len(paths)
    if n == 0:
        return []

    results: list[np.ndarray | None] = [None] * n

    def _load(idx: int, p: Path | str) -> tuple[int, np.ndarray]:
        return idx, read_tiff(p)

    effective_workers = min(workers, n)
    logger.info("Reading %d TIFF planes with %d workers", n, effective_workers)

    with ThreadPoolExecutor(max_workers=effective_workers) as pool:
        futures = {
            pool.submit(_load, i, p): i for i, p in enumerate(paths)
        }
        for future in as_completed(futures):
            idx, data = future.result()
            results[idx] = data

    return results  # type: ignore[return-value]


def merge_stack(
    plane_paths: list[Path | str],
    output_path: str | Path,
    dtype: np.dtype | type = np.uint16,
    workers: int = 4,
) -> None:
    """Read plane images in parallel and write a merged TIFF stack.

    The individual plane images are loaded concurrently, cast to *dtype*,
    stacked along a new leading axis, and written as a single TIFF.

    Parameters
    ----------
    plane_paths:
        Ordered list of paths to individual plane images.
    output_path:
        Where to write the merged stack.
    dtype:
        Output data type (default ``np.uint16``).
    workers:
        Number of parallel read threads.

    Raises
    ------
    ValueError
        If *plane_paths* is empty.
    """
    if not plane_paths:
        raise ValueError("plane_paths must not be empty")

    planes = read_planes_parallel(plane_paths, workers=workers)

    logger.info(
        "Merging %d planes -> %s (dtype=%s)",
        len(planes),
        output_path,
        np.dtype(dtype).name,
    )
    stack = np.stack([p.astype(dtype) for p in planes], axis=0)
    write_tiff(stack, output_path, bigtiff=True)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_tiff_stack(path: str | Path, expected_frames: int) -> bool:
    """Check that a TIFF file has the expected number of frames.

    Only metadata is read -- no pixel data is loaded, making this a very
    fast sanity check.

    Parameters
    ----------
    path:
        Path to the TIFF file.
    expected_frames:
        The number of frames (leading dimension) expected.

    Returns
    -------
    bool
        *True* if the first dimension of the shape equals *expected_frames*.
    """
    path = Path(path)
    if not path.exists():
        logger.warning("TIFF file does not exist: %s", path)
        return False

    shape = read_tiff_shape(path)
    if len(shape) < 2:
        logger.warning(
            "TIFF %s has unexpected shape %s (expected at least 2 dimensions)",
            path,
            shape,
        )
        return False

    actual_frames = shape[0]
    if actual_frames != expected_frames:
        logger.warning(
            "TIFF %s has %d frames, expected %d",
            path,
            actual_frames,
            expected_frames,
        )
        return False

    logger.debug("TIFF %s validated: %d frames", path, actual_frames)
    return True
