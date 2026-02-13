"""Cross-platform path utilities.

On WSL (Windows Subsystem for Linux) the project frequently needs to translate
between Windows-style paths (``D:\\Data\\experiment``) and POSIX mount points
(``/mnt/d/Data/experiment``).  The helpers in this module handle that
conversion transparently and also provide a few convenience wrappers around
common filesystem operations.
"""

from __future__ import annotations

import logging
import platform
import re
from pathlib import Path, PurePosixPath, PureWindowsPath

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Windows <-> WSL path conversion
# ---------------------------------------------------------------------------

# Matches a Windows absolute path, e.g. "D:\foo\bar" or "D:/foo/bar"
_WIN_ABS_RE = re.compile(r"^([A-Za-z]):[\\\/]")


def _is_wsl() -> bool:
    """Return *True* when running inside WSL."""
    try:
        return "microsoft" in platform.uname().release.lower()
    except Exception:
        return False


def windows_to_wsl_path(path: str) -> Path:
    """Convert a Windows path to its WSL ``/mnt/`` equivalent.

    Examples
    --------
    >>> windows_to_wsl_path(r"D:\\Data\\experiment")
    PosixPath('/mnt/d/Data/experiment')

    Parameters
    ----------
    path:
        A Windows-style absolute path (e.g. ``D:\\foo\\bar``).

    Returns
    -------
    Path
        The corresponding WSL POSIX path.

    Raises
    ------
    ValueError
        If *path* is not a recognised Windows absolute path.
    """
    m = _WIN_ABS_RE.match(path)
    if not m:
        raise ValueError(f"Not a Windows absolute path: {path!r}")
    drive = m.group(1).lower()
    rest = path[m.end():].replace("\\", "/")
    return Path(f"/mnt/{drive}/{rest}")


def wsl_to_windows_path(path: str | Path) -> str:
    """Convert a WSL ``/mnt/<drive>/...`` path to a Windows path.

    Examples
    --------
    >>> wsl_to_windows_path(Path("/mnt/d/Data/experiment"))
    'D:\\\\Data\\\\experiment'

    Parameters
    ----------
    path:
        A WSL POSIX path rooted at ``/mnt/<drive_letter>/``.

    Returns
    -------
    str
        The equivalent Windows path string.

    Raises
    ------
    ValueError
        If *path* does not match the ``/mnt/<letter>/`` pattern.
    """
    posix = PurePosixPath(path)
    parts = posix.parts  # ('/', 'mnt', 'd', 'Data', ...)
    if len(parts) < 3 or parts[0] != "/" or parts[1] != "mnt" or len(parts[2]) != 1:
        raise ValueError(f"Not a WSL mount path: {path}")
    drive = parts[2].upper()
    rest = PureWindowsPath(*parts[3:]) if len(parts) > 3 else PureWindowsPath("")
    return f"{drive}:\\{rest}"


def normalize_path(path: str | Path) -> Path:
    """Normalize *path* for the current platform.

    If running on WSL and *path* looks like a Windows absolute path it is
    automatically converted.  Otherwise, the path is resolved and returned
    as a :class:`~pathlib.Path`.

    Parameters
    ----------
    path:
        A filesystem path in either Windows or POSIX format.

    Returns
    -------
    Path
    """
    s = str(path)

    # Detect Windows-style absolute paths and convert when on WSL.
    if _WIN_ABS_RE.match(s):
        if _is_wsl():
            logger.debug("Auto-converting Windows path to WSL: %s", s)
            return windows_to_wsl_path(s)
        # On native Windows, just resolve directly.
        return Path(s).resolve()

    return Path(s).resolve()


# ---------------------------------------------------------------------------
# Filesystem helpers
# ---------------------------------------------------------------------------


def ensure_dir(path: str | Path) -> Path:
    """Create *path* (and parents) if it does not exist, then return it.

    Parameters
    ----------
    path:
        Directory path to create.

    Returns
    -------
    Path
        The same *path* as a resolved :class:`~pathlib.Path`.
    """
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def find_files_matching(directory: str | Path, pattern: str) -> list[Path]:
    """Find files that match a glob *pattern* under *directory*.

    Results are sorted by filename to guarantee deterministic ordering
    regardless of filesystem enumeration order.

    Parameters
    ----------
    directory:
        Root directory to search.
    pattern:
        A glob pattern (e.g. ``"*.tif"`` or ``"**/*.csv"``).

    Returns
    -------
    list[Path]
        Sorted list of matching paths.

    Raises
    ------
    FileNotFoundError
        If *directory* does not exist.
    """
    d = Path(directory)
    if not d.is_dir():
        raise FileNotFoundError(f"Directory not found: {d}")
    return sorted(d.glob(pattern), key=lambda p: p.name)
