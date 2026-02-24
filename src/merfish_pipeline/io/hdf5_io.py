"""HDF5 / IMS reader for ANDOR microscope files.

Imaris **.ims** files are HDF5 containers that follow a specific hierarchy.
The key paths this module understands are:

* ``/DataSet/ResolutionLevel 0/TimePoint 0/Channel {n}/Data`` -- image data
* ``/DataSetInfo/Channel {n}`` -- channel metadata (name, colour, range)
* ``/DataSetInfo/Image`` -- global image dimensions and physical extents

All public functions accept a *path* argument pointing at an ``.ims`` file
and handle the ``h5py.File`` context internally so callers never need to
manage HDF5 handles directly.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import h5py
import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_DS_ROOT = "DataSet/ResolutionLevel 0/TimePoint 0"
_INFO_ROOT = "DataSetInfo"


def _open(path: Path | str) -> h5py.File:
    """Open an IMS file for reading, with a clear error on failure."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"IMS file not found: {path}")
    return h5py.File(str(path), "r")


def _attr_str(group: h5py.Group | h5py.Dataset, key: str, default: str = "") -> str:
    """Read an HDF5 attribute and decode to *str* if it is bytes."""
    val = group.attrs.get(key, default)
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="replace")
    if isinstance(val, np.ndarray):
        # Some IMS writers store single-element byte arrays.
        try:
            return val.flat[0].decode("utf-8", errors="replace")
        except (AttributeError, UnicodeDecodeError):
            return str(val.flat[0])
    return str(val) if val is not default else default


def _attr_float(group: h5py.Group | h5py.Dataset, key: str, default: float = 0.0) -> float:
    """Read an HDF5 attribute and convert to *float*."""
    raw = _attr_str(group, key, "")
    if not raw:
        return default
    try:
        return float(raw)
    except (ValueError, TypeError):
        return default


def _count_channels(f: h5py.File) -> int:
    """Return the number of channels present in the dataset."""
    tp_group = f.get(_DS_ROOT)
    if tp_group is None:
        return 0
    count = 0
    while f"Channel {count}" in tp_group:
        count += 1
    return count


def _count_z_slices(f: h5py.File, channel: int = 0) -> int:
    """Return the number of Z slices for a given channel."""
    ds_path = f"{_DS_ROOT}/Channel {channel}/Data"
    ds = f.get(ds_path)
    if ds is None:
        return 0
    # Data shape is (Z, H, W) for a single-channel stack.
    return int(ds.shape[0]) if ds.ndim >= 3 else 1


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def read_ims_metadata(path: str | Path) -> dict[str, Any]:
    """Read metadata from an IMS (Imaris HDF5) file.

    Parameters
    ----------
    path:
        Path to the ``.ims`` file.

    Returns
    -------
    dict
        Keys:

        * **n_channels** (*int*) -- number of fluorescence channels
        * **n_z_slices** (*int*) -- number of Z slices (from channel 0)
        * **image_shape** (*tuple[int, int]*) -- *(height, width)*
        * **channel_names** (*list[str]*) -- human-readable channel labels
        * **voxel_size** (*tuple[float, float, float]*) -- *(x, y, z)* in
          microns
        * **stage_position** (*dict*) -- ``{x, y, z_min, z_max, z_step}``
    """
    path = Path(path)
    logger.debug("Reading IMS metadata: %s", path)

    with _open(path) as f:
        n_channels = _count_channels(f)
        n_z_slices = _count_z_slices(f, channel=0) if n_channels > 0 else 0

        # Prefer ImageSizeZ metadata over array shape: IMS arrays may be
        # padded beyond the actual number of acquired z-slices (e.g. array
        # shape = 48 but only 41 slices were acquired).
        if n_channels > 0:
            ch0_group = f.get(f"{_DS_ROOT}/Channel 0")
            if ch0_group is not None:
                image_size_z = _attr_str(ch0_group, "ImageSizeZ", "")
                if image_size_z:
                    try:
                        n_z_slices = int(image_size_z)
                    except (ValueError, TypeError):
                        pass  # keep array-based count

        # --- image shape (height, width) from actual data -----------------
        image_shape: tuple[int, int] = (0, 0)
        ds = f.get(f"{_DS_ROOT}/Channel 0/Data")
        if ds is not None:
            if ds.ndim >= 3:
                image_shape = (int(ds.shape[1]), int(ds.shape[2]))
            elif ds.ndim == 2:
                image_shape = (int(ds.shape[0]), int(ds.shape[1]))

        # --- channel names ------------------------------------------------
        channel_names: list[str] = []
        for ch in range(n_channels):
            info_grp = f.get(f"{_INFO_ROOT}/Channel {ch}")
            if info_grp is not None:
                name = _attr_str(info_grp, "Name", f"Channel {ch}")
                channel_names.append(name)
            else:
                channel_names.append(f"Channel {ch}")

        # --- voxel size / extents from DataSetInfo/Image ------------------
        img_info = f.get(f"{_INFO_ROOT}/Image")
        voxel_size = (0.0, 0.0, 0.0)
        stage_position: dict[str, Any] = {}

        if img_info is not None:
            # Image dimensions stored as attributes
            x_size = _attr_float(img_info, "X", 0)
            y_size = _attr_float(img_info, "Y", 0)
            z_size = _attr_float(img_info, "Z", 0)

            ext_min_x = _attr_float(img_info, "ExtMin0", 0)
            ext_min_y = _attr_float(img_info, "ExtMin1", 0)
            ext_min_z = _attr_float(img_info, "ExtMin2", 0)
            ext_max_x = _attr_float(img_info, "ExtMax0", 0)
            ext_max_y = _attr_float(img_info, "ExtMax1", 0)
            ext_max_z = _attr_float(img_info, "ExtMax2", 0)

            # Voxel size = extent_range / pixel_count
            vx = (ext_max_x - ext_min_x) / x_size if x_size > 0 else 0.0
            vy = (ext_max_y - ext_min_y) / y_size if y_size > 0 else 0.0
            vz = (ext_max_z - ext_min_z) / z_size if z_size > 0 else 0.0
            voxel_size = (vx, vy, vz)

            z_step = (ext_max_z - ext_min_z) / (n_z_slices - 1) if n_z_slices > 1 else 0.0

            stage_position = {
                "x": ext_min_x,
                "y": ext_min_y,
                "z_min": ext_min_z,
                "z_max": ext_max_z,
                "z_step": z_step,
            }
        else:
            logger.warning("No DataSetInfo/Image group found in %s", path)

    return {
        "n_channels": n_channels,
        "n_z_slices": n_z_slices,
        "image_shape": image_shape,
        "channel_names": channel_names,
        "voxel_size": voxel_size,
        "stage_position": stage_position,
    }


def read_ims_channel(path: str | Path, channel: int, z_slice: int) -> np.ndarray:
    """Read a single channel / Z-slice from an IMS file.

    Parameters
    ----------
    path:
        Path to the ``.ims`` file.
    channel:
        Zero-based channel index.
    z_slice:
        Zero-based Z-slice index.

    Returns
    -------
    np.ndarray
        2-D image array *(H, W)*.
    """
    path = Path(path)
    ds_path = f"{_DS_ROOT}/Channel {channel}/Data"

    with _open(path) as f:
        ds = f.get(ds_path)
        if ds is None:
            raise KeyError(
                f"Channel {channel} not found in {path}. "
                f"Available channels: {_count_channels(f)}"
            )
        n_z = int(ds.shape[0]) if ds.ndim >= 3 else 1
        if z_slice < 0 or z_slice >= n_z:
            raise IndexError(
                f"z_slice={z_slice} out of range [0, {n_z}) for channel {channel}"
            )
        if ds.ndim >= 3:
            data: np.ndarray = ds[z_slice, :, :]
        else:
            data = ds[:]
    return data


def read_ims_stack(path: str | Path, channel: int) -> np.ndarray:
    """Read all Z-slices for a single channel from an IMS file.

    Parameters
    ----------
    path:
        Path to the ``.ims`` file.
    channel:
        Zero-based channel index.

    Returns
    -------
    np.ndarray
        3-D array with shape *(Z, H, W)*.
    """
    path = Path(path)
    ds_path = f"{_DS_ROOT}/Channel {channel}/Data"
    logger.debug("Reading IMS stack: %s channel=%d", path, channel)

    with _open(path) as f:
        ds = f.get(ds_path)
        if ds is None:
            raise KeyError(
                f"Channel {channel} not found in {path}. "
                f"Available channels: {_count_channels(f)}"
            )
        data: np.ndarray = ds[:]
    # Ensure 3-D even if there is only one Z-slice.
    if data.ndim == 2:
        data = data[np.newaxis, :, :]
    return data


def extract_stage_positions(path: str | Path) -> dict[str, Any]:
    """Extract stage X, Y, Z positions from IMS metadata.

    Parameters
    ----------
    path:
        Path to the ``.ims`` file.

    Returns
    -------
    dict
        Keys:

        * **stage_pos_x** (*float*)
        * **stage_pos_y** (*float*)
        * **z_positions** (*list[float]*) -- one value per Z-slice,
          linearly interpolated from ``ExtMin2`` to ``ExtMax2``.
    """
    path = Path(path)
    logger.debug("Extracting stage positions: %s", path)

    with _open(path) as f:
        img_info = f.get(f"{_INFO_ROOT}/Image")
        if img_info is None:
            raise KeyError(f"No DataSetInfo/Image group in {path}")

        ext_min_x = _attr_float(img_info, "ExtMin0", 0)
        ext_min_y = _attr_float(img_info, "ExtMin1", 0)
        ext_min_z = _attr_float(img_info, "ExtMin2", 0)
        ext_max_z = _attr_float(img_info, "ExtMax2", 0)

        n_z = _count_z_slices(f, channel=0)

    if n_z > 1:
        z_positions = list(np.linspace(ext_min_z, ext_max_z, n_z))
    elif n_z == 1:
        z_positions = [ext_min_z]
    else:
        z_positions = []

    return {
        "stage_pos_x": ext_min_x,
        "stage_pos_y": ext_min_y,
        "z_positions": z_positions,
    }


def list_ims_contents(path: str | Path) -> dict[str, Any]:
    """Walk the HDF5 hierarchy and log it for debugging.

    Parameters
    ----------
    path:
        Path to the ``.ims`` file.

    Returns
    -------
    dict
        Summary with keys:

        * **groups** (*list[str]*) -- all group paths
        * **datasets** (*list[dict]*) -- each with ``path``, ``shape``,
          ``dtype``
    """
    path = Path(path)
    groups: list[str] = []
    datasets: list[dict[str, Any]] = []

    def _visitor(name: str, obj: h5py.HLObject) -> None:
        if isinstance(obj, h5py.Group):
            groups.append(name)
            logger.info("  [group] %s", name)
        elif isinstance(obj, h5py.Dataset):
            info = {"path": name, "shape": obj.shape, "dtype": str(obj.dtype)}
            datasets.append(info)
            logger.info("  [dataset] %s  shape=%s  dtype=%s", name, obj.shape, obj.dtype)

    logger.info("Contents of IMS file: %s", path)
    with _open(path) as f:
        f.visititems(_visitor)

    return {"groups": groups, "datasets": datasets}
