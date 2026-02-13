"""NIKON microscope adapter.

NIKON raw data follows the same directory layout as ONI — TIFF files in
wavelength-specific subdirectories.  The adapter reuses the ONI parsing
logic but applies NIKON-specific defaults.
"""

from __future__ import annotations

from merfish_pipeline.microscopes.oni import ONIAdapter


class NIKONAdapter(ONIAdapter):
    """Adapter for NIKON raw data.

    NIKON uses the same file organization as ONI (wavelength subdirectories
    containing ``merFISH_{ir}_{fov}_{z}.TIFF``), so this adapter inherits
    all ONI logic and only overrides the ``name``.

    Any NIKON-specific behaviour (e.g. different position file formats)
    can be added here in the future.
    """

    name = "nikon"
