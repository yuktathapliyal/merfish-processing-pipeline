"""Microscope adapter package.

Provides a factory function to get the correct adapter for a given microscope
type, and re-exports all adapter classes.
"""

from __future__ import annotations

from typing import Any

from merfish_pipeline.microscopes.andor import ANDORAdapter
from merfish_pipeline.microscopes.base import MicroscopeAdapter
from merfish_pipeline.microscopes.nikon import NIKONAdapter
from merfish_pipeline.microscopes.oni import ONIAdapter

_ADAPTERS: dict[str, type[MicroscopeAdapter]] = {
    "oni": ONIAdapter,
    "nikon": NIKONAdapter,
    "andor": ANDORAdapter,
}


def get_adapter(config: Any) -> MicroscopeAdapter:
    """Return the appropriate microscope adapter for the given config.

    Parameters
    ----------
    config:
        A ``PipelineConfig`` instance (must have ``experiment.microscope``).

    Returns
    -------
    MicroscopeAdapter
        Instantiated adapter for the configured microscope type.

    Raises
    ------
    ValueError
        If the microscope type is not recognized.
    """
    microscope_name = config.experiment.microscope
    cls = _ADAPTERS.get(microscope_name)
    if cls is None:
        available = ", ".join(sorted(_ADAPTERS))
        raise ValueError(
            f"Unknown microscope type: {microscope_name!r}. Available: {available}"
        )
    return cls(config)


__all__ = [
    "MicroscopeAdapter",
    "ONIAdapter",
    "NIKONAdapter",
    "ANDORAdapter",
    "get_adapter",
]
