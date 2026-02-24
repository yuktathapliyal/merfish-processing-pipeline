"""Shared CLI utilities."""

from __future__ import annotations

import time

import click


class StatusBanner:
    """Always-visible loading banner with elapsed time.

    Prints ``{message}...`` on entry and `` done ({elapsed}s)`` on exit.

    Usage::

        with StatusBanner("Loading pipeline"):
            heavy_imports()
            load_config()

    Output::

        Loading pipeline... done (1.2s)
    """

    def __init__(self, message: str = "Loading") -> None:
        self._message = message
        self._start: float = 0.0

    def __enter__(self) -> "StatusBanner":
        self._start = time.monotonic()
        click.echo(f"{self._message}... ", nl=False, err=True)
        return self

    def __exit__(self, *_: object) -> None:
        elapsed = time.monotonic() - self._start
        click.echo(f"done ({elapsed:.1f}s)", err=True)
