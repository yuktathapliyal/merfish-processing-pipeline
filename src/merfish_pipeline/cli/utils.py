"""Shared CLI utilities."""

from __future__ import annotations

import itertools
import sys
import threading
import time


class Spinner:
    """Simple threaded spinner for CLI feedback during slow operations.

    Only displays the spinner when stderr is a real terminal (TTY).
    In non-TTY environments (pipes, CI, etc.) this is a silent no-op.

    Usage::

        with Spinner("Loading pipeline"):
            heavy_imports()
            load_config()
    """

    def __init__(self, message: str = "Loading") -> None:
        self._message = message
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._active = sys.stderr.isatty()

    def __enter__(self) -> "Spinner":
        if self._active:
            self._thread = threading.Thread(target=self._spin, daemon=True)
            self._thread.start()
        return self

    def __exit__(self, *_: object) -> None:
        if not self._active:
            return
        self._stop.set()
        if self._thread is not None:
            self._thread.join()
        # Clear the spinner line
        sys.stderr.write("\r" + " " * (len(self._message) + 4) + "\r")
        sys.stderr.flush()

    def _spin(self) -> None:
        for ch in itertools.cycle("|/-\\"):
            if self._stop.is_set():
                break
            sys.stderr.write(f"\r{self._message} {ch}")
            sys.stderr.flush()
            time.sleep(0.1)
