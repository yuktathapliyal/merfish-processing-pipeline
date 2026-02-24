"""Root CLI group for ``merfish-pipe``."""

from __future__ import annotations

import sys
import time

# Print immediately so the user knows the tool is starting, before any
# heavy imports (pydantic, yaml, etc.) that can take 30-60s on cold servers.
_boot_start = time.monotonic()
print("Starting merfish-pipe... ", end="", flush=True, file=sys.stderr)

import click  # noqa: E402

from merfish_pipeline import __version__  # noqa: E402

_boot_elapsed = time.monotonic() - _boot_start
print(f"ready ({_boot_elapsed:.1f}s)", flush=True, file=sys.stderr)


@click.group()
@click.version_option(version=__version__, prog_name="merfish-pipe")
def cli() -> None:
    """merFISH Processing Pipeline — standardized pipeline for ANDOR, NIKON, and ONI microscopes."""


# Lazy imports to keep --help fast and avoid circular imports.

def _register_commands() -> None:
    from merfish_pipeline.cli.config_cmd import config  # noqa: F811
    from merfish_pipeline.cli.run import run  # noqa: F811
    from merfish_pipeline.cli.status import status  # noqa: F811

    cli.add_command(run)
    cli.add_command(config)
    cli.add_command(status)


_register_commands()
