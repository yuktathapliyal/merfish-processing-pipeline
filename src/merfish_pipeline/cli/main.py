"""Root CLI group for ``merfish-pipe``."""

from __future__ import annotations

import click

from merfish_pipeline import __version__


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
