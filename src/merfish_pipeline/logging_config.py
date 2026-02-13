"""Logging configuration for the merFISH processing pipeline.

Provides dual-output logging: file handler at DEBUG level and console handler
at a configurable level (INFO by default).
"""

from __future__ import annotations

import logging
import sys
from datetime import datetime
from pathlib import Path

LOGGER_NAME = "merfish_pipeline"
LOG_FORMAT = "[%(asctime)s] %(levelname)s %(name)s: %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logging(
    log_dir: Path | None = None,
    experiment_name: str = "",
    level: str = "INFO",
    console: bool = True,
    file: bool = True,
) -> logging.Logger:
    """Configure and return the root pipeline logger.

    Parameters
    ----------
    log_dir:
        Directory for log files. Required when *file* is True.
    experiment_name:
        Experiment identifier embedded in the log filename.
    level:
        Minimum log level for the console handler (e.g. "DEBUG", "INFO").
    console:
        Whether to attach a console (stderr) handler.
    file:
        Whether to attach a file handler.

    Returns
    -------
    logging.Logger
        The configured ``merfish_pipeline`` root logger.
    """
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.DEBUG)

    # Avoid duplicate handlers on repeated calls.
    logger.handlers.clear()

    formatter = logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

    if console:
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(getattr(logging, level.upper(), logging.INFO))
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    if file:
        if log_dir is None:
            raise ValueError("log_dir must be provided when file logging is enabled")
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        exp_part = f"_{experiment_name}" if experiment_name else ""
        log_filename = f"merfish_pipe{exp_part}_{timestamp}.log"

        file_handler = logging.FileHandler(log_dir / log_filename, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger
