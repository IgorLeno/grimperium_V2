"""
Logging utilities for Grimperium.

This module provides configurable logging with:
    - Console and file handlers
    - Colored output for terminal
    - Progress tracking for long operations
    - Integration with tqdm for progress bars

Example:
    >>> from grimperium.utils import setup_logging, get_logger
    >>> setup_logging(level="INFO")
    >>> logger = get_logger(__name__)
    >>> logger.info("Training started...")

"""

import logging
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path
from types import TracebackType

# Default format
DEFAULT_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def _resolve_level(level: str | int) -> int:
    """Convert string level name to int."""
    if isinstance(level, int):
        return level
    numeric = getattr(logging, level.upper(), None)
    if isinstance(numeric, int):
        return numeric
    raise ValueError(f"Invalid log level: {level}")


def setup_logging(
    level: str | int = "INFO",
    log_file: str | Path | None = None,
    format_string: str = DEFAULT_FORMAT,
    date_format: str = DEFAULT_DATE_FORMAT,
    max_bytes: int = 5 * 1024 * 1024,
    backup_count: int = 3,
    console_level: str | int = "WARNING",
) -> None:
    """
    Configure global logging for Grimperium.

    Args:
        level: Logging level for file handler ('DEBUG', 'INFO', 'WARNING', 'ERROR')
        log_file: Optional file path for log output
        format_string: Log message format
        date_format: Date/time format
        max_bytes: Max size per log file before rotation (default 5MB)
        backup_count: Number of rotated backup files to keep
        console_level: Logging level for console/stderr handler

    Example:
        >>> setup_logging(level="DEBUG", log_file="grimperium.log")

    """
    resolved_level = _resolve_level(level)
    resolved_console_level = _resolve_level(console_level)
    formatter = logging.Formatter(format_string, datefmt=date_format)

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers.clear()

    if log_file is not None:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = RotatingFileHandler(
            str(log_path),
            maxBytes=max_bytes,
            backupCount=backup_count,
            encoding="utf-8",
        )
        file_handler.setLevel(resolved_level)
        file_handler.setFormatter(formatter)
        root.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(resolved_console_level)
    console_handler.setFormatter(formatter)
    root.addHandler(console_handler)


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance for a module.

    Args:
        name: Logger name (typically __name__)

    Returns:
        Configured logger instance

    Example:
        >>> logger = get_logger(__name__)
        >>> logger.info("Processing molecule...")

    """
    return logging.getLogger(name)


class ProgressLogger:
    """
    Context manager for logging progress of long operations.

    Example:
        >>> with ProgressLogger("Training", total=100) as progress:
        ...     for i in range(100):
        ...         # do work
        ...         progress.update(1)

    """

    def __init__(
        self,
        description: str,
        total: int | None = None,
        unit: str = "it",
        disable: bool = False,
    ) -> None:
        """
        Initialize ProgressLogger.

        Args:
            description: Description of the operation
            total: Total number of iterations
            unit: Unit name for progress display
            disable: Disable progress output

        """
        self.description = description
        self.total = total
        self.unit = unit
        self.disable = disable
        self._progress = 0

    def __enter__(self) -> "ProgressLogger":
        """Enter context."""
        raise NotImplementedError("Will be implemented in Batch 5")

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        """Exit context."""
        raise NotImplementedError("Will be implemented in Batch 5")

    def update(self, n: int = 1) -> None:
        """Update progress by n steps."""
        raise NotImplementedError("Will be implemented in Batch 5")

    def set_description(self, description: str) -> None:
        """Update progress description."""
        raise NotImplementedError("Will be implemented in Batch 5")


def log_metrics(
    metrics: dict[str, float],
    prefix: str = "",
    logger: logging.Logger | None = None,
) -> None:
    """
    Log a dictionary of metrics in a formatted way.

    Args:
        metrics: Dictionary of metric names to values
        prefix: Prefix for log message
        logger: Logger instance (uses default if None)

    Example:
        >>> metrics = {"rmse": 0.5, "mae": 0.3, "r2": 0.95}
        >>> log_metrics(metrics, prefix="Test")
        # Output: Test | RMSE: 0.5000 | MAE: 0.3000 | R²: 0.9500

    """
    raise NotImplementedError("Will be implemented in Batch 5")
