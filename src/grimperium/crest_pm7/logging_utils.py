"""Structured logging utilities for CREST-PM7 Pipeline.

Provides JSONL logging for analysis and text logging for debugging.
"""

import json
import logging
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config import PM7Config

# Reserved keys that should not be overwritten by extra_data
_RESERVED_KEYS = {"timestamp", "level", "logger", "message"}


def _sanitize_session_name(session_name: str) -> str:
    """Sanitize session name to avoid logger hierarchy issues.

    Args:
        session_name: Raw session name (may contain dots, special chars)

    Returns:
        Sanitized session name (alphanumeric, underscore, hyphen only)
    """
    # Replace non-alphanumeric characters (except underscore and hyphen) with underscores
    sanitized = re.sub(r"[^a-zA-Z0-9_-]", "_", session_name).strip("_-")
    # Return "default" if empty after sanitization
    return sanitized if sanitized else "default"


class StructuredLogHandler(logging.Handler):
    """Handler that writes structured JSONL logs for analysis."""

    def __init__(self, log_file: Path) -> None:
        """Initialize the handler.

        Args:
            log_file: Path to the JSONL log file
        """
        super().__init__()
        self.log_file = log_file
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

    def emit(self, record: logging.LogRecord) -> None:
        """Write a log record as JSONL."""
        try:
            log_entry = {
                "timestamp": datetime.now(timezone.utc).isoformat(
                    timespec="milliseconds"
                ),
                "level": record.levelname,
                "logger": record.name,
                "message": record.getMessage(),
            }

            # Add extra fields if present
            if hasattr(record, "mol_id"):
                log_entry["mol_id"] = record.mol_id
            if hasattr(record, "smiles"):
                log_entry["smiles"] = record.smiles
            if hasattr(record, "grade"):
                log_entry["grade"] = record.grade
            if hasattr(record, "hof"):
                log_entry["hof"] = record.hof
            if hasattr(record, "success"):
                log_entry["success"] = record.success

            # Merge extra_data filtering out reserved keys
            if hasattr(record, "extra_data"):
                filtered = {
                    k: v
                    for k, v in record.extra_data.items()
                    if k not in _RESERVED_KEYS and k not in log_entry
                }
                log_entry.update(filtered)

            with open(self.log_file, "a", encoding="utf-8") as f:
                f.write(json.dumps(log_entry, ensure_ascii=False) + "\n")
        except Exception:
            self.handleError(record)


def setup_logging(
    config: PM7Config,
    session_name: str | None = None,
) -> logging.Logger:
    """Set up logging for a pipeline session.

    Creates both JSONL (for analysis) and text (for debugging) handlers.
    Uses a session-specific child logger to avoid clearing shared handlers.

    Args:
        config: Pipeline configuration
        session_name: Optional session identifier

    Returns:
        Configured logger for the pipeline
    """
    if session_name is None:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        session_name = f"phase_{config.phase.value if hasattr(config.phase, 'value') else config.phase}_{timestamp}"

    # Sanitize session name to avoid logger hierarchy issues
    sanitized_name = _sanitize_session_name(session_name)

    log_dir = config.output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # JSONL log for analysis
    jsonl_file = log_dir / f"{session_name}.jsonl"

    # Text log for debugging
    text_file = log_dir / f"{session_name}.log"

    # Create a child logger specific to this session to avoid clearing shared logger handlers
    logger = logging.getLogger(f"grimperium.crest_pm7.{sanitized_name}")
    logger.setLevel(logging.DEBUG)

    # Disable propagation to parent loggers (prevents duplicate entries)
    logger.propagate = False

    # Clear existing handlers only on this session-specific logger
    logger.handlers.clear()

    # JSONL handler
    jsonl_handler = StructuredLogHandler(jsonl_file)
    jsonl_handler.setLevel(logging.INFO)
    logger.addHandler(jsonl_handler)

    # Text file handler
    text_handler = logging.FileHandler(text_file, encoding="utf-8")
    text_handler.setLevel(logging.DEBUG)
    text_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    text_handler.setFormatter(text_formatter)
    logger.addHandler(text_handler)

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter("%(levelname)s - %(message)s")
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    return logger


def log_molecule_start(
    logger: logging.Logger,
    mol_id: str,
    smiles: str,
) -> None:
    """Log the start of molecule processing."""
    extra = {"mol_id": mol_id, "smiles": smiles}
    record = logger.makeRecord(
        logger.name,
        logging.INFO,
        "",
        0,
        "mol_start",
        (),
        None,
    )
    for key, value in extra.items():
        setattr(record, key, value)
    logger.handle(record)


def log_molecule_complete(
    logger: logging.Logger,
    mol_id: str,
    grade: str,
    hof: float | None,
    success: bool,
) -> None:
    """Log the completion of molecule processing."""
    extra = {"mol_id": mol_id, "grade": grade, "hof": hof, "success": success}
    record = logger.makeRecord(
        logger.name,
        logging.INFO,
        "",
        0,
        "mol_complete",
        (),
        None,
    )
    for key, value in extra.items():
        setattr(record, key, value)
    logger.handle(record)


def log_with_extra(
    logger: logging.Logger,
    level: int,
    message: str,
    extra_data: dict[str, Any],
) -> None:
    """Log a message with additional structured data."""
    record = logger.makeRecord(
        logger.name,
        level,
        "",
        0,
        message,
        (),
        None,
    )
    record.extra_data = extra_data
    logger.handle(record)
