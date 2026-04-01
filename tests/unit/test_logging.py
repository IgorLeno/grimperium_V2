"""Tests for grimperium.utils.logging module."""

import logging
import re
import sys
from collections.abc import Generator
from logging.handlers import RotatingFileHandler
from pathlib import Path

import pytest

from grimperium.utils.logging import get_logger, setup_logging


@pytest.fixture(autouse=True)
def _clean_root_logger() -> Generator[None, None, None]:
    """Isolate root logger state between tests."""
    root = logging.getLogger()
    original_handlers = root.handlers[:]
    original_level = root.level
    root.handlers.clear()
    yield
    root.handlers.clear()
    root.handlers.extend(original_handlers)
    root.level = original_level


class TestSetupLogging:
    """Tests for setup_logging() configuration."""

    def test_creates_console_handler(self) -> None:
        """Root logger gets a StreamHandler after setup."""
        setup_logging()
        root = logging.getLogger()
        stream_handlers = [
            h
            for h in root.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, RotatingFileHandler)
        ]
        assert len(stream_handlers) >= 1

    def test_console_handler_level_default_warning(self) -> None:
        """Console handler defaults to WARNING level."""
        setup_logging()
        root = logging.getLogger()
        console = next(
            h
            for h in root.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, RotatingFileHandler)
        )
        assert console.level == logging.WARNING

    def test_console_handler_writes_to_stderr(self) -> None:
        """Console handler stream is sys.stderr."""
        setup_logging()
        root = logging.getLogger()
        console = next(
            h
            for h in root.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, RotatingFileHandler)
        )
        assert console.stream is sys.stderr

    def test_with_log_file_creates_rotating_handler(self, tmp_path: Path) -> None:
        """RotatingFileHandler is added when log_file is provided."""
        log_file = tmp_path / "test.log"
        setup_logging(log_file=log_file)
        root = logging.getLogger()
        rotating = [h for h in root.handlers if isinstance(h, RotatingFileHandler)]
        assert len(rotating) == 1

    def test_file_handler_level(self, tmp_path: Path) -> None:
        """File handler level matches the level parameter."""
        log_file = tmp_path / "test.log"
        setup_logging(level="DEBUG", log_file=log_file)
        root = logging.getLogger()
        rotating = next(h for h in root.handlers if isinstance(h, RotatingFileHandler))
        assert rotating.level == logging.DEBUG

    def test_file_handler_rotation_config(self, tmp_path: Path) -> None:
        """File handler uses default rotation: 5MB, 3 backups."""
        log_file = tmp_path / "test.log"
        setup_logging(log_file=log_file)
        root = logging.getLogger()
        rotating = next(h for h in root.handlers if isinstance(h, RotatingFileHandler))
        assert rotating.maxBytes == 5 * 1024 * 1024
        assert rotating.backupCount == 3

    def test_creates_log_directory(self, tmp_path: Path) -> None:
        """Parent directory for log file is created automatically."""
        log_file = tmp_path / "subdir" / "deep" / "test.log"
        setup_logging(log_file=log_file)
        assert log_file.parent.exists()

    def test_custom_rotation_params(self, tmp_path: Path) -> None:
        """Custom max_bytes and backup_count are respected."""
        log_file = tmp_path / "test.log"
        setup_logging(log_file=log_file, max_bytes=1024, backup_count=7)
        root = logging.getLogger()
        rotating = next(h for h in root.handlers if isinstance(h, RotatingFileHandler))
        assert rotating.maxBytes == 1024
        assert rotating.backupCount == 7

    def test_idempotent_no_duplicate_handlers(self, tmp_path: Path) -> None:
        """Calling setup_logging twice does not duplicate handlers."""
        log_file = tmp_path / "test.log"
        setup_logging(log_file=log_file)
        setup_logging(log_file=log_file)
        root = logging.getLogger()
        assert len(root.handlers) == 2  # 1 file + 1 console

    def test_custom_format_string(self) -> None:
        """Formatter uses the provided format string."""
        custom_fmt = "%(levelname)s - %(message)s"
        setup_logging(format_string=custom_fmt)
        root = logging.getLogger()
        console = next(
            h
            for h in root.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, RotatingFileHandler)
        )
        assert console.formatter is not None
        assert console.formatter._fmt == custom_fmt

    def test_custom_console_level(self) -> None:
        """console_level='ERROR' sets console handler to ERROR."""
        setup_logging(console_level="ERROR")
        root = logging.getLogger()
        console = next(
            h
            for h in root.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, RotatingFileHandler)
        )
        assert console.level == logging.ERROR

    def test_no_file_handler_when_log_file_none(self) -> None:
        """No RotatingFileHandler when log_file is None."""
        setup_logging(log_file=None)
        root = logging.getLogger()
        rotating = [h for h in root.handlers if isinstance(h, RotatingFileHandler)]
        assert len(rotating) == 0


class TestGetLogger:
    """Tests for get_logger()."""

    def test_returns_logger_instance(self) -> None:
        """Returns a logging.Logger instance."""
        logger = get_logger("test.module")
        assert isinstance(logger, logging.Logger)

    def test_logger_name_matches(self) -> None:
        """Logger name matches the input."""
        logger = get_logger("grimperium.utils.foo")
        assert logger.name == "grimperium.utils.foo"

    def test_same_logger_for_same_name(self) -> None:
        """Same name returns the same logger object."""
        logger1 = get_logger("grimperium.test")
        logger2 = get_logger("grimperium.test")
        assert logger1 is logger2


class TestLoggingIntegration:
    """Integration tests for the full logging pipeline."""

    def test_debug_message_written_to_file(self, tmp_path: Path) -> None:
        """DEBUG message appears in the log file."""
        log_file = tmp_path / "test.log"
        setup_logging(level="DEBUG", log_file=log_file)
        logger = get_logger("integration.test")
        logger.debug("test debug message")

        # Flush handlers
        for h in logging.getLogger().handlers:
            h.flush()

        content = log_file.read_text()
        assert "test debug message" in content

    def test_console_filters_debug(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """DEBUG message does NOT reach the WARNING-level console handler."""
        setup_logging(console_level="WARNING")
        logger = get_logger("integration.filter")
        logger.debug("should not appear")

        captured = capsys.readouterr()
        assert "should not appear" not in captured.err

    def test_format_includes_timestamp_and_level(self, tmp_path: Path) -> None:
        """Log line matches expected format with timestamp and level."""
        log_file = tmp_path / "test.log"
        setup_logging(level="DEBUG", log_file=log_file)
        logger = get_logger("integration.format")
        logger.info("format check")

        for h in logging.getLogger().handlers:
            h.flush()

        content = log_file.read_text()
        # Expected: "2026-03-18 12:00:00 | INFO     | integration.format | format check"
        pattern = r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} \| INFO\s+\| integration\.format \| format check"
        assert re.search(
            pattern, content
        ), f"Log content did not match expected format:\n{content}"
