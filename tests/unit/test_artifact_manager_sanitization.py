"""
Unit tests for ArtifactManager._sanitize_id() path injection protection.

Tests verify that dangerous characters, path traversal attempts, and edge cases
are properly sanitized before being used in filesystem operations.
"""

import pytest
from pathlib import Path
from src.grimperium.crest_pm7.batch.artifact_manager import ArtifactManager


@pytest.fixture
def manager():
    """Create a temporary ArtifactManager for testing."""
    temp_dir = Path("/tmp/test_artifacts")
    temp_dir.mkdir(exist_ok=True)
    return ArtifactManager(artifact_dir=temp_dir)


class TestSanitizeId:
    """Test suite for _sanitize_id() method."""

    def test_sanitize_id_normal_input(self, manager):
        """Normal alphanumeric ID should pass through unchanged."""
        result = manager._sanitize_id("mol_12345")
        assert result == "mol_12345"

    def test_sanitize_id_removes_forward_slash(self, manager):
        """Forward slashes should be replaced with underscores."""
        result = manager._sanitize_id("mol/test")
        assert "/" not in result
        assert "_" in result

    def test_sanitize_id_removes_backslash(self, manager):
        """Backslashes should be replaced with underscores."""
        result = manager._sanitize_id("mol\\test")
        assert "\\" not in result
        assert "_" in result

    def test_sanitize_id_blocks_directory_traversal_dots(self, manager):
        """Double dots (..) should be replaced to prevent directory traversal."""
        result = manager._sanitize_id("mol..passwd")
        assert ".." not in result
        # Should contain underscores instead
        assert "_" in result

    def test_sanitize_id_blocks_parent_directory_reference(self, manager):
        """Path like ../etc/passwd should be neutralized."""
        result = manager._sanitize_id("../../../etc/passwd")
        # Should not contain ".." or slashes
        assert ".." not in result
        assert "/" not in result

    def test_sanitize_id_removes_null_byte(self, manager):
        """Null bytes should be removed."""
        result = manager._sanitize_id("mol\x00test")
        assert "\x00" not in result

    def test_sanitize_id_handles_windows_reserved_con(self, manager):
        """Windows reserved name 'CON' should be prefixed with underscore."""
        result = manager._sanitize_id("CON")
        assert result == "_CON"

    def test_sanitize_id_handles_windows_reserved_prn(self, manager):
        """Windows reserved name 'PRN' should be prefixed with underscore."""
        result = manager._sanitize_id("PRN")
        assert result == "_PRN"

    def test_sanitize_id_handles_windows_reserved_nul(self, manager):
        """Windows reserved name 'NUL' should be prefixed with underscore."""
        result = manager._sanitize_id("NUL")
        assert result == "_NUL"

    def test_sanitize_id_handles_windows_reserved_com1(self, manager):
        """Windows reserved name 'COM1' should be prefixed with underscore."""
        result = manager._sanitize_id("COM1")
        assert result == "_COM1"

    def test_sanitize_id_handles_windows_reserved_lpt1(self, manager):
        """Windows reserved name 'LPT1' should be prefixed with underscore."""
        result = manager._sanitize_id("LPT1")
        assert result == "_LPT1"

    def test_sanitize_id_truncates_long_ids(self, manager):
        """IDs longer than 255 chars should be truncated."""
        long_id = "a" * 300
        result = manager._sanitize_id(long_id)
        assert len(result) <= 255

    def test_sanitize_id_preserves_valid_long_id(self, manager):
        """Valid long ID should be truncated to exactly 255 chars."""
        long_id = "valid_mol_" + "x" * 250
        result = manager._sanitize_id(long_id)
        assert len(result) == 255
        assert result.startswith("valid_mol_")

    def test_sanitize_id_removes_multiple_dots(self, manager):
        """Multiple dots should be replaced."""
        result = manager._sanitize_id("mol...test")
        assert "..." not in result

    def test_sanitize_id_handles_leading_dot(self, manager):
        """Leading dot should be handled (hidden files on Unix)."""
        result = manager._sanitize_id(".hidden")
        # Should not be a hidden file after sanitization
        assert not result.startswith(".")

    def test_sanitize_id_handles_trailing_dot(self, manager):
        """Trailing dot should be handled."""
        result = manager._sanitize_id("mol_test.")
        # Trailing dot removed or replaced
        assert not result.endswith(".")

    def test_sanitize_id_colon_removal(self, manager):
        """Colons (invalid on Windows) should be removed."""
        result = manager._sanitize_id("mol:test")
        assert ":" not in result

    def test_sanitize_id_asterisk_removal(self, manager):
        """Asterisks (wildcard char) should be removed."""
        result = manager._sanitize_id("mol*test")
        assert "*" not in result

    def test_sanitize_id_question_mark_removal(self, manager):
        """Question marks (wildcard char) should be removed."""
        result = manager._sanitize_id("mol?test")
        assert "?" not in result

    def test_sanitize_id_angle_brackets_removal(self, manager):
        """Angle brackets should be removed."""
        result = manager._sanitize_id("mol<test>")
        assert "<" not in result
        assert ">" not in result

    def test_sanitize_id_pipe_removal(self, manager):
        """Pipes should be removed."""
        result = manager._sanitize_id("mol|test")
        assert "|" not in result

    def test_sanitize_id_quote_removal(self, manager):
        """Quotes should be removed."""
        result = manager._sanitize_id('mol"test')
        assert '"' not in result

    def test_sanitize_id_real_world_injection_attempt_1(self, manager):
        """Real attack: ../../../etc/passwd should be blocked."""
        result = manager._sanitize_id("../../../etc/passwd")
        # Should not be able to traverse up directories
        assert ".." not in result
        assert "/" not in result

    def test_sanitize_id_real_world_injection_attempt_2(self, manager):
        """Real attack: mol_id with null byte should be blocked."""
        result = manager._sanitize_id("mol_1\x00.txt")
        assert "\x00" not in result

    def test_sanitize_id_returns_string(self, manager):
        """Result should always be a string."""
        result = manager._sanitize_id("test")
        assert isinstance(result, str)

    def test_sanitize_id_empty_string_handling(self, manager):
        """Empty string should return something safe."""
        result = manager._sanitize_id("")
        # Should not be empty or be "_empty_"
        assert result is not None
        assert len(result) > 0
