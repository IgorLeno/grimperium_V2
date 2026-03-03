"""Unit tests for MOPAC optimizer input generation.

Tests PRECISE keyword enforcement and .mop file creation.
"""

import logging
from pathlib import Path

import pytest

from grimperium.crest_pm7.config import PM7Config
from grimperium.crest_pm7.mopac_optimizer import _create_mopac_input


@pytest.fixture()
def sample_xyz(tmp_path: Path) -> Path:
    """Create a minimal XYZ file for testing."""
    xyz = tmp_path / "test.xyz"
    xyz.write_text("2\ntest molecule\nH  0.0  0.0  0.0\nH  0.0  0.0  0.74\n")
    return xyz


class TestPreciseKeyword:
    """Tests that PRECISE keyword is always present in .mop files."""

    def test_precise_present_default_config(
        self, sample_xyz: Path, tmp_path: Path
    ) -> None:
        """PRECISE is present with default config (mopac_precise_scf=True)."""
        config = PM7Config()
        mop = tmp_path / "out.mop"
        assert _create_mopac_input(sample_xyz, mop, config=config)
        keywords = mop.read_text().splitlines()[0]
        assert "PRECISE" in keywords

    def test_precise_present_when_config_false(
        self, sample_xyz: Path, tmp_path: Path
    ) -> None:
        """PRECISE is still present even when config says False."""
        config = PM7Config(mopac_precise_scf=False)
        mop = tmp_path / "out.mop"
        assert _create_mopac_input(sample_xyz, mop, config=config)
        keywords = mop.read_text().splitlines()[0]
        assert "PRECISE" in keywords

    def test_warning_logged_when_config_false(
        self, sample_xyz: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """WARNING is logged when config.mopac_precise_scf=False."""
        config = PM7Config(mopac_precise_scf=False)
        mop = tmp_path / "out.mop"
        with caplog.at_level(logging.WARNING, logger="grimperium.crest_pm7.mopac_optimizer"):
            _create_mopac_input(sample_xyz, mop, config=config)
        assert any("overridden" in r.message for r in caplog.records)

    def test_no_warning_when_config_true(
        self, sample_xyz: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """No warning when config.mopac_precise_scf=True (default)."""
        config = PM7Config()
        mop = tmp_path / "out.mop"
        with caplog.at_level(logging.WARNING, logger="grimperium.crest_pm7.mopac_optimizer"):
            _create_mopac_input(sample_xyz, mop, config=config)
        assert not any("overridden" in r.message for r in caplog.records)

    def test_precise_present_without_config(
        self, sample_xyz: Path, tmp_path: Path
    ) -> None:
        """PRECISE is present even when config=None."""
        mop = tmp_path / "out.mop"
        assert _create_mopac_input(sample_xyz, mop, config=None)
        keywords = mop.read_text().splitlines()[0]
        assert "PRECISE" in keywords

    def test_scfcrt_present_with_config(
        self, sample_xyz: Path, tmp_path: Path
    ) -> None:
        """SCFCRT= is appended when config is provided."""
        config = PM7Config()
        mop = tmp_path / "out.mop"
        assert _create_mopac_input(sample_xyz, mop, config=config)
        keywords = mop.read_text().splitlines()[0]
        assert "SCFCRT=" in keywords
