"""Unit tests for MOPAC optimizer input generation.

Tests PRECISE keyword enforcement and .mop file creation.
"""

import logging
import subprocess
from pathlib import Path

import pytest

from grimperium.crest_pm7.config import PM7Config
from grimperium.crest_pm7.mopac_optimizer import _create_mopac_input, run_mopac


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
        with caplog.at_level(
            logging.WARNING, logger="grimperium.crest_pm7.mopac_optimizer"
        ):
            _create_mopac_input(sample_xyz, mop, config=config)
        assert any("overridden" in r.message for r in caplog.records)

    def test_no_warning_when_config_true(
        self, sample_xyz: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """No warning when config.mopac_precise_scf=True (default)."""
        config = PM7Config()
        mop = tmp_path / "out.mop"
        with caplog.at_level(
            logging.WARNING, logger="grimperium.crest_pm7.mopac_optimizer"
        ):
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

    def test_scfcrt_present_with_config(self, sample_xyz: Path, tmp_path: Path) -> None:
        """SCFCRT= is appended when config is provided."""
        config = PM7Config()
        mop = tmp_path / "out.mop"
        assert _create_mopac_input(sample_xyz, mop, config=config)
        keywords = mop.read_text().splitlines()[0]
        assert "SCFCRT=" in keywords


class TestParameterizedMopacInput:
    """Tests for Hamiltonian and charge-state keyword generation."""

    def test_custom_hamiltonian_replaces_pm7_without_aux_by_default(
        self, sample_xyz: Path, tmp_path: Path
    ) -> None:
        """A caller can request AM1 without inheriting the Mini AUX default."""
        mop = tmp_path / "out.mop"

        assert _create_mopac_input(
            sample_xyz,
            mop,
            config=None,
            hamiltonian="AM1",
            extra_keywords=["GNORM=0.01"],
        )

        keywords = mop.read_text().splitlines()[0].split()
        assert keywords[0] == "AM1"
        assert "PM7" not in keywords
        assert "EF" in keywords
        assert "PRECISE" in keywords
        assert "GNORM=0.01" in keywords
        assert "AUX" not in keywords

    def test_run_mopac_passes_hamiltonian_charge_and_multiplicity_to_input(
        self, sample_xyz: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """run_mopac propagates calculation parameters into the generated .mop."""
        config = PM7Config(temp_dir=tmp_path, mopac_executable="mopac")
        work_dir = tmp_path / "mopac-run"

        def fake_run(
            cmd: list[str],
            cwd: Path,
            capture_output: bool,
            text: bool,
            timeout: float,
        ) -> subprocess.CompletedProcess[str]:
            mop_file = Path(cmd[1])
            mop_file.with_suffix(".out").write_text(
                "FINAL HEAT OF FORMATION = -10.0 KCAL/MOL\n",
                encoding="utf-8",
            )
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        monkeypatch.setattr(
            "grimperium.crest_pm7.mopac_optimizer.subprocess.run",
            fake_run,
        )

        result = run_mopac(
            mol_id="mol",
            xyz_file=sample_xyz,
            config=config,
            timeout=10.0,
            work_dir=work_dir,
            hamiltonian="PM3",
            charge=-1,
            multiplicity=3,
        )

        assert result.hof == -10.0
        keywords = (work_dir / "mol_conf000.mop").read_text().splitlines()[0].split()
        assert keywords[0] == "PM3"
        assert "CHARGE=-1" in keywords
        assert "TRIPLET" in keywords
        assert "AUX" not in keywords
