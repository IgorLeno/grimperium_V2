"""Unit tests for MOPAC optimizer input generation.

Tests PRECISE keyword enforcement and .mop file creation.
"""

import logging
import subprocess
from pathlib import Path

import pytest

from grimperium.crest_pm7.config import PM7Config
from grimperium.crest_pm7.mopac_optimizer import (
    MOPACResult,
    _create_mopac_input,
    optimize_conformer,
    run_mopac,
)


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
            *,
            cwd: Path,
            capture_output: bool,
            text: bool,
            timeout: float,
            env: object = None,
            **kwargs: object,
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

    def test_optimize_conformer_propagates_hamiltonian(
        self, sample_xyz: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """optimize_conformer deve repassar hamiltonian/charge/multiplicity."""
        config = PM7Config(temp_dir=tmp_path, mopac_executable="mopac")
        captured: dict[str, object] = {}

        def fake_run_mopac(**kwargs: object) -> MOPACResult:
            captured.update(kwargs)
            return MOPACResult()

        monkeypatch.setattr(
            "grimperium.crest_pm7.mopac_optimizer.run_mopac",
            fake_run_mopac,
        )

        optimize_conformer(
            mol_id="mol",
            xyz_file=sample_xyz,
            config=config,
            timeout=10.0,
            hamiltonian="AM1",
            charge=1,
            multiplicity=2,
        )

        assert captured["hamiltonian"] == "AM1"
        assert captured["charge"] == 1
        assert captured["multiplicity"] == 2

    def test_invalid_multiplicity_raises_value_error(
        self, sample_xyz: Path, tmp_path: Path
    ) -> None:
        """Multiplicidade fora do mapa deve levantar ValueError imediatamente."""
        mop = tmp_path / "out.mop"
        with pytest.raises(ValueError, match="Unsupported multiplicity"):
            _create_mopac_input(sample_xyz, mop, multiplicity=7)

    def test_invalid_hamiltonian_raises_value_error(
        self, sample_xyz: Path, tmp_path: Path
    ) -> None:
        """Hamiltoniano fora do conjunto suportado deve levantar ValueError."""
        mop = tmp_path / "out.mop"
        with pytest.raises(ValueError, match="Unsupported Hamiltonian"):
            _create_mopac_input(sample_xyz, mop, hamiltonian="MNDO")

    def test_default_hamiltonian_is_pm7(
        self, sample_xyz: Path, tmp_path: Path
    ) -> None:
        """Sem hamiltonian explícito, o primeiro token deve ser PM7."""
        mop = tmp_path / "out.mop"
        assert _create_mopac_input(sample_xyz, mop)
        keywords = mop.read_text().splitlines()[0].split()
        assert keywords[0] == "PM7"
        assert "AUX" not in keywords
