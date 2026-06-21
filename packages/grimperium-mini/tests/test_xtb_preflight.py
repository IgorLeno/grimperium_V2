"""Tests for xTB preflight check and per-molecule status attribution."""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from grimperium_mini.config import MiniConfig
from grimperium_mini.io import load_molecules
from grimperium_mini.pipeline import process_molecule, run_pipeline
from grimperium_mini.xtb import verify_xtb_runtime

DATA = Path("data/grimperium_mini_pipeline_tcc.xlsx")


# ── verify_xtb_runtime unit tests ────────────────────────────────────────────


def test_verify_xtb_runtime_success(monkeypatch: pytest.MonkeyPatch) -> None:
    def _ok(cmd: object, *, cwd: object, capture_output: object, text: object, timeout: object, check: object) -> MagicMock:
        (Path(str(cwd)) / "xtbopt.xyz").write_text("1\ntest\nO 0 0 0\n")
        result = MagicMock()
        result.returncode = 0
        result.stderr = ""
        return result

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _ok)
    ok, reason = verify_xtb_runtime("xtb", threads=1)
    assert ok is True
    assert reason == ""


def test_verify_xtb_runtime_rc_nonzero(monkeypatch: pytest.MonkeyPatch) -> None:
    def _fail(cmd: object, *, cwd: object, capture_output: object, text: object, timeout: object, check: object) -> MagicMock:
        result = MagicMock()
        result.returncode = 2
        result.stderr = "Fortran runtime error: Missing comma between descriptors"
        return result

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _fail)
    ok, reason = verify_xtb_runtime("xtb", threads=1)
    assert ok is False
    assert "rc=2" in reason


def test_verify_xtb_runtime_missing_output(monkeypatch: pytest.MonkeyPatch) -> None:
    """rc=0 but xtbopt.xyz not produced → failure."""

    def _no_output(cmd: object, *, cwd: object, capture_output: object, text: object, timeout: object, check: object) -> MagicMock:
        result = MagicMock()
        result.returncode = 0
        result.stderr = ""
        return result

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _no_output)
    ok, reason = verify_xtb_runtime("xtb", threads=1)
    assert ok is False
    assert "xtbopt.xyz" in reason


def test_verify_xtb_runtime_not_found(monkeypatch: pytest.MonkeyPatch) -> None:
    def _not_found(*args: object, **kwargs: object) -> None:
        raise FileNotFoundError("xtb: command not found")

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _not_found)
    ok, reason = verify_xtb_runtime("xtb_missing")
    assert ok is False
    assert "not found" in reason


def test_verify_xtb_runtime_timeout(monkeypatch: pytest.MonkeyPatch) -> None:
    def _timeout(*args: object, **kwargs: object) -> None:
        raise subprocess.TimeoutExpired("xtb", 60)

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _timeout)
    ok, reason = verify_xtb_runtime("xtb")
    assert ok is False
    assert "timed out" in reason


# ── run_pipeline preflight integration ───────────────────────────────────────

_DRY_ROW: dict[str, object] = {
    "mol_id": "m0",
    "molecule": "",
    "smiles": "",
    "n_conformers_founded": 0,
    "hf_am1_kJmol": "",
    "hf_pm3_kJmol": "",
    "hf_pm7_kJmol": "",
    "status": "dry_run",
    "crest_status": "dry_run",
    "mopac_status": "dry_run",
    "hf_exp_kJmol": "",
    "total_time_s": 0.0,
    "total_time_min": 0.0,
}


@pytest.mark.skipif(not DATA.exists(), reason="xlsx fixture not available")
def test_run_pipeline_disables_xtb_on_preflight_failure(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """When preflight fails, config.xtb_enabled is set to False and warning logged."""
    import logging

    monkeypatch.setattr(
        "grimperium_mini.pipeline.verify_xtb_runtime",
        lambda *a, **kw: (False, "Fortran runtime error"),
    )
    monkeypatch.setattr(
        "grimperium_mini.pipeline.process_molecule",
        lambda mol, config, dry_run=False: _DRY_ROW,
    )
    monkeypatch.setattr("grimperium_mini.pipeline.append_summary_row", lambda *a, **kw: None)

    config = MiniConfig(
        work_root=tmp_path / "runs",
        results_dir=tmp_path / "results",
        xtb_enabled=True,
    )

    with caplog.at_level(logging.WARNING, logger="grimperium_mini.pipeline"):
        run_pipeline(DATA, config, limit=1, dry_run=False)

    assert config.xtb_enabled is False
    assert any("xTB pre-optimization disabled" in r.message for r in caplog.records)


# ── process_molecule status attribution ──────────────────────────────────────


@pytest.mark.skipif(not DATA.exists(), reason="xlsx fixture not available")
def test_process_molecule_xtb_failure_crest_not_attempted(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """When xTB fails, crest_status reflects xTB failure and MOPAC is not_attempted."""

    def _xtb_fail(cmd: object, *, cwd: object, capture_output: object, text: object, timeout: object, check: object) -> MagicMock:
        result = MagicMock()
        result.returncode = 2
        result.stderr = "Fortran runtime error"
        return result

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _xtb_fail)

    mols = load_molecules(DATA, limit=1)
    config = MiniConfig(
        work_root=tmp_path / "runs",
        results_dir=tmp_path / "results",
        xtb_enabled=True,
    )
    row = process_molecule(mols[0], config, dry_run=False)

    assert row["crest_status"] == "xtb_failed"
    assert row["mopac_status"] == "not_attempted"
    assert row["status"] == "failed"
