"""Tests for xTB preflight check and per-task status attribution."""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock

from grimperium_mini.config import MiniConfig
from grimperium_mini.io import load_pipeline_tasks
from grimperium_mini.pipeline import process_task, run_pipeline
from grimperium_mini.xtb import verify_xtb_runtime

DATA = Path("data/grimperium_mini_pipeline_tcc.xlsx")


# ── verify_xtb_runtime unit tests ────────────────────────────────────────────


def test_verify_xtb_runtime_success(monkeypatch):
    def _ok(cmd, *, cwd, capture_output, text, timeout, check):
        (Path(cwd) / "xtbopt.xyz").write_text("1\ntest\nO 0 0 0\n")
        result = MagicMock()
        result.returncode = 0
        result.stderr = ""
        return result

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _ok)
    ok, reason = verify_xtb_runtime("xtb", threads=1)
    assert ok is True
    assert reason == ""


def test_verify_xtb_runtime_rc_nonzero(monkeypatch):
    def _fail(cmd, *, cwd, capture_output, text, timeout, check):
        result = MagicMock()
        result.returncode = 2
        result.stderr = "Fortran runtime error: Missing comma between descriptors"
        return result

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _fail)
    ok, reason = verify_xtb_runtime("xtb", threads=1)
    assert ok is False
    assert "rc=2" in reason


def test_verify_xtb_runtime_missing_output(monkeypatch):
    """rc=0 but xtbopt.xyz not produced → failure."""

    def _no_output(cmd, *, cwd, capture_output, text, timeout, check):
        result = MagicMock()
        result.returncode = 0
        result.stderr = ""
        return result

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _no_output)
    ok, reason = verify_xtb_runtime("xtb", threads=1)
    assert ok is False
    assert "xtbopt.xyz" in reason


def test_verify_xtb_runtime_not_found(monkeypatch):
    def _not_found(*args, **kwargs):
        raise FileNotFoundError("xtb: command not found")

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _not_found)
    ok, reason = verify_xtb_runtime("xtb_missing")
    assert ok is False
    assert "not found" in reason


def test_verify_xtb_runtime_timeout(monkeypatch):
    def _timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired("xtb", 60)

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _timeout)
    ok, reason = verify_xtb_runtime("xtb")
    assert ok is False
    assert "timed out" in reason


# ── run_pipeline preflight integration ───────────────────────────────────────


def test_run_pipeline_disables_xtb_on_preflight_failure(monkeypatch, tmp_path, caplog):
    """When preflight fails, config.xtb_enabled is set to False and warning logged."""
    import logging

    monkeypatch.setattr(
        "grimperium_mini.pipeline.verify_xtb_runtime",
        lambda *a, **kw: (False, "Fortran runtime error"),
    )
    # Stub process_task so the test doesn't need real external tools.
    monkeypatch.setattr(
        "grimperium_mini.pipeline.process_task",
        lambda task, config, dry_run=False: ({}, []),
    )

    config = MiniConfig(
        work_root=tmp_path / "runs",
        results_dir=tmp_path / "results",
        xtb_enabled=True,
    )

    with caplog.at_level(logging.WARNING, logger="grimperium_mini.pipeline"):
        run_pipeline(DATA, ["AM1"], config, limit=1, dry_run=False)

    assert config.xtb_enabled is False
    assert any("xTB pre-optimization disabled" in r.message for r in caplog.records)


# ── process_task status attribution ──────────────────────────────────────────


def test_process_task_xtb_failure_leaves_crest_not_attempted(monkeypatch, tmp_path):
    """When xTB fails, crest_status must stay 'not_attempted' (CREST never ran)."""

    def _xtb_fail(cmd, *, cwd, capture_output, text, timeout, check):
        result = MagicMock()
        result.returncode = 2
        result.stderr = "Fortran runtime error"
        return result

    monkeypatch.setattr("grimperium_mini.xtb.subprocess.run", _xtb_fail)

    task = load_pipeline_tasks(DATA, methods=["AM1"], limit=1)[0]
    config = MiniConfig(
        work_root=tmp_path / "runs",
        results_dir=tmp_path / "results",
        xtb_enabled=True,
    )
    row, _ = process_task(task, config, dry_run=False)

    assert row["xtb_status"] == "failed"
    assert row["crest_status"] == "not_attempted"
    assert row["mopac_status"] == "not_attempted"
    assert row["status"] == "failed"
