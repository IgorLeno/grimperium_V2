"""Regression tests for baseline reference-comparison logic.

These tests cover the logic moved out of the operational core
(result_evaluator.py) into baseline_evaluation.py.  They use stub
PM7Result objects and in-memory baseline dicts.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from grimperium.crest_pm7.config import QualityGrade
from tests.regression.baseline_evaluation import (
    TOLERANCE_ABSOLUTE,
    BaselineEvaluator,
    PhaseABaselineResult,
)

# ── Helpers ───────────────────────────────────────────────────────────────────


def _pm7_result(
    mol_id: str,
    success: bool = True,
    hof: float | None = -65.0,
    grade: QualityGrade = QualityGrade.A,
) -> MagicMock:
    """Build a minimal PM7Result stub."""
    r = MagicMock()
    r.mol_id = mol_id
    r.success = success
    r.most_stable_hof = hof
    r.quality_grade = grade
    return r


# ── TOLERANCE_ABSOLUTE ────────────────────────────────────────────────────────


def test_tolerance_absolute_value() -> None:
    assert TOLERANCE_ABSOLUTE == 2.5


# ── load_baseline ─────────────────────────────────────────────────────────────


def test_load_baseline_reads_molecules_and_criteria(tmp_path: Path) -> None:
    baseline = {
        "molecules": {"m1": {"hof_value": -65.0}},
        "phase_a_success_criteria": {"min_success_rate": 0.9},
        "tolerance_kcal_mol": 3.0,
    }
    p = tmp_path / "baseline.json"
    p.write_text(json.dumps(baseline), encoding="utf-8")

    ev = BaselineEvaluator()
    assert ev.load_baseline(p) is True
    assert ev.baseline_loaded is True
    assert "m1" in ev.baseline
    assert ev.criteria["min_success_rate"] == 0.9
    assert ev.tolerance == 3.0


def test_load_baseline_returns_false_on_missing_file() -> None:
    ev = BaselineEvaluator()
    assert ev.load_baseline(Path("/nonexistent/baseline.json")) is False
    assert ev.baseline_loaded is False


# ── evaluate_molecule ─────────────────────────────────────────────────────────


def test_evaluate_molecule_passes_when_hof_in_range() -> None:
    ev = BaselineEvaluator(tolerance=2.5)
    result = _pm7_result("m1", success=True, hof=-65.0, grade=QualityGrade.A)
    expected = {"hof_value": -65.0, "success_expected": True}

    mol_ev = ev.evaluate_molecule(result, expected)

    assert mol_ev.passed is True
    assert mol_ev.hof_in_range is True
    assert mol_ev.issues == []


def test_evaluate_molecule_fails_when_hof_out_of_range() -> None:
    ev = BaselineEvaluator(tolerance=2.5)
    result = _pm7_result("m1", success=True, hof=-70.0, grade=QualityGrade.A)
    expected = {"hof_value": -65.0}

    mol_ev = ev.evaluate_molecule(result, expected)

    assert mol_ev.passed is False
    assert mol_ev.hof_in_range is False
    assert any("hof_out_of_range" in issue for issue in mol_ev.issues)


def test_evaluate_molecule_fails_when_processing_failed() -> None:
    ev = BaselineEvaluator()
    result = _pm7_result("m1", success=False, hof=None, grade=QualityGrade.FAILED)

    mol_ev = ev.evaluate_molecule(result, {})

    assert mol_ev.passed is False
    assert "processing_failed" in mol_ev.issues


def test_evaluate_molecule_uses_baseline_lookup_when_expected_is_none() -> None:
    ev = BaselineEvaluator(tolerance=2.5)
    ev.baseline["m1"] = {"hof_value": -65.0}
    result = _pm7_result("m1", success=True, hof=-65.0)

    mol_ev = ev.evaluate_molecule(result, expected=None)

    assert mol_ev.hof_expected == -65.0
    assert mol_ev.hof_in_range is True


# ── evaluate_phase_a ──────────────────────────────────────────────────────────


def test_evaluate_phase_a_computes_baseline_pass_rate() -> None:
    ev = BaselineEvaluator(tolerance=2.5)
    ev.baseline = {
        "m1": {"hof_value": -65.0},
        "m2": {"hof_value": -65.0},
    }
    results = [
        _pm7_result("m1", success=True, hof=-65.0),
        _pm7_result("m2", success=True, hof=-70.0),  # out of range
    ]

    phase_ev = ev.evaluate_phase_a(results)

    assert phase_ev.baseline_pass_rate == pytest.approx(0.5)


def test_evaluate_phase_a_passes_when_all_criteria_met() -> None:
    ev = BaselineEvaluator(tolerance=2.5)
    ev.baseline = {"m1": {"hof_value": -65.0}}
    ev.criteria = {
        "min_success_rate": 1.0,
        "min_hof_extraction_rate": 1.0,
        "min_baseline_pass_rate": 1.0,
        "min_grade_ab_rate": 0.5,
    }
    results = [_pm7_result("m1", success=True, hof=-65.0)]

    phase_ev = ev.evaluate_phase_a(results)

    assert phase_ev.passed is True
    assert phase_ev.issues == []


def test_evaluate_phase_a_fails_when_criteria_not_met() -> None:
    ev = BaselineEvaluator(tolerance=2.5)
    ev.baseline = {"m1": {"hof_value": -65.0}}
    ev.criteria = {"min_success_rate": 1.0}
    results = [_pm7_result("m1", success=False, hof=None, grade=QualityGrade.FAILED)]

    phase_ev = ev.evaluate_phase_a(results)

    assert phase_ev.passed is False
    assert any("success_rate" in issue for issue in phase_ev.issues)


def test_evaluate_phase_a_empty_results_adds_issue() -> None:
    ev = BaselineEvaluator()
    phase_ev = ev.evaluate_phase_a([])
    assert "no_results" in phase_ev.issues


def test_phase_a_baseline_result_to_dict_round_trips() -> None:
    ev = BaselineEvaluator(tolerance=2.5)
    ev.baseline = {"m1": {"hof_value": -65.0}}
    results = [_pm7_result("m1", success=True, hof=-65.0)]

    phase_ev: PhaseABaselineResult = ev.evaluate_phase_a(results)
    d = phase_ev.to_dict()

    assert d["baseline_pass_rate"] == phase_ev.baseline_pass_rate
    assert len(d["molecules"]) == 1
    assert d["molecules"][0]["mol_id"] == "m1"
