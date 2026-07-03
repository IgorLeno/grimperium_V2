"""Baseline comparison helper for CREST-PM7 regression testing.

This module contains the reference-comparison logic that was removed from
the operational core (src/grimperium/crest_pm7/result_evaluator.py) as part
of the PR7 boundary cleanup.  It compares pipeline output against stored
expected HOF values and quality grades.

Intended use: call BaselineEvaluator inside pytest regression tests that
load a baseline JSON fixture and assert pass rates.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from grimperium import DictStrAny
from grimperium.crest_pm7.config import QualityGrade
from grimperium.crest_pm7.molecule_processor import PM7Result

LOG = logging.getLogger("grimperium.regression.baseline_evaluation")

# Default HOF tolerance used when the baseline JSON does not specify one.
TOLERANCE_ABSOLUTE: float = 2.5  # kcal/mol


@dataclass
class MoleculeBaselineResult:
    """Reference-comparison result for a single molecule.

    Attributes:
        mol_id: Molecule identifier.
        passed: True when all reference checks pass.
        hof_expected: Expected HOF from baseline.
        hof_actual: Actual extracted HOF.
        hof_min: Lower HOF bound (expected − tolerance).
        hof_max: Upper HOF bound (expected + tolerance).
        hof_in_range: Whether actual HOF is within [hof_min, hof_max].
        grade_expected: Grade recorded in baseline.
        grade_actual: Grade produced by the pipeline.
        grade_acceptable: True when grade is A or B.
        success_expected: Whether this molecule was expected to succeed.
        success_actual: Whether it actually succeeded.
        issues: Short issue tags.
    """

    mol_id: str
    passed: bool = False
    hof_expected: float | None = None
    hof_actual: float | None = None
    hof_min: float | None = None
    hof_max: float | None = None
    hof_in_range: bool = False
    grade_expected: str | None = None
    grade_actual: str | None = None
    grade_acceptable: bool = False
    success_expected: bool = True
    success_actual: bool = False
    issues: list[str] = field(default_factory=list)


@dataclass
class PhaseABaselineResult:
    """Reference-comparison result for a Phase A batch.

    Attributes:
        passed: True when all criteria thresholds are met.
        molecules: Per-molecule baseline results.
        success_rate: Fraction of molecules that succeeded.
        hof_extraction_rate: Fraction with a valid HOF value.
        baseline_pass_rate: Fraction that passed all reference checks.
        grade_ab_rate: Fraction graded A or B.
        crash_count: Molecules with fatal errors (reserved).
        criteria: Criteria dict loaded from the baseline JSON.
        issues: Threshold violations and other batch issues.
    """

    passed: bool = False
    molecules: list[MoleculeBaselineResult] = field(default_factory=list)
    success_rate: float = 0.0
    hof_extraction_rate: float = 0.0
    baseline_pass_rate: float = 0.0
    grade_ab_rate: float = 0.0
    crash_count: int = 0
    criteria: DictStrAny = field(default_factory=dict)
    issues: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a dictionary for logging or assertion output."""
        return {
            "passed": self.passed,
            "molecules": [
                {
                    "mol_id": m.mol_id,
                    "passed": m.passed,
                    "hof_expected": m.hof_expected,
                    "hof_actual": m.hof_actual,
                    "hof_min": m.hof_min,
                    "hof_max": m.hof_max,
                    "hof_in_range": m.hof_in_range,
                    "grade_expected": m.grade_expected,
                    "grade_actual": m.grade_actual,
                    "grade_acceptable": m.grade_acceptable,
                    "success_expected": m.success_expected,
                    "success_actual": m.success_actual,
                    "issues": m.issues,
                }
                for m in self.molecules
            ],
            "success_rate": self.success_rate,
            "hof_extraction_rate": self.hof_extraction_rate,
            "baseline_pass_rate": self.baseline_pass_rate,
            "grade_ab_rate": self.grade_ab_rate,
            "crash_count": self.crash_count,
            "criteria": self.criteria,
            "issues": self.issues,
        }


class BaselineEvaluator:
    """Compares PM7 results against stored baseline expectations.

    Load a baseline JSON with :meth:`load_baseline`, then call
    :meth:`evaluate_phase_a` to get a :class:`PhaseABaselineResult`.

    Baseline JSON shape::

        {
          "molecules": {
            "mol_id": {
              "hof_value": -65.0,        # optional
              "hof_min": -67.5,          # optional, overrides tolerance
              "hof_max": -62.5,          # optional, overrides tolerance
              "quality_grade_expected": "A",  # optional
              "success_expected": true   # optional, default true
            }
          },
          "phase_a_success_criteria": {
            "min_success_rate": 1.0,
            "min_hof_extraction_rate": 1.0,
            "min_baseline_pass_rate": 1.0,
            "min_grade_ab_rate": 0.67,
            "require_zero_crashes": true
          },
          "tolerance_kcal_mol": 2.5     # optional, overrides TOLERANCE_ABSOLUTE
        }
    """

    def __init__(self, tolerance: float = TOLERANCE_ABSOLUTE) -> None:
        """Initialise with optional tolerance override.

        Args:
            tolerance: HOF tolerance in kcal/mol used when the baseline
                JSON does not specify hof_min/hof_max per molecule.
        """
        self.tolerance = tolerance
        self.baseline: DictStrAny = {}
        self.criteria: DictStrAny = {}
        self.baseline_loaded: bool = False

    def load_baseline(self, path: Path) -> bool:
        """Load baseline expectations from a JSON file.

        Args:
            path: Path to baseline JSON.

        Returns:
            True on success; False if the file could not be read.
        """
        try:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)

            self.baseline = data.get("molecules", {})
            self.criteria = data.get("phase_a_success_criteria", {})
            self.tolerance = data.get("tolerance_kcal_mol", self.tolerance)
            self.baseline_loaded = True

            LOG.info("Loaded baseline with %d molecules", len(self.baseline))
            return True
        except Exception as exc:
            LOG.warning("Failed to load baseline: %s", exc)
            self.baseline_loaded = False
            return False

    def evaluate_molecule(
        self,
        result: PM7Result,
        expected: DictStrAny | None = None,
    ) -> MoleculeBaselineResult:
        """Compare a single molecule result against the baseline.

        Args:
            result: PM7Result to evaluate.
            expected: Expected-value dict; falls back to ``self.baseline``
                lookup when None.

        Returns:
            :class:`MoleculeBaselineResult`.
        """
        if expected is None:
            expected = self.baseline.get(result.mol_id, {})

        ev = MoleculeBaselineResult(mol_id=result.mol_id)

        # Success
        ev.success_expected = expected.get("success_expected", True)
        ev.success_actual = result.success
        if ev.success_expected and not ev.success_actual:
            ev.issues.append("processing_failed")

        # HOF comparison
        hof_expected = expected.get("hof_value")
        if hof_expected is not None:
            ev.hof_expected = hof_expected
            ev.hof_min = expected.get("hof_min", hof_expected - self.tolerance)
            ev.hof_max = expected.get("hof_max", hof_expected + self.tolerance)
            ev.hof_actual = result.most_stable_hof

            if ev.hof_actual is not None:
                ev.hof_in_range = ev.hof_min <= ev.hof_actual <= ev.hof_max
                if not ev.hof_in_range:
                    ev.issues.append(
                        f"hof_out_of_range: {ev.hof_actual:.2f} "
                        f"not in [{ev.hof_min:.2f}, {ev.hof_max:.2f}]"
                    )
            else:
                ev.issues.append("hof_not_extracted")

        # Grade
        ev.grade_expected = expected.get("quality_grade_expected")
        ev.grade_actual = result.quality_grade.value
        ev.grade_acceptable = result.quality_grade in (QualityGrade.A, QualityGrade.B)
        if not ev.grade_acceptable:
            ev.issues.append(f"grade_not_acceptable: {ev.grade_actual}")

        # Overall pass
        ev.passed = (
            ev.success_actual
            and (ev.hof_actual is not None)
            and ev.hof_in_range
            and ev.grade_acceptable
        )

        return ev

    def evaluate_phase_a(self, results: list[PM7Result]) -> PhaseABaselineResult:
        """Run baseline comparison over a full batch.

        Args:
            results: PM7Result objects from a processing run.

        Returns:
            :class:`PhaseABaselineResult` with pass/fail verdict.
        """
        ev = PhaseABaselineResult(criteria=self.criteria)

        if not results:
            ev.issues.append("no_results")
            return ev

        for result in results:
            ev.molecules.append(self.evaluate_molecule(result))

        n_total = len(results)
        n_success = sum(1 for r in results if r.success)
        ev.success_rate = n_success / n_total

        n_hof = sum(1 for r in results if r.most_stable_hof is not None)
        ev.hof_extraction_rate = n_hof / n_total

        n_passed = sum(1 for m in ev.molecules if m.passed)
        ev.baseline_pass_rate = n_passed / n_total

        n_grade_ab = sum(
            1 for r in results if r.quality_grade in (QualityGrade.A, QualityGrade.B)
        )
        ev.grade_ab_rate = n_grade_ab / n_total

        # Threshold checks from the loaded criteria
        min_success = self.criteria.get("min_success_rate", 1.0)
        min_hof = self.criteria.get("min_hof_extraction_rate", 1.0)
        min_baseline = self.criteria.get("min_baseline_pass_rate", 1.0)
        min_grade_ab = self.criteria.get("min_grade_ab_rate", 0.67)
        require_zero_crashes = self.criteria.get("require_zero_crashes", True)

        if ev.success_rate < min_success:
            ev.issues.append(f"success_rate {ev.success_rate:.1%} < {min_success:.1%}")
        if ev.hof_extraction_rate < min_hof:
            ev.issues.append(
                f"hof_extraction_rate {ev.hof_extraction_rate:.1%} < {min_hof:.1%}"
            )
        if ev.baseline_pass_rate < min_baseline:
            ev.issues.append(
                f"baseline_pass_rate {ev.baseline_pass_rate:.1%} < {min_baseline:.1%}"
            )
        if ev.grade_ab_rate < min_grade_ab:
            ev.issues.append(
                f"grade_ab_rate {ev.grade_ab_rate:.1%} < {min_grade_ab:.1%}"
            )
        if require_zero_crashes and ev.crash_count > 0:
            ev.issues.append(f"crashes: {ev.crash_count}")

        ev.passed = len(ev.issues) == 0

        LOG.info(
            "Baseline evaluation: passed=%s, success=%.1f%%, baseline=%.1f%%",
            ev.passed,
            ev.success_rate * 100,
            ev.baseline_pass_rate * 100,
        )

        return ev
