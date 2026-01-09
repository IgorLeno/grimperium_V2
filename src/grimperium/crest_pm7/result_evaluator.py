"""Result evaluation and baseline validation for CREST-PM7 Pipeline.

Evaluates results against expected baseline values.
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from .config import QualityGrade
from .molecule_processor import PM7Result

LOG = logging.getLogger("grimperium.crest_pm7.result_evaluator")

# Baseline tolerance
TOLERANCE_ABSOLUTE = 2.5  # kcal/mol


@dataclass
class MoleculeEvaluation:
    """Evaluation result for a single molecule.

    Attributes:
        mol_id: Molecule identifier
        passed: Whether the molecule passed all criteria
        hof_expected: Expected HOF value
        hof_actual: Actual HOF value
        hof_min: Minimum acceptable HOF
        hof_max: Maximum acceptable HOF
        hof_in_range: Whether HOF is within tolerance
        grade_expected: Expected quality grade
        grade_actual: Actual quality grade
        grade_acceptable: Whether grade is A or B
        success_expected: Expected success status
        success_actual: Actual success status
        issues: List of evaluation issues
    """

    mol_id: str
    passed: bool = False
    hof_expected: Optional[float] = None
    hof_actual: Optional[float] = None
    hof_min: Optional[float] = None
    hof_max: Optional[float] = None
    hof_in_range: bool = False
    grade_expected: Optional[str] = None
    grade_actual: Optional[str] = None
    grade_acceptable: bool = False
    success_expected: bool = True
    success_actual: bool = False
    issues: list[str] = field(default_factory=list)


@dataclass
class PhaseAEvaluation:
    """Evaluation result for Phase A.

    Attributes:
        passed: Whether Phase A passed all criteria
        molecules: Individual molecule evaluations
        success_rate: Overall success rate
        hof_extraction_rate: HOF extraction rate
        baseline_pass_rate: Baseline validation pass rate
        grade_ab_rate: Rate of A/B grades
        crash_count: Number of crashes
        criteria: Criteria used for evaluation
        issues: Overall evaluation issues
    """

    passed: bool = False
    molecules: list[MoleculeEvaluation] = field(default_factory=list)
    success_rate: float = 0.0
    hof_extraction_rate: float = 0.0
    baseline_pass_rate: float = 0.0
    grade_ab_rate: float = 0.0
    crash_count: int = 0
    criteria: dict = field(default_factory=dict)
    issues: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary with full diagnostic fields."""
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


class ResultEvaluator:
    """Evaluates PM7 results against baseline expectations."""

    def __init__(
        self,
        baseline_path: Optional[Path] = None,
        tolerance: float = TOLERANCE_ABSOLUTE,
    ) -> None:
        """Initialize evaluator.

        Args:
            baseline_path: Path to baseline JSON file
            tolerance: HOF tolerance in kcal/mol

        Attributes:
            baseline_loaded: True when load_baseline succeeded, False otherwise.
                Downstream callers should check baseline_loaded before relying on baseline validation.
        """
        self.tolerance = tolerance
        self.baseline: dict = {}
        self.criteria: dict = {}
        self.baseline_loaded: bool = False

        if baseline_path and baseline_path.exists():
            result = self.load_baseline(baseline_path)
            if not result:
                LOG.error(f"Failed to load baseline from {baseline_path}")

    def load_baseline(self, path: Path) -> bool:
        """Load baseline expectations from JSON file.

        Args:
            path: Path to baseline JSON

        Returns:
            True if loaded successfully
        """
        try:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)

            self.baseline = data.get("molecules", {})
            self.criteria = data.get("phase_a_success_criteria", {})
            self.tolerance = data.get("tolerance_kcal_mol", TOLERANCE_ABSOLUTE)
            self.baseline_loaded = True

            LOG.info(f"Loaded baseline with {len(self.baseline)} molecules")
            return True

        except Exception as e:
            LOG.warning(f"Failed to load baseline: {e}")
            self.baseline_loaded = False
            return False

    def evaluate_molecule(
        self,
        result: PM7Result,
        expected: Optional[dict] = None,
    ) -> MoleculeEvaluation:
        """Evaluate a single molecule result.

        Args:
            result: PM7Result to evaluate
            expected: Expected values (or None to use baseline)

        Returns:
            MoleculeEvaluation
        """
        if expected is None:
            expected = self.baseline.get(result.mol_id, {})

        eval_result = MoleculeEvaluation(mol_id=result.mol_id)

        # Success check
        eval_result.success_expected = expected.get("success_expected", True)
        eval_result.success_actual = result.success

        if eval_result.success_expected and not eval_result.success_actual:
            eval_result.issues.append("processing_failed")

        # HOF check
        hof_expected = expected.get("hof_value")
        if hof_expected is not None:
            eval_result.hof_expected = hof_expected
            eval_result.hof_min = expected.get("hof_min", hof_expected - self.tolerance)
            eval_result.hof_max = expected.get("hof_max", hof_expected + self.tolerance)
            eval_result.hof_actual = result.most_stable_hof

            if eval_result.hof_actual is not None:
                eval_result.hof_in_range = (
                    eval_result.hof_min <= eval_result.hof_actual <= eval_result.hof_max
                )
                if not eval_result.hof_in_range:
                    eval_result.issues.append(
                        f"hof_out_of_range: {eval_result.hof_actual:.2f} "
                        f"not in [{eval_result.hof_min:.2f}, {eval_result.hof_max:.2f}]"
                    )
            else:
                eval_result.issues.append("hof_not_extracted")

        # Grade check
        eval_result.grade_expected = expected.get("quality_grade_expected")
        eval_result.grade_actual = result.quality_grade.value
        eval_result.grade_acceptable = result.quality_grade in (QualityGrade.A, QualityGrade.B)

        if not eval_result.grade_acceptable:
            eval_result.issues.append(f"grade_not_acceptable: {eval_result.grade_actual}")

        # Determine pass/fail
        eval_result.passed = (
            eval_result.success_actual
            and (eval_result.hof_actual is not None)
            and eval_result.hof_in_range
            and eval_result.grade_acceptable
        )

        return eval_result

    def evaluate_phase_a(
        self,
        results: list[PM7Result],
    ) -> PhaseAEvaluation:
        """Evaluate Phase A results.

        Args:
            results: List of PM7Result objects

        Returns:
            PhaseAEvaluation with overall assessment
        """
        eval_result = PhaseAEvaluation()
        eval_result.criteria = self.criteria

        if not results:
            eval_result.issues.append("no_results")
            return eval_result

        # Evaluate each molecule
        for result in results:
            mol_eval = self.evaluate_molecule(result)
            eval_result.molecules.append(mol_eval)

        n_total = len(results)

        # Calculate rates
        n_success = sum(1 for r in results if r.success)
        eval_result.success_rate = n_success / n_total

        n_hof = sum(1 for r in results if r.most_stable_hof is not None)
        eval_result.hof_extraction_rate = n_hof / n_total

        n_baseline_pass = sum(1 for m in eval_result.molecules if m.passed)
        eval_result.baseline_pass_rate = n_baseline_pass / n_total

        n_grade_ab = sum(
            1 for r in results
            if r.quality_grade in (QualityGrade.A, QualityGrade.B)
        )
        eval_result.grade_ab_rate = n_grade_ab / n_total

        # Check criteria
        min_success = self.criteria.get("min_success_rate", 1.0)
        min_hof = self.criteria.get("min_hof_extraction_rate", 1.0)
        min_baseline = self.criteria.get("min_baseline_pass_rate", 1.0)
        min_grade_ab = self.criteria.get("min_grade_ab_rate", 0.67)
        require_zero_crashes = self.criteria.get("require_zero_crashes", True)

        if eval_result.success_rate < min_success:
            eval_result.issues.append(
                f"success_rate {eval_result.success_rate:.1%} < {min_success:.1%}"
            )

        if eval_result.hof_extraction_rate < min_hof:
            eval_result.issues.append(
                f"hof_extraction_rate {eval_result.hof_extraction_rate:.1%} < {min_hof:.1%}"
            )

        if eval_result.baseline_pass_rate < min_baseline:
            eval_result.issues.append(
                f"baseline_pass_rate {eval_result.baseline_pass_rate:.1%} < {min_baseline:.1%}"
            )

        if eval_result.grade_ab_rate < min_grade_ab:
            eval_result.issues.append(
                f"grade_ab_rate {eval_result.grade_ab_rate:.1%} < {min_grade_ab:.1%}"
            )

        if require_zero_crashes and eval_result.crash_count > 0:
            eval_result.issues.append(f"crashes: {eval_result.crash_count}")

        # Overall pass/fail
        eval_result.passed = len(eval_result.issues) == 0

        LOG.info(
            f"Phase A evaluation: passed={eval_result.passed}, "
            f"success={eval_result.success_rate:.1%}, "
            f"baseline={eval_result.baseline_pass_rate:.1%}"
        )

        return eval_result
