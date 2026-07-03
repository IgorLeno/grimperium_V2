"""Operational result evaluation for the CREST-PM7 pipeline.

Checks whether molecules succeeded, whether HOF was extracted, and
whether the quality grade is acceptable (A or B).  No comparison
against expected reference values — that logic lives in
tests/regression/baseline_evaluation.py.
"""

import logging
from dataclasses import dataclass, field

from .config import QualityGrade
from .molecule_processor import PM7Result

LOG = logging.getLogger("grimperium.crest_pm7.result_evaluator")


@dataclass
class MoleculeEvaluation:
    """Operational evaluation for a single molecule.

    Attributes:
        mol_id: Molecule identifier.
        success_actual: Whether processing succeeded.
        hof_actual: Extracted HOF value, or None if not available.
        grade_actual: Quality grade string.
        grade_acceptable: True when grade is A or B.
        issues: Short issue tags for rapid filtering.
    """

    mol_id: str
    success_actual: bool = False
    hof_actual: float | None = None
    grade_actual: str | None = None
    grade_acceptable: bool = False
    issues: list[str] = field(default_factory=list)


@dataclass
class PhaseAEvaluation:
    """Operational metrics for a Phase A batch.

    Rates are descriptive — no pass/fail threshold without a reference
    baseline.  Use BaselineEvaluator in tests/regression/ for that.

    Attributes:
        molecules: Per-molecule operational evaluations.
        success_rate: Fraction of molecules that succeeded.
        hof_extraction_rate: Fraction with a valid HOF value.
        grade_ab_rate: Fraction graded A or B.
        crash_count: Molecules with fatal errors (currently unused; reserved).
        issues: Batch-level issues, e.g. "no_results".
    """

    molecules: list[MoleculeEvaluation] = field(default_factory=list)
    success_rate: float = 0.0
    hof_extraction_rate: float = 0.0
    grade_ab_rate: float = 0.0
    crash_count: int = 0
    issues: list[str] = field(default_factory=list)


class ResultEvaluator:
    """Evaluates PM7 results operationally (no reference baseline)."""

    def evaluate_molecule(self, result: PM7Result) -> MoleculeEvaluation:
        """Evaluate a single molecule operationally.

        Args:
            result: PM7Result from the pipeline.

        Returns:
            MoleculeEvaluation with operational fields.
        """
        ev = MoleculeEvaluation(mol_id=result.mol_id)
        ev.success_actual = result.success
        ev.hof_actual = result.most_stable_hof
        ev.grade_actual = result.quality_grade.value
        ev.grade_acceptable = result.quality_grade in (QualityGrade.A, QualityGrade.B)

        if not ev.success_actual:
            ev.issues.append("processing_failed")
        if ev.hof_actual is None:
            ev.issues.append("hof_not_extracted")
        if not ev.grade_acceptable:
            ev.issues.append(f"grade_not_acceptable: {ev.grade_actual}")

        return ev

    def evaluate_phase_a(self, results: list[PM7Result]) -> PhaseAEvaluation:
        """Compute operational metrics for a batch of results.

        Args:
            results: PM7Result objects from a processing run.

        Returns:
            PhaseAEvaluation with aggregate rates.
        """
        ev = PhaseAEvaluation()

        if not results:
            ev.issues.append("no_results")
            return ev

        for result in results:
            ev.molecules.append(self.evaluate_molecule(result))

        n_total = len(results)
        ev.success_rate = sum(1 for r in results if r.success) / n_total
        ev.hof_extraction_rate = (
            sum(1 for r in results if r.most_stable_hof is not None) / n_total
        )
        ev.grade_ab_rate = (
            sum(
                1
                for r in results
                if r.quality_grade in (QualityGrade.A, QualityGrade.B)
            )
            / n_total
        )

        LOG.info(
            "Phase A metrics: success=%.1f%%, hof_extracted=%.1f%%, grade_ab=%.1f%%",
            ev.success_rate * 100,
            ev.hof_extraction_rate * 100,
            ev.grade_ab_rate * 100,
        )

        return ev
