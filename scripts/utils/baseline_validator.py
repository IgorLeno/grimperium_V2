#!/usr/bin/env python3
"""Baseline Validator for CREST-PM7 Pipeline.

Validates results against baseline expectations.

Usage:
    python scripts/utils/baseline_validator.py results.json baseline.json
"""

import argparse
import json
import sys
from pathlib import Path

# Add src to path for local development
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from grimperium.crest_pm7 import (
    PM7Result,
    ResultEvaluator,
    QualityGrade,
    TOLERANCE_ABSOLUTE,
)


def load_results(results_path: Path) -> list[PM7Result]:
    """Load results from JSON file.

    Args:
        results_path: Path to results JSON

    Returns:
        List of PM7Result objects
    """
    with open(results_path, encoding="utf-8") as f:
        data = json.load(f)

    results = []
    for r in data.get("results", []):
        # Reconstruct PM7Result from dict
        result = PM7Result(
            mol_id=r["mol_id"],
            smiles=r["smiles"],
        )
        result.success = r.get("success", False)
        result.quality_grade = QualityGrade(r.get("quality_grade", "FAILED"))

        # We need most_stable_hof but it's a property, so we need to hack this
        # For validation purposes, we'll use the value directly
        result._cached_hof = r.get("most_stable_hof")
        results.append(result)

    return results


def main() -> int:
    """Run baseline validation.

    Returns:
        Exit code (0 for pass, 1 for fail)
    """
    parser = argparse.ArgumentParser(description="Validate results against baseline")
    parser.add_argument("results", type=Path, help="Results JSON file")
    parser.add_argument("baseline", type=Path, help="Baseline JSON file")
    parser.add_argument(
        "--tolerance",
        type=float,
        default=TOLERANCE_ABSOLUTE,
        help=f"HOF tolerance in kcal/mol (default: {TOLERANCE_ABSOLUTE})",
    )
    args = parser.parse_args()

    if not args.results.exists():
        print(f"ERROR: Results file not found: {args.results}")
        return 1

    if not args.baseline.exists():
        print(f"ERROR: Baseline file not found: {args.baseline}")
        return 1

    print(f"Loading results from: {args.results}")
    print(f"Loading baseline from: {args.baseline}")
    print(f"Tolerance: {args.tolerance} kcal/mol")

    # Load results
    with open(args.results, encoding="utf-8") as f:
        results_data = json.load(f)

    # Load baseline
    with open(args.baseline, encoding="utf-8") as f:
        baseline = json.load(f)

    molecules_baseline = baseline.get("molecules", {})
    criteria = baseline.get("phase_a_success_criteria", {})

    print(f"\nBaseline contains {len(molecules_baseline)} molecules")
    print(f"Results contain {len(results_data.get('results', []))} molecules")

    # Validate each molecule
    print("\n--- Validation Results ---")
    passes = 0
    fails = 0

    for r in results_data.get("results", []):
        mol_id = r["mol_id"]
        expected = molecules_baseline.get(mol_id, {})

        if not expected:
            print(f"\n{mol_id}: SKIP (not in baseline)")
            continue

        hof_actual = r.get("most_stable_hof")
        hof_expected = expected.get("hof_value")
        hof_min = expected.get("hof_min", hof_expected - args.tolerance if hof_expected else None)
        hof_max = expected.get("hof_max", hof_expected + args.tolerance if hof_expected else None)

        success = r.get("success", False)
        grade = r.get("quality_grade", "FAILED")

        issues = []

        if not success:
            issues.append("processing_failed")

        if hof_actual is None:
            issues.append("hof_not_extracted")
        elif hof_min is not None and hof_max is not None:
            if not (hof_min <= hof_actual <= hof_max):
                issues.append(f"hof_out_of_range: {hof_actual:.2f} not in [{hof_min:.2f}, {hof_max:.2f}]")

        if grade not in ("A", "B"):
            issues.append(f"grade_not_acceptable: {grade}")

        passed = len(issues) == 0

        if passed:
            passes += 1
            status = "PASS"
        else:
            fails += 1
            status = "FAIL"

        print(f"\n{mol_id}: {status}")
        if hof_actual is not None and hof_expected is not None:
            print(f"  HOF: {hof_actual:.2f} kcal/mol (expected: {hof_expected:.2f})")
        print(f"  Grade: {grade}")
        if issues:
            for issue in issues:
                print(f"  Issue: {issue}")

    # Summary
    total = passes + fails
    pass_rate = passes / total if total > 0 else 0

    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"Total molecules: {total}")
    print(f"Passed: {passes}")
    print(f"Failed: {fails}")
    print(f"Pass rate: {pass_rate:.1%}")

    min_pass_rate = criteria.get("min_baseline_pass_rate", 1.0)
    if pass_rate >= min_pass_rate:
        print(f"\nVALIDATION: PASSED (>= {min_pass_rate:.1%})")
        return 0
    else:
        print(f"\nVALIDATION: FAILED (< {min_pass_rate:.1%})")
        return 1


if __name__ == "__main__":
    sys.exit(main())
