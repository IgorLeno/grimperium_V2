#!/usr/bin/env python3
"""Phase A Quick Test for CREST-PM7 Pipeline.

Runs the pipeline on baseline molecules and validates results.

Usage:
    python scripts/phase_a_quick_test.py

Exit codes:
    0: All Phase A criteria passed
    1: Phase A criteria failed
"""

import csv
import json
import sys
from pathlib import Path

# Add src to path for local development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from grimperium.crest_pm7 import (
    CRESTPM7Pipeline,
    PM7Config,
    validate_environment,
)


def load_test_molecules(csv_path: Path) -> list[tuple[str, str]]:
    """Load test molecules from CSV file.

    Args:
        csv_path: Path to CSV file

    Returns:
        List of (mol_id, smiles) tuples
    """
    molecules = []
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            molecules.append((row["mol_id"], row["smiles"]))
    return molecules


def main() -> int:
    """Run Phase A quick test.

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    print("=" * 60)
    print("CREST-PM7 Pipeline - Phase A Quick Test")
    print("=" * 60)

    # Paths
    base_dir = Path(__file__).parent.parent
    molecules_csv = base_dir / "data/molecules_pm7/testing/baselines/phase_a_molecules.csv"
    expected_json = base_dir / "data/molecules_pm7/testing/baselines/phase_a_expected.json"
    output_dir = base_dir / "data/molecules_pm7/computed"
    results_json = output_dir / "phase_a_results.json"

    # Check files exist
    if not molecules_csv.exists():
        print(f"ERROR: Molecules CSV not found: {molecules_csv}")
        return 1

    if not expected_json.exists():
        print(f"ERROR: Expected JSON not found: {expected_json}")
        return 1

    # Load test molecules
    molecules = load_test_molecules(molecules_csv)
    print(f"\nLoaded {len(molecules)} test molecules:")
    for mol_id, smiles in molecules:
        print(f"  - {mol_id}: {smiles}")

    # Configure pipeline
    config = PM7Config(
        phase="A",
        max_conformers=3,
        output_dir=output_dir,
    )

    pipeline = CRESTPM7Pipeline(config)

    # Validate environment
    print("\n--- Environment Validation ---")
    if not pipeline.validate():
        print("ERROR: Environment validation failed")
        return 1
    print("Environment OK")

    # Setup pipeline
    print("\n--- Pipeline Setup ---")
    pipeline.setup(session_name="phase_a_quick_test")
    pipeline.load_baseline(expected_json)

    # Process molecules
    print("\n--- Processing Molecules ---")
    for mol_id, smiles in molecules:
        print(f"\nProcessing: {mol_id}")
        result = pipeline.process_molecule(mol_id, smiles)
        print(f"  Success: {result.success}")
        print(f"  Grade: {result.quality_grade.value}")
        print(f"  HOF: {result.most_stable_hof}")
        if result.error_message:
            print(f"  Error: {result.error_message}")

    # Save results
    print("\n--- Saving Results ---")
    pipeline.save_results(results_json)
    print(f"Results saved to: {results_json}")

    # Evaluate Phase A
    print("\n--- Phase A Evaluation ---")
    evaluation = pipeline.evaluate_phase_a()

    print(f"\nSuccess Rate: {evaluation.success_rate:.1%}")
    print(f"HOF Extraction Rate: {evaluation.hof_extraction_rate:.1%}")
    print(f"Baseline Pass Rate: {evaluation.baseline_pass_rate:.1%}")
    print(f"Grade A/B Rate: {evaluation.grade_ab_rate:.1%}")

    print("\nMolecule Results:")
    for mol_eval in evaluation.molecules:
        status = "PASS" if mol_eval.passed else "FAIL"
        print(f"  {mol_eval.mol_id}: {status}")
        if mol_eval.hof_actual is not None:
            print(f"    HOF: {mol_eval.hof_actual:.2f} (expected: {mol_eval.hof_expected:.2f})")
        if mol_eval.issues:
            for issue in mol_eval.issues:
                print(f"    Issue: {issue}")

    print("\n" + "=" * 60)
    if evaluation.passed:
        print("PHASE A: PASSED")
        print("=" * 60)
        return 0
    else:
        print("PHASE A: FAILED")
        print("Issues:")
        for issue in evaluation.issues:
            print(f"  - {issue}")
        print("=" * 60)
        return 1


if __name__ == "__main__":
    sys.exit(main())
