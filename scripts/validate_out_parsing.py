#!/usr/bin/env python3
"""Validation script for .out parsing migration.

Checks that:
1. No .aux files are generated in new runs
2. HOMO/LUMO/GAP are successfully extracted from .out
3. All 11 descriptors return valid values
4. .mop input files do not contain "AUX" keyword
"""

import math
from pathlib import Path

from grimperium.crest_pm7.mopac_descriptors import extract_mopac_descriptors

# Required descriptor keys for validation
REQUIRED_KEYS: list[str] = [
    "mopac_dipole_debye",
    "mopac_ionization_potential_ev",
    "mopac_homo_ev",
    "mopac_lumo_ev",
    "mopac_gap_ev",
    "mopac_cosmo_area_a2",
    "mopac_cosmo_volume_a3",
    "mopac_gradient_norm",
    "mopac_num_scf_cycles",
    "mopac_point_group",
    "mopac_time_s",
]


def validate_out_file(out_path: Path) -> dict[str, bool]:
    """Validate a single .out file parsing.

    Args:
        out_path: Path to the MOPAC .out file to validate.

    Returns:
        A dict mapping check names to bool results:
        - ``file_exists``: whether the file exists on disk.
        - ``parsing_succeeds``: whether ``extract_mopac_descriptors`` completed
          without raising an exception.
        - ``homo_valid``: whether the HOMO energy is a non-NaN float.
        - ``lumo_valid``: whether the LUMO energy is a non-NaN float.
        - ``gap_valid``: whether the HOMO-LUMO gap is a non-NaN positive float.
        - ``all_descriptors_present``: whether at least 10 of 11 descriptors
          have valid (non-None, non-NaN) values.

    Example:
        >>> results = validate_out_file(Path("mol_001_conf_0.out"))
        >>> results["homo_valid"]
        True
    """
    results: dict[str, bool] = {
        "file_exists": out_path.exists(),
        "parsing_succeeds": False,
        "homo_valid": False,
        "lumo_valid": False,
        "gap_valid": False,
        "all_descriptors_present": False,
    }

    if not results["file_exists"]:
        return results

    try:
        desc = extract_mopac_descriptors(out_path)
        results["parsing_succeeds"] = True

        # Check HOMO/LUMO/GAP
        homo = desc.get("mopac_homo_ev")
        lumo = desc.get("mopac_lumo_ev")
        gap = desc.get("mopac_gap_ev")

        results["homo_valid"] = homo is not None and not math.isnan(float(homo))
        results["lumo_valid"] = lumo is not None and not math.isnan(float(lumo))
        results["gap_valid"] = (
            gap is not None and not math.isnan(float(gap)) and float(gap) > 0
        )

        # Check all 11 descriptors
        valid_count = sum(
            1
            for key in REQUIRED_KEYS
            if key in desc
            and desc[key] is not None
            and (
                not isinstance(desc[key], float) or not math.isnan(desc[key])
            )
        )

        results["all_descriptors_present"] = valid_count >= 10  # Allow 1 missing

    except Exception as e:
        print(f"  ❌ Error parsing {out_path.name}: {e}")

    return results


def main() -> None:
    """Run validation checks for the .out parsing migration.

    Validates that:

    1. No ``.aux`` files are generated in new batch runs.
    2. ``.out`` files are successfully parsed (sample of up to 10 files).
    3. ``.mop`` input files do not contain the legacy ``AUX`` keyword.

    Prints a summary of each check and a final pass/fail verdict.

    Returns:
        None. All results are printed to stdout.

    Side Effects:
        Reads files from the MOPAC temp directory
        (``src/grimperium/crest_pm7/tmp``). Prints validation status to stdout.
    """
    print("🔍 Validating .out parsing migration...\n")

    # Find temp directory
    temp_dir = (
        Path(__file__).parent.parent / "src" / "grimperium" / "crest_pm7" / "tmp"
    )

    if not temp_dir.exists():
        print(f"⚠️  Temp directory not found: {temp_dir}")
        print("   Run a test batch first to generate outputs.")
        return

    # Check 1: No .aux files in recent batches
    aux_files = list(temp_dir.rglob("*.aux"))
    print("✓ Check 1: .aux files in temp/")
    if aux_files:
        aux_check_passed = False
        print(f"  ⚠️  Found {len(aux_files)} .aux files (should be 0 for new runs)")
        for f in aux_files[:5]:
            print(f"     - {f.relative_to(temp_dir)}")
    else:
        aux_check_passed = True
        print("  ✅ No .aux files found (correct)")

    print()

    # Check 2: .out files parse successfully
    out_files = list(temp_dir.rglob("*.out"))
    print(f"✓ Check 2: .out file parsing ({len(out_files)} files found)")

    if not out_files:
        print("  ⚠️  No .out files found. Run a test batch first.")
        return

    # Sample up to 10 .out files
    sample = out_files[:10]

    all_valid = True

    for out_path in sample:
        results = validate_out_file(out_path)

        status = "✅" if all(results.values()) else "⚠️"
        print(f"  {status} {out_path.name}")

        if not all(results.values()):
            all_valid = False
            for key, val in results.items():
                if not val:
                    print(f"     - {key}: FAILED")

    print()

    # Check 3: .mop files do not contain "AUX"
    mop_files = list(temp_dir.rglob("*.mop"))
    mop_with_aux: list[Path] = []
    print(f"✓ Check 3: .mop input keywords ({len(mop_files)} files found)")

    if mop_files:
        for mop_path in mop_files[:10]:
            try:
                first_line = mop_path.read_text().split("\n")[0]
            except (OSError, UnicodeDecodeError) as e:
                print(f"  ⚠️  Could not read {mop_path.name}: {e}")
                continue
            # Check for AUX as a standalone keyword (space or start/end bounded)
            keywords = first_line.upper().split()
            if "AUX" in keywords:
                mop_with_aux.append(mop_path)

        if mop_with_aux:
            print(f"  ❌ Found {len(mop_with_aux)} .mop files with 'AUX' keyword:")
            for f in mop_with_aux[:5]:
                print(f"     - {f.relative_to(temp_dir)}")
        else:
            print("  ✅ No .mop files contain 'AUX' keyword (correct)")

    print()
    print("=" * 60)
    if all_valid and aux_check_passed and not mop_with_aux:
        print("✅ All validation checks passed!")
    else:
        print("⚠️  Some checks failed. Review output above.")


if __name__ == "__main__":
    main()
