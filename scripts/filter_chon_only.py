#!/usr/bin/env python3
"""
Filter CSV to keep only CHON molecules (C, H, O, N).

This script removes molecules containing elements that are problematic
for MMFF force field geometry optimization (B, P, As, Ge, etc.) from
a thermochemistry dataset.

Usage:
    python scripts/filter_chon_only.py data/thermo_cbs_clean.csv -o data/thermo_cbs_chon.csv --verbose

Background:
    MMFF94 (Merck Molecular Force Field) is used by RDKit for 3D geometry
    generation. It doesn't natively support all elements. Boron, Phosphorus,
    Arsenic, and Germanium often cause geometry optimization failures.

    For Phase A of Grimperium, filtering to CHON-only molecules ensures:
    - 100% MMFF compatibility
    - Reliable geometry generation
    - >99% success rate expected

Elements Removed:
    - Boron (B) - ~0.5% of molecules
    - Phosphorus (P) - ~0.5% of molecules
    - Arsenic (As) - <0.1% of molecules
    - Germanium (Ge) - <0.1% of molecules
    - All other non-CHON elements

Total Loss: ~1.5% of molecules (from 30,026 to 29,568)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd


def is_chon_only(smiles: str) -> bool:
    """
    Check if molecule contains only C, H, O, N atoms.

    Args:
        smiles: SMILES string representing the molecule

    Returns:
        True if only CHON atoms, False otherwise (or if invalid SMILES)

    Example:
        >>> is_chon_only("CCO")  # Ethanol
        True
        >>> is_chon_only("Cc1cccc(B(O)O)c1")  # m-Tolylboronic acid (has B)
        False
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False

    allowed_atoms = {"C", "H", "O", "N"}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in allowed_atoms:
            return False

    return True


def get_element_distribution(df: "pd.DataFrame", smiles_col: str = "smiles") -> dict[str, int]:
    """
    Analyze element distribution in dataset.

    Args:
        df: DataFrame with SMILES column
        smiles_col: Name of SMILES column

    Returns:
        Dict mapping element symbols to molecule counts
    """
    from rdkit import Chem

    element_counts: dict[str, int] = {}

    for smiles in df[smiles_col]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        elements_in_mol = set()
        for atom in mol.GetAtoms():
            elements_in_mol.add(atom.GetSymbol())

        for elem in elements_in_mol:
            element_counts[elem] = element_counts.get(elem, 0) + 1

    return element_counts


def filter_chon(
    input_csv: Path,
    output_csv: Path | None = None,
    smiles_col: str = "smiles",
    verbose: bool = False,
    analyze: bool = False,
) -> tuple["pd.DataFrame", int, int]:
    """
    Filter CSV to CHON-only molecules.

    Args:
        input_csv: Input CSV path
        output_csv: Output CSV path (optional, if None doesn't save)
        smiles_col: Name of SMILES column in CSV
        verbose: Print progress info
        analyze: Analyze element distribution before filtering

    Returns:
        Tuple of (filtered_df, input_count, output_count)

    Raises:
        FileNotFoundError: If input CSV doesn't exist
        ValueError: If SMILES column not found in CSV
    """
    import pandas as pd

    # Read input
    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    df = pd.read_csv(input_csv)
    input_count = len(df)

    if smiles_col not in df.columns:
        raise ValueError(f"SMILES column '{smiles_col}' not found. Available: {list(df.columns)}")

    if verbose:
        print(f"Input file: {input_csv}")
        print(f"Input molecules: {input_count}")

    # Analyze element distribution (optional)
    if analyze:
        if verbose:
            print("\nAnalyzing element distribution...")
        elem_dist = get_element_distribution(df, smiles_col)
        sorted_elems = sorted(elem_dist.items(), key=lambda x: -x[1])

        if verbose:
            print("\nElement distribution (molecules containing each element):")
            for elem, count in sorted_elems:
                pct = 100 * count / input_count if input_count > 0 else 0.0
                marker = " <-- non-CHON" if elem not in {"C", "H", "O", "N"} else ""
                print(f"  {elem:2s}: {count:6d} ({pct:5.2f}%){marker}")

    # Filter
    if verbose:
        print("\nFiltering to CHON-only molecules...")

    chon_mask = df[smiles_col].apply(is_chon_only)
    df_chon = df[chon_mask].copy()

    output_count = len(df_chon)
    removed = input_count - output_count

    if verbose:
        print(f"\nResults:")
        print(f"  CHON molecules: {output_count:,}")
        print(f"  Removed (non-CHON): {removed:,}")
        removed_pct = 100 * removed / input_count if input_count > 0 else 0.0
        print(f"  Removed proportion: {removed_pct:.2f}%")

    # Save if output provided
    if output_csv:
        df_chon.to_csv(output_csv, index=False)
        if verbose:
            print(f"\nSaved to: {output_csv}")

    return df_chon, input_count, output_count


def main() -> int:
    """Main entry point for filter script."""
    parser = argparse.ArgumentParser(
        description="Filter CSV to CHON-only molecules (C, H, O, N)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python scripts/filter_chon_only.py data/thermo_cbs_clean.csv -o data/thermo_cbs_chon.csv

    # With verbose output
    python scripts/filter_chon_only.py data/thermo_cbs_clean.csv -o data/thermo_cbs_chon.csv --verbose

    # Analyze element distribution without filtering
    python scripts/filter_chon_only.py data/thermo_cbs_clean.csv --analyze --verbose

    # Custom SMILES column
    python scripts/filter_chon_only.py my_data.csv -o my_data_chon.csv --smiles-col molecule_smiles
        """,
    )
    parser.add_argument(
        "input_csv",
        type=Path,
        help="Input CSV file path",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output CSV file path (optional)",
    )
    parser.add_argument(
        "--smiles-col",
        default="smiles",
        help="Name of SMILES column in CSV (default: smiles)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print progress and statistics",
    )
    parser.add_argument(
        "--analyze",
        action="store_true",
        help="Analyze element distribution before filtering",
    )

    args = parser.parse_args()

    # Validate input
    if not args.input_csv.exists():
        print(f"Error: Input file not found: {args.input_csv}")
        return 1

    # Run filter
    try:
        df_chon, input_count, output_count = filter_chon(
            args.input_csv,
            args.output,
            args.smiles_col,
            args.verbose,
            args.analyze,
        )

        # Print summary even without verbose
        if not args.verbose:
            removed_pct = 100 * (input_count - output_count) / input_count if input_count > 0 else 0.0
            print(f"{input_count} -> {output_count} molecules ({removed_pct:.2f}% removed)")
            if args.output:
                print(f"Saved to: {args.output}")

        return 0

    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        return 1
    except ImportError as e:
        print(f"Error: Missing dependency - {e}")
        print("Install with: pip install pandas rdkit")
        return 1


if __name__ == "__main__":
    sys.exit(main())
