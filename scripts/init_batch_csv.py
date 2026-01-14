#!/usr/bin/env python3
"""Initialize batch tracking CSV from source data.

This script creates a batch tracking CSV file from thermo_cbs_clean.csv,
adding all necessary columns for batch processing status tracking.

Usage:
    python scripts/init_batch_csv.py --input data/thermo_cbs_clean.csv \
        --output data/batch_tracking.csv
    python scripts/init_batch_csv.py --input data/thermo_cbs_clean.csv \
        --output data/test_batch.csv --limit 10

The output CSV will have all columns defined in BatchRowCSV model.
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
LOG = logging.getLogger(__name__)


def compute_descriptors(smiles: str) -> dict:
    """Compute molecular descriptors from SMILES.

    Args:
        smiles: SMILES string

    Returns:
        Dict with nrotbonds, tpsa, aromatic_rings, has_heteroatoms.
        If SMILES is invalid, all values are None.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "nrotbonds": None,
            "tpsa": None,
            "aromatic_rings": None,
            "has_heteroatoms": None,
        }

    return {
        "nrotbonds": Descriptors.NumRotatableBonds(mol),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "aromatic_rings": Descriptors.NumAromaticRings(mol),
        "has_heteroatoms": any(
            atom.GetAtomicNum() not in (1, 6) for atom in mol.GetAtoms()
        ),
    }


def initialize_batch_csv(
    input_path: Path,
    output_path: Path,
    limit: int | None = None,
    max_retries: int = 3,
) -> int:
    """Initialize batch tracking CSV from source data.

    Args:
        input_path: Path to source CSV (thermo_cbs_clean.csv)
        output_path: Path to output batch tracking CSV
        limit: Optional limit on number of molecules
        max_retries: Max retry attempts before Skip

    Returns:
        Number of molecules written
    """
    LOG.info(f"Reading source data from {input_path}")

    # Read source CSV
    df_source = pd.read_csv(input_path)

    # Validate required columns
    required_cols = ["smiles", "nheavy"]
    missing_cols = set(required_cols) - set(df_source.columns)

    if missing_cols:
        missing_sorted = sorted(missing_cols)
        expected_sorted = sorted(required_cols)
        found_sorted = sorted(df_source.columns)

        raise ValueError(
            f"Input CSV missing required columns: {missing_sorted}\n"
            f"Expected columns: {expected_sorted}\n"
            f"Found columns: {found_sorted}"
        )

    # Determine mol_id column
    if "Unnamed: 0" in df_source.columns:
        df_source["mol_id"] = df_source["Unnamed: 0"].astype(str)
    elif "mol_id" in df_source.columns:
        df_source["mol_id"] = df_source["mol_id"].astype(str)
    else:
        # Generate mol_ids
        df_source["mol_id"] = [f"mol_{i:06d}" for i in range(len(df_source))]

    # Apply limit if specified
    if limit is not None:
        df_source = df_source.head(limit)
        LOG.info(f"Limited to {limit} molecules")

    LOG.info(f"Processing {len(df_source)} molecules")

    # Compute additional descriptors
    LOG.info("Computing molecular descriptors...")
    df_descriptors = pd.DataFrame([
        compute_descriptors(smiles)
        for smiles in df_source["smiles"]
    ])

    # Create output DataFrame with all BatchRowCSV columns
    df_output = pd.DataFrame()

    # === Identification ===
    df_output["mol_id"] = df_source["mol_id"]
    df_output["smiles"] = df_source["smiles"]

    # === Molecular Descriptors ===
    df_output["nheavy"] = df_source["nheavy"]
    df_output["nrotbonds"] = df_descriptors["nrotbonds"]
    df_output["tpsa"] = df_descriptors["tpsa"]
    df_output["aromatic_rings"] = df_descriptors["aromatic_rings"]
    df_output["has_heteroatoms"] = df_descriptors["has_heteroatoms"]

    # === Reference Data ===
    if "H298_cbs" in df_source.columns:
        df_output["reference_hof"] = df_source["H298_cbs"]
    else:
        df_output["reference_hof"] = None

    # === Batch Status (initialized) ===
    df_output["status"] = "Pending"
    df_output["retry_count"] = 0
    df_output["max_retries"] = max_retries

    # === Batch Assignment (NULL initially) ===
    df_output["batch_id"] = None
    df_output["batch_order"] = None
    df_output["batch_failure_policy"] = None

    # === Timeout Configuration (NULL until assigned) ===
    df_output["assigned_crest_timeout"] = None
    df_output["assigned_mopac_timeout"] = None

    # === CREST Execution Results (NULL until processed) ===
    df_output["crest_status"] = None
    df_output["crest_conformers_generated"] = None
    df_output["crest_time"] = None
    df_output["crest_error"] = None

    # === MOPAC Execution Results (NULL until processed) ===
    df_output["num_conformers_selected"] = None
    df_output["most_stable_hof"] = None
    df_output["quality_grade"] = None

    # === Delta-E (NULL until processed) ===
    df_output["delta_e_12"] = None
    df_output["delta_e_13"] = None
    df_output["delta_e_15"] = None

    # === Final Status (NULL until processed) ===
    df_output["success"] = None
    df_output["error_message"] = None
    df_output["total_execution_time"] = None

    # === Actual Timeouts Used (NULL until processed) ===
    df_output["actual_crest_timeout_used"] = None
    df_output["actual_mopac_timeout_used"] = None

    # === Timestamps (NULL until processed) ===
    df_output["timestamp"] = None

    # === Error Tracking ===
    df_output["last_error_message"] = None

    # Write output CSV
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_output.to_csv(output_path, index=False)

    LOG.info(f"Wrote {len(df_output)} molecules to {output_path}")
    LOG.info(f"Columns: {list(df_output.columns)}")

    # Summary statistics
    LOG.info("Summary:")
    LOG.info(f"  Total molecules: {len(df_output)}")
    LOG.info(
        f"  nheavy range: {df_output['nheavy'].min()} - "
        f"{df_output['nheavy'].max()}"
    )
    LOG.info(
        f"  nrotbonds range: {df_output['nrotbonds'].min()} - "
        f"{df_output['nrotbonds'].max()}"
    )
    if df_output["reference_hof"].notna().any():
        LOG.info(
            f"  reference_hof range: {df_output['reference_hof'].min():.2f} - "
            f"{df_output['reference_hof'].max():.2f}"
        )

    return len(df_output)


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Initialize batch tracking CSV from source data."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        default=Path("data/thermo_cbs_clean.csv"),
        help="Path to source CSV file",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("data/batch_tracking.csv"),
        help="Path to output batch tracking CSV",
    )
    parser.add_argument(
        "--limit",
        "-l",
        type=int,
        default=None,
        help="Limit number of molecules (for testing)",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Max retry attempts before Skip (default: 3)",
    )

    args = parser.parse_args()

    if not args.input.exists():
        LOG.error(f"Input file not found: {args.input}")
        return 1

    try:
        count = initialize_batch_csv(
            input_path=args.input,
            output_path=args.output,
            limit=args.limit,
            max_retries=args.max_retries,
        )
        LOG.info(f"Successfully initialized {count} molecules")
        return 0
    except Exception as e:
        LOG.error(f"Failed to initialize CSV: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
