"""Dataset manager for source-of-truth CSV synchronization.

Manages synchronization between the canonical source CSV (thermo_cbs_chon.csv)
and the working CSV used for batch processing.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

LOG = logging.getLogger(__name__)


class DatasetManager:
    """Manage source-of-truth synchronization with thermo_cbs_chon.csv.

    Follows BatchCSVManager patterns for consistency in CSV handling.
    The source CSV contains reference molecular data, while the working CSV
    contains the same molecules plus tracking columns for batch processing.

    Attributes:
        source_csv: Path to the canonical source CSV file.
        working_csv: Path to the working CSV for batch processing.
    """

    # Columns to copy from source to working CSV
    SOURCE_COLUMNS = [
        "smiles",
        "nheavy",
        "multiplicity",
        "charge",
        "H298_cbs",
        "H298_b3",
    ]

    def __init__(self, source_csv: Path, working_csv: Path) -> None:
        """Initialize dataset manager.

        Args:
            source_csv: Path to canonical source CSV file.
            working_csv: Path to working CSV for batch processing.
        """
        self.source_csv = Path(source_csv)
        self.working_csv = Path(working_csv)

    def get_molecule_count(self) -> int:
        """Get total molecule count in source CSV.

        Returns:
            Number of molecules in the source CSV.

        Raises:
            FileNotFoundError: If source CSV does not exist.
        """
        if not self.source_csv.exists():
            msg = f"Source CSV not found: {self.source_csv}"
            LOG.error(msg)
            raise FileNotFoundError(msg)

        df = pd.read_csv(self.source_csv)
        count = len(df)
        LOG.info("Source CSV contains %d molecules", count)
        return count

    def sync_from_source(self, n_molecules: int) -> pd.DataFrame:
        """Read exactly first N molecules from source CSV.

        Args:
            n_molecules: Number of molecules to read from source.

        Returns:
            DataFrame with first N molecules, indexed with mol_id format.

        Raises:
            FileNotFoundError: If source CSV does not exist.
            ValueError: If n_molecules exceeds available molecules.
        """
        if not self.source_csv.exists():
            msg = f"Source CSV not found: {self.source_csv}"
            LOG.error(msg)
            raise FileNotFoundError(msg)

        df = pd.read_csv(self.source_csv)
        total = len(df)

        if n_molecules > total:
            msg = f"Requested {n_molecules} molecules but source has only {total}"
            LOG.error(msg)
            raise ValueError(msg)

        # Select first N molecules
        df = df.head(n_molecules).copy()

        # Generate mol_id with 5-digit zero-padded format
        df["mol_id"] = [f"mol_{i:05d}" for i in range(1, n_molecules + 1)]

        # Add status column
        df["status"] = "PENDING"

        # Reorder columns: mol_id first, then status, then source columns
        cols_order = ["mol_id", "status"]
        for col in self.SOURCE_COLUMNS:
            if col in df.columns:
                cols_order.append(col)

        # Add any remaining columns not in SOURCE_COLUMNS
        for col in df.columns:
            if col not in cols_order and col not in ["Unnamed: 0"]:
                cols_order.append(col)

        df = df[cols_order]

        LOG.info("Synced %d molecules from source", n_molecules)
        return df

    def create_working_csv(
        self, n_molecules: int, settings: dict[str, Any] | None = None
    ) -> None:
        """Create working CSV with first N source molecules.

        Creates a new working CSV file with molecules from the source,
        formatted for batch processing with mol_id and status columns.

        Args:
            n_molecules: Number of molecules to include in working CSV.
            settings: Optional settings dict (reserved for future use).

        Raises:
            FileNotFoundError: If source CSV does not exist.
            ValueError: If n_molecules exceeds available molecules.
        """
        _ = settings  # Reserved for future use
        df = self.sync_from_source(n_molecules)

        # Ensure parent directory exists
        self.working_csv.parent.mkdir(parents=True, exist_ok=True)

        # Save working CSV
        df.to_csv(self.working_csv, index=False)
        LOG.info(
            "Created working CSV with %d molecules: %s", n_molecules, self.working_csv
        )

    def validate_sync(self) -> bool:
        """Check if working CSV matches source (same SMILES, order).

        Validates that the working CSV SMILES column matches the corresponding
        entries in the source CSV, preserving order.

        Returns:
            True if working CSV matches source, False otherwise.
        """
        if not self.source_csv.exists():
            LOG.error("Source CSV not found: %s", self.source_csv)
            return False

        if not self.working_csv.exists():
            LOG.error("Working CSV not found: %s", self.working_csv)
            return False

        try:
            source_df = pd.read_csv(self.source_csv)
            working_df = pd.read_csv(self.working_csv)
        except Exception as e:
            LOG.error("Error reading CSV files: %s", e)
            return False

        n_working = len(working_df)

        if n_working > len(source_df):
            LOG.error(
                "Working CSV has %d molecules but source has only %d",
                n_working,
                len(source_df),
            )
            return False

        # Compare SMILES columns (first n_working entries)
        source_smiles = source_df["smiles"].head(n_working).tolist()
        working_smiles = working_df["smiles"].tolist()

        if source_smiles != working_smiles:
            LOG.error("SMILES mismatch between working and source CSV")
            return False

        LOG.info("Working CSV is in sync with source (%d molecules)", n_working)
        return True

    def refresh_database(self) -> None:
        """Resync working CSV from source.

        Recreates the working CSV using the same number of molecules
        as currently in the working CSV. If working CSV doesn't exist,
        uses the full source CSV.

        Raises:
            FileNotFoundError: If source CSV does not exist.
        """
        if not self.source_csv.exists():
            msg = f"Source CSV not found: {self.source_csv}"
            LOG.error(msg)
            raise FileNotFoundError(msg)

        # Determine how many molecules to sync
        if self.working_csv.exists():
            try:
                working_df = pd.read_csv(self.working_csv)
                n_molecules = len(working_df)
            except Exception:
                # If working CSV is corrupted, use source count
                n_molecules = self.get_molecule_count()
        else:
            n_molecules = self.get_molecule_count()

        # Recreate working CSV
        self.create_working_csv(n_molecules)
        LOG.info("Database refreshed with %d molecules", n_molecules)
