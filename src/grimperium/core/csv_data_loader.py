"""
CSV data loader with strict/permissive validation modes.

This module provides flexible CSV loading with two validation modes:
- strict: Raise on any validation error (development)
- permissive: Skip invalid rows, continue processing (production)

**AJUSTE v2.2:** _check_duplicates() integrated into validation flow
- Duplicates always fail, even in permissive mode (data corruption)
- Other validation errors can be skipped in permissive mode
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from grimperium.core.molecule import Molecule

import logging

import pandas as pd

logger = logging.getLogger(__name__)


class CSVDataLoaderError(Exception):
    """CSV data loader error."""

    pass


@dataclass
class ValidationError:
    """A single validation error."""

    row: int
    mol_id: str
    error: str


@dataclass
class ValidationReport:
    """Validation report with all errors."""

    errors: list[ValidationError] = field(default_factory=list)

    @property
    def total_errors(self) -> int:
        """Total number of validation errors."""
        return len(self.errors)

    def add_error(self, row: int, mol_id: str, error: str) -> None:
        """Add a validation error."""
        self.errors.append(ValidationError(row=row, mol_id=mol_id, error=error))


class CSVDataLoader:
    """
    Load molecules from CSV with flexible validation.

    Two modes:
    - strict=True: Raise CSVDataLoaderError on any invalid row
    - strict=False: Skip invalid rows, return valid ones only

    **AJUSTE v2.2:** Duplicate mol_ids always fail (even in permissive mode)

    Usage:
        >>> loader = CSVDataLoader(Path("data/thermo_pm7.csv"), strict=False)
        >>> df = loader.load_dataframe()
        >>> report = loader.get_validation_report()
        >>> if report.total_errors > 0:
        ...     print(f"Skipped {report.total_errors} rows")
    """

    # Required columns for a valid CSV
    REQUIRED_COLUMNS: frozenset[str] = frozenset(
        {
            "mol_id",
            "smiles",
            "nheavy",
            "status",
        }
    )

    # Optional columns (will be added if missing)
    OPTIONAL_COLUMNS: dict[str, Any] = {
        "charge": 0,
        "multiplicity": 1,
        "H298_cbs": None,
        "batch_id": "",
        "reruns": 0,
        "nrotbonds": 0,
        "timestamp_added": None,
        "timestamp_started": None,
        "timestamp_completed": None,
        "crest_status": None,
        "crest_conformers_generated": None,
        "crest_time": None,
        "mopac_status": None,
        "mopac_time": None,
        "delta_1": None,
        "delta_2": None,
        "delta_3": None,
        "most_stable_hof": None,
        "error_message": None,
    }

    def __init__(
        self,
        csv_path: Path,
        strict: bool = True,
    ):
        """
        Initialize CSV loader.

        Args:
            csv_path: Path to CSV file
            strict: If True, raise on any validation error.
                   If False, skip invalid rows and continue.
        """
        self.csv_path = Path(csv_path)
        self.strict = strict
        self.validation_report = ValidationReport()

        logger.info(f"CSVDataLoader initialized: {csv_path} (strict={strict})")

    def load_dataframe(self) -> pd.DataFrame:
        """
        Load and validate CSV data.

        Returns:
            Validated DataFrame with all columns.

        Raises:
            CSVDataLoaderError: If file not found, invalid format,
                               or (strict mode) any validation error.
        """
        # Check file exists
        if not self.csv_path.exists():
            raise CSVDataLoaderError(f"CSV file not found: {self.csv_path}")

        # Load CSV
        try:
            df = pd.read_csv(self.csv_path)
        except Exception as e:
            raise CSVDataLoaderError(f"Failed to parse CSV: {e}") from e

        logger.info(f"Loaded {len(df)} rows from {self.csv_path}")

        # Validate columns
        self._validate_columns(df)

        # Add missing optional columns
        df = self._add_missing_columns(df)

        # Validate data
        if self.strict:
            self._validate_data_strict(df)
        else:
            df = self._validate_data_permissive(df)

        return df

    def get_validation_report(self) -> ValidationReport:
        """Get validation report with all errors."""
        return self.validation_report

    def _validate_columns(self, df: pd.DataFrame) -> None:
        """Check required columns exist."""
        missing = self.REQUIRED_COLUMNS - set(df.columns)
        if missing:
            raise CSVDataLoaderError(
                f"Missing required columns: {sorted(missing)}\n"
                f"Required: {sorted(self.REQUIRED_COLUMNS)}\n"
                f"Found: {sorted(df.columns)}"
            )

    def _add_missing_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add missing optional columns with default values."""
        for col, default in self.OPTIONAL_COLUMNS.items():
            if col not in df.columns:
                df[col] = default
                logger.debug(f"Added missing column '{col}' with default: {default}")
        return df

    def _validate_data_strict(self, df: pd.DataFrame) -> None:
        """
        Strict validation: raise on any error.

        **AJUSTE v2.2:** Check duplicates first (critical error)
        """
        # AJUSTE #3: Integrar _check_duplicates first
        self._check_duplicates(df)

        # Then row-by-row validation
        for idx, row in df.iterrows():
            error = self._validate_row(row, int(idx))
            if error:
                raise CSVDataLoaderError(error)

    def _validate_data_permissive(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Permissive validation: skip invalid rows, continue.

        **AJUSTE v2.2:** Check duplicates first (always fail, even permissive)
        """
        # AJUSTE #3: Duplicates always fail (data corruption, not transient)
        try:
            self._check_duplicates(df)
        except CSVDataLoaderError as e:
            logger.error(f"Duplicate check failed (critical error): {e}")
            # In permissive mode, we still fail on duplicates
            # because it indicates data corruption, not validation error
            raise

        # Then row-by-row validation (can skip)
        valid_indices = []

        for idx, row in df.iterrows():
            error = self._validate_row(row, int(idx))

            if error:
                mol_id = row.get("mol_id", "unknown")
                logger.warning(f"Row {idx}: {error}")
                self.validation_report.add_error(int(idx), str(mol_id), error)
            else:
                valid_indices.append(idx)

        # Filter to valid rows
        df_valid = df.loc[valid_indices].reset_index(drop=True)
        skipped = len(df) - len(df_valid)

        logger.info(
            f"Permissive validation: {len(df_valid)} valid rows, "
            f"{skipped} rows skipped"
        )

        return df_valid

    def _validate_row(self, row: pd.Series, idx: int) -> str | None:  # noqa: ARG002
        """
        Validate a single row.

        Args:
            row: DataFrame row
            idx: Row index

        Returns:
            Error message if invalid, None if valid.
        """
        from grimperium.core.value_converter import MoleculeValueConverter
        from grimperium.crest_pm7.batch.enums import MoleculeStatus

        # mol_id is required
        mol_id, err = MoleculeValueConverter.to_string(
            row.get("mol_id"), "mol_id", allow_empty=False
        )
        if err:
            return f"mol_id is {err}: {row.get('mol_id')}"

        # smiles is required
        smiles, err = MoleculeValueConverter.to_string(
            row.get("smiles"), "smiles", allow_empty=False
        )
        if err:
            return f"smiles is {err}: {row.get('smiles')}"

        # nheavy is required and must be positive
        nheavy, err = MoleculeValueConverter.to_int(row.get("nheavy"), "nheavy")
        if err:
            return f"nheavy is {err}: {row.get('nheavy')}"
        if nheavy is not None and nheavy <= 0:
            return f"nheavy must be positive: {nheavy}"

        # status must be valid enum
        status, err = MoleculeValueConverter.to_enum(
            row.get("status"), MoleculeStatus, "status"
        )
        if err:
            return f"status is {err}: {row.get('status')}"

        # timestamp_added must be valid ISO 8601 if present
        timestamp_added = row.get("timestamp_added")
        if timestamp_added and not pd.isna(timestamp_added):
            _, err = MoleculeValueConverter.to_datetime(
                timestamp_added, "timestamp_added"
            )
            if err:
                return f"timestamp_added is {err}: {timestamp_added}"

        # All validations passed
        return None

    def _check_duplicates(self, df: pd.DataFrame) -> None:
        """
        Check for duplicate mol_ids.

        **AJUSTE v2.2:** Now integrated, always called in validation

        Duplicates always fail (even in permissive mode).
        Indicates data corruption, not transient error.

        Raises:
            CSVDataLoaderError: if duplicates found
        """
        duplicates = df[df.duplicated(subset=["mol_id"], keep=False)]

        if not duplicates.empty:
            dup_ids = duplicates["mol_id"].unique()
            dup_count = len(duplicates)

            raise CSVDataLoaderError(
                f"Duplicate mol_id found: {list(dup_ids)[:5]} "
                f"({dup_count} total duplicate rows)\n"
                f"This indicates data corruption. "
                f"Remove duplicates manually and retry.\n"
                f"Action: \n"
                f"1. Backup CSV: cp {self.csv_path} {self.csv_path}.bak\n"
                f"2. Remove duplicates: keep first occurrence only\n"
                f"3. Run grimperium validate --strict"
            )


class BatchDataManager:
    """
    Manage batch data loading and conversion to Molecule objects.

    Uses CSVDataLoader for validation and converts DataFrame rows
    to Molecule dataclass instances.

    Usage:
        >>> manager = BatchDataManager(Path("data/thermo_pm7.csv"), strict=False)
        >>> molecules = manager.load_batch()
        >>> report = manager.loader.get_validation_report()
    """

    def __init__(
        self,
        csv_path: Path,
        strict: bool = True,
    ):
        """
        Initialize batch data manager.

        Args:
            csv_path: Path to CSV file
            strict: Validation mode (strict or permissive)
        """
        self.csv_path = Path(csv_path)
        self.loader = CSVDataLoader(csv_path, strict=strict)
        self._df: pd.DataFrame | None = None

    def load_batch(self) -> list[Molecule]:
        """
        Load CSV and convert to Molecule objects.

        Returns:
            List of Molecule instances.

        Raises:
            CSVDataLoaderError: If loading fails.
        """
        from grimperium.core.molecule import Molecule

        # Load and validate
        self._df = self.loader.load_dataframe()

        # Convert to Molecule objects
        molecules = []
        for idx, row in self._df.iterrows():
            try:
                mol = Molecule.from_csv_dict(row.to_dict())
                molecules.append(mol)
            except (ValueError, KeyError) as e:
                logger.error(f"Row {idx}: Failed to create Molecule: {e}")
                if self.loader.strict:
                    raise CSVDataLoaderError(f"Row {idx}: {e}") from e
                # In permissive mode, skip this row
                continue

        logger.info(f"Created {len(molecules)} Molecule objects")
        return molecules

    def count_by_status(self) -> dict[str, int]:
        """
        Count molecules by status.

        Returns:
            Dictionary with status counts.
        """
        if self._df is None:
            return {}

        counts = self._df["status"].value_counts().to_dict()

        # Normalize keys to lowercase
        return {str(k).lower(): int(v) for k, v in counts.items()}

    def get_dataframe(self) -> pd.DataFrame | None:
        """Get the loaded DataFrame (or None if not loaded)."""
        return self._df
