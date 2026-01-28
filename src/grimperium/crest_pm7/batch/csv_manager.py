"""CSV manager for batch processing status tracking.

This module provides BatchCSVManager for:
- Loading and saving CSV tracking files
- Selecting molecules for batches
- Managing status transitions
- Mapping PM7Result to CSV columns

FUTURE PARALLELIZATION: Add threading.Lock() to protect DataFrame operations.
"""

import logging
import math
from pathlib import Path
from typing import Any

import pandas as pd

from grimperium.crest_pm7.batch.enums import (
    BatchFailurePolicy,
    BatchSortingStrategy,
    MoleculeStatus,
)
from grimperium.crest_pm7.batch.models import Batch, BatchMolecule
from grimperium.crest_pm7.progress import (
    CREST_STATUS_NOT_ATTEMPTED,
    CREST_STATUS_PREOPT,
    CREST_STATUS_SEARCH,
    MOPAC_STATUS_NOT_ATTEMPTED,
    MOPAC_STATUS_RUNNING,
)

LOG = logging.getLogger("grimperium.crest_pm7.batch.csv_manager")


class BatchCSVManager:
    """Manages CSV tracking for batch processing.

    Handles:
    - CSV file I/O with status tracking
    - Batch selection with configurable strategies
    - Status transitions (PENDING -> SELECTED -> RUNNING -> OK/RERUN/SKIP)
    - PM7Result to CSV column mapping

    Attributes:
        csv_path: Path to the tracking CSV file
        df: Pandas DataFrame with molecule data (loaded lazily)
    """

    # Result columns that are cleared when resetting a batch
    RESULT_COLUMNS = [
        "crest_status",
        "crest_conformers_generated",
        "crest_time",
        "crest_error",
        "mopac_status",  # NEW: Overall MOPAC status
        "num_conformers_selected",
        "mopac_time",  # NEW: Aggregated MOPAC execution time
        "H298_pm7",  # Renamed from most_stable_hof
        "abs_diff",  # NEW
        "abs_diff_%",  # NEW
        "quality_grade",
        "success",
        "error_message",
        "total_execution_time",
        "assigned_crest_timeout",  # Renamed from actual_crest_timeout_used
        "assigned_mopac_timeout",  # Renamed from actual_mopac_timeout_used
        "delta_1",  # Renamed from delta_e_12
        "delta_2",  # Renamed from delta_e_13
        "delta_3",  # Renamed from delta_e_15
        "conformer_selected",  # NEW: Which conformer was selected (1-3)
        "reruns",  # NEW
        "timestamp",
    ]

    # Identity column
    IDENTITY_COLUMNS = ["mol_id"]

    # Molecular properties columns
    MOLECULAR_PROPERTIES_COLUMNS = [
        "smiles",
        "nheavy",
        "nrotbonds",
        "tpsa",
        "aromatic_rings",
        "has_heteroatoms",
        "reference_hof",
    ]

    # Batch info columns
    BATCH_INFO_COLUMNS = [
        "status",
        "batch_id",
        "batch_order",
        "batch_failure_policy",
    ]

    # CREST configuration columns
    CREST_CONFIG_COLUMNS = [
        "v3",  # Renamed from crest_v3
        "qm",  # Renamed from crest_quick
        "nci",  # Renamed from crest_nci
        "c_method",  # Renamed from crest_gfnff
        "energy_window",  # Renamed from crest_ewin
        "rmsd_threshold",  # Renamed from crest_rthr
        "crest_optlev",  # Keep old name for now
        "threads",  # Renamed from crest_threads
        "xtb",  # Renamed from crest_xtb_preopt
    ]

    # MOPAC configuration columns
    MOPAC_CONFIG_COLUMNS = [
        "precise_scf",  # Renamed from mopac_precise
        "scf_threshold",  # Renamed from mopac_scfcrt
        # Note: mopac_itry, mopac_pulay, mopac_prtall, mopac_archive not in Phase A CSV
    ]

    # Retry tracking columns
    RETRY_TRACKING_COLUMNS = [
        "retry_count",
        "last_error_message",
        "max_retries",
    ]

    # Phase B reserved columns (for future ML features)
    PHASE_B_RESERVED_COLUMNS = [
        "reserved_42",
        "reserved_43",
        "reserved_44",
        "reserved_45",
        "reserved_46",
        "reserved_47",
        "reserved_48",
        "reserved_49",
    ]

    def __init__(self, csv_path: Path | None) -> None:
        """Initialize CSV manager.

        Args:
            csv_path: Path to CSV tracking file (can be None for schema-only use)
        """
        self.csv_path = Path(csv_path) if csv_path is not None else None
        self.df: pd.DataFrame | None = None
        # FUTURE PARALLELIZATION: Add threading.Lock() here
        # self._lock = threading.Lock()

    def get_schema(self) -> list[str]:
        """Get the full CSV schema with all column names.

        Returns:
            List of 49 column names in order (Phase A schema):
            - Identity (1): mol_id
            - Molecular properties (7): smiles, nheavy, nrotbonds, tpsa, etc.
            - Batch info (4): status, batch_id, batch_order, batch_failure_policy
            - CREST configuration (9): v3, qm, nci, c_method, energy_window, etc.
            - MOPAC configuration (2): precise_scf, scf_threshold
            - Results (18): crest_status, H298_pm7, deltas, mopac_time, etc.
            - Retry tracking (2): retry_count, last_error_message
            - Phase B reserved (8): reserved_42 through reserved_49
        """
        return (
            self.IDENTITY_COLUMNS
            + self.MOLECULAR_PROPERTIES_COLUMNS
            + self.BATCH_INFO_COLUMNS
            + self.CREST_CONFIG_COLUMNS
            + self.MOPAC_CONFIG_COLUMNS
            + self.RESULT_COLUMNS
            + self.RETRY_TRACKING_COLUMNS
            + self.PHASE_B_RESERVED_COLUMNS
        )

    def load_csv(self) -> pd.DataFrame:
        """Load CSV file into DataFrame.

        Returns:
            Loaded DataFrame

        Raises:
            FileNotFoundError: If CSV file doesn't exist
            RuntimeError: If csv_path is None
        """
        if self.csv_path is None:
            raise RuntimeError("csv_path is None - cannot load CSV")
        if not self.csv_path.exists():
            raise FileNotFoundError(f"CSV file not found: {self.csv_path}")

        self.df = pd.read_csv(
            self.csv_path,
            dtype={
                "mol_id": str,
                "smiles": str,
                "status": str,
                "batch_id": str,
                "batch_failure_policy": str,
                "crest_status": str,
                "quality_grade": str,
                "crest_error": str,
                "error_message": str,
                "last_error_message": str,
                # Force string type for timestamp to avoid float64 issues
                "timestamp": str,
            },
        )

        # Normalize column names (trim whitespace) to avoid schema mismatches
        original_columns = list(self.df.columns)
        normalized_columns = [
            col.strip() if isinstance(col, str) else col for col in original_columns
        ]
        if normalized_columns != original_columns:
            if len(set(normalized_columns)) != len(normalized_columns):
                duplicates = [
                    col
                    for col in set(normalized_columns)
                    if normalized_columns.count(col) > 1
                ]
                raise ValueError(
                    "CSV has duplicate columns after normalization: "
                    f"{', '.join(sorted(duplicates))}"
                )
            self.df.columns = normalized_columns
            LOG.warning(
                "Normalized CSV column names (trimmed whitespace) to match schema."
            )

        # Ensure boolean columns use object dtype to avoid FutureWarning
        # when assigning True/False to columns that may contain NaN
        bool_columns = ["success"]
        for col in bool_columns:
            if col in self.df.columns:
                # Convert to object dtype to allow mixed bool/NaN values
                self.df[col] = self.df[col].astype(object)

        # Validate required columns
        required_cols = {"mol_id", "smiles", "nheavy", "status"}
        missing = required_cols - set(self.df.columns)
        if missing:
            raise ValueError(f"CSV missing required columns: {missing}")

        # Normalize status column (case-insensitive mapping to enum values)
        status_map = {s.value.lower(): s.value for s in MoleculeStatus}
        self.df["status"] = self.df["status"].str.strip().str.lower().map(status_map)
        invalid_status = self.df["status"].isna()
        if invalid_status.any():
            first_invalid_idx = self.df[invalid_status].index[0]
            LOG.warning(
                f"Found {invalid_status.sum()} rows with invalid status values "
                f"(first at row {first_invalid_idx}). Setting to PENDING."
            )
            self.df.loc[invalid_status, "status"] = MoleculeStatus.PENDING.value

        # Validate unique mol_ids
        if self.df["mol_id"].duplicated().any():
            duplicates = self.df[self.df["mol_id"].duplicated()]["mol_id"].tolist()
            display_list = duplicates[:5]
            msg = f"Duplicate mol_ids found: {display_list}"
            if len(duplicates) > 5:
                msg += f" ... and {len(duplicates) - 5} more"
            raise ValueError(msg)

        LOG.info(f"Loaded {len(self.df)} molecules from {self.csv_path}")
        return self.df

    def save_csv(self) -> None:
        """Save DataFrame to CSV file."""
        if self.csv_path is None:
            raise RuntimeError("csv_path is None - cannot save CSV")
        if self.df is None:
            raise RuntimeError("No DataFrame loaded - call load_csv() first")

        self.df.to_csv(self.csv_path, index=False)
        LOG.debug(f"Saved {len(self.df)} molecules to {self.csv_path}")

    def _ensure_loaded(self) -> pd.DataFrame:
        """Ensure DataFrame is loaded."""
        if self.df is None:
            self.load_csv()
        return self.df

    def _safe_int(self, val: Any, default: int = 0) -> int:
        """Safely convert value to int, handling NaN and invalid types.

        WARNING: Truncates decimal values (e.g., 3.9 → 3).
        Use with caution on CSV fields that may contain floats.
        A LOG.warning is emitted when truncation occurs.
        Uses math.isclose() to avoid spurious warnings from floating-point precision.

        Args:
            val: Value to convert (Any type)
            default: Default to return if conversion fails

        Returns:
            int: Converted value (truncated if float), or default

        Behavior:
            1. If NaN/None → return default
            2. Try int(val) → return if succeeds
            3. Try int(float(val)) → check for fractional part:
               - If has fraction (e.g., 3.9) → LOG.warning + return truncated int
               - If no fraction (e.g., 3.0) → return int silently
            4. If all fail → LOG.warning + return default
        """
        if pd.isna(val):
            return default

        try:
            # Direct conversion works for integers
            return int(val)
        except (ValueError, TypeError):
            # Try float→int conversion (handles "3.0", 3.5, etc.)
            try:
                if isinstance(val, str):
                    val = val.strip()  # Remove whitespace

                float_val = float(val)
                int_val = int(float_val)

                # Warn if truncating non-integer float
                # Use math.isclose to avoid spurious float precision warnings
                if not math.isclose(float_val, int_val, rel_tol=1e-9, abs_tol=1e-9):
                    LOG.warning(
                        f"Truncating float {float_val} to int {int_val}. "
                        f"Original value: {val}. This may cause data loss. "
                        f"Consider validating or rounding before CSV import."
                    )

                return int_val

            except (ValueError, TypeError):
                # Last resort: log and return default
                LOG.warning(
                    f"Cannot convert '{val}' to int, using default {default}. "
                    f"Type: {type(val).__name__}"
                )
                return default

    def _get_row_index(self, mol_id: str) -> int:
        """Get DataFrame index for mol_id.

        Args:
            mol_id: Molecule identifier

        Returns:
            Row index

        Raises:
            KeyError: If mol_id not found
        """
        df = self._ensure_loaded()
        mask = df["mol_id"] == mol_id
        if not mask.any():
            raise KeyError(f"mol_id not found: {mol_id}")
        return int(df[mask].index[0])

    def get_reference_hof(self, mol_id: str) -> float | None:
        """Get reference HOF (H298_cbs) value for a molecule from CSV.

        Args:
            mol_id: Molecule identifier

        Returns:
            Reference HOF value in kcal/mol, or None if not available/invalid
        """
        try:
            df = self._ensure_loaded()
            idx = self._get_row_index(mol_id)
            if "reference_hof" in df.columns:
                val = df.at[idx, "reference_hof"]
            elif "H298_cbs" in df.columns:
                val = df.at[idx, "H298_cbs"]
            else:
                LOG.debug(f"[{mol_id}] CSV missing reference_hof/H298_cbs columns")
                return None
            if pd.isna(val):
                return None
            return float(val)
        except (KeyError, ValueError, TypeError) as e:
            LOG.debug(f"[{mol_id}] Could not get reference_hof: {e}")
            return None

    def generate_batch_id(self) -> str:
        """Generate unique batch ID with sequential numbering.

        Returns:
            Batch ID in format 'batch_NNNN' where NNNN is 4-digit counter
            (e.g., 'batch_0001', 'batch_0002', ...)
        """
        df = self._ensure_loaded()

        # Extract existing batch numbers from batch_id column
        existing_batches = df["batch_id"].dropna().unique()
        batch_numbers = []

        for batch_id in existing_batches:
            # Extract number from format 'batch_NNNN'
            if isinstance(batch_id, str) and batch_id.startswith("batch_"):
                try:
                    num = int(batch_id.replace("batch_", ""))
                    batch_numbers.append(num)
                except ValueError:
                    # Skip malformed batch IDs
                    pass

        # Get next number
        next_number = max(batch_numbers, default=0) + 1
        return f"batch_{next_number:04d}"

    def select_batch(
        self,
        batch_id: str,
        batch_size: int,
        crest_timeout_minutes: int,
        mopac_timeout_minutes: int,
        strategy: BatchSortingStrategy = BatchSortingStrategy.RERUN_FIRST_THEN_EASY,
        failure_policy: BatchFailurePolicy = BatchFailurePolicy.PARTIAL_OK,
    ) -> Batch:
        """Select molecules for a new batch.

        Selects up to batch_size molecules from PENDING/RERUN status,
        applies sorting strategy, and marks them as SELECTED.

        Args:
            batch_id: Unique batch identifier
            batch_size: Maximum molecules to select
            crest_timeout_minutes: CREST timeout for all molecules
            mopac_timeout_minutes: MOPAC timeout for all molecules
            strategy: How to select and order molecules
            failure_policy: How to handle failures

        Returns:
            Batch with selected molecules

        Raises:
            ValueError: If parameters are invalid
            RuntimeError: If molecules are currently RUNNING
        """
        # Validate parameters
        if batch_size < 1:
            raise ValueError("batch_size must be >= 1")
        if crest_timeout_minutes <= 0:
            raise ValueError("crest_timeout_minutes must be > 0")
        if mopac_timeout_minutes <= 0:
            raise ValueError("mopac_timeout_minutes must be > 0")
        if not batch_id:
            raise ValueError("batch_id cannot be empty")

        df = self._ensure_loaded()

        # Check for running molecules (prevent concurrent batches)
        running_mask = df["status"] == MoleculeStatus.RUNNING.value
        if running_mask.any():
            raise RuntimeError(
                f"Cannot create batch: {running_mask.sum()} molecules still RUNNING"
            )

        # Find available molecules (PENDING or RERUN)
        available_mask = df["status"].isin(
            [
                MoleculeStatus.PENDING.value,
                MoleculeStatus.RERUN.value,
            ]
        )
        available_df = df[available_mask].copy()

        if available_df.empty:
            LOG.warning("No pending/rerun molecules available for batch")
            return Batch(
                batch_id=batch_id,
                molecules=[],
                crest_timeout_minutes=crest_timeout_minutes,
                mopac_timeout_minutes=mopac_timeout_minutes,
                sorting_strategy=strategy,
                failure_policy=failure_policy,
            )

        # Apply sorting strategy
        available_df = self._apply_sorting_strategy(available_df, strategy)

        # Limit to batch_size
        actual_size = min(batch_size, len(available_df))
        if actual_size < batch_size:
            LOG.info(
                f"Requested {batch_size} molecules, only {len(available_df)} available"
            )

        selected_df = available_df.head(actual_size)

        # Create BatchMolecule list
        molecules: list[BatchMolecule] = []
        for batch_order, (idx, row) in enumerate(selected_df.iterrows(), start=1):
            mol = BatchMolecule(
                mol_id=row["mol_id"],
                smiles=row["smiles"],
                batch_order=batch_order,
                nheavy=self._safe_int(row["nheavy"], 0),
                nrotbonds=self._safe_int(row.get("nrotbonds", 0), 0),
                retry_count=self._safe_int(row.get("retry_count", 0), 0),
            )
            molecules.append(mol)

            # Update DataFrame: mark as SELECTED
            df.at[idx, "status"] = MoleculeStatus.SELECTED.value
            df.at[idx, "batch_id"] = batch_id
            df.at[idx, "batch_order"] = batch_order
            df.at[idx, "batch_failure_policy"] = failure_policy.value
            df.at[idx, "assigned_crest_timeout"] = crest_timeout_minutes
            df.at[idx, "assigned_mopac_timeout"] = mopac_timeout_minutes

        self.save_csv()
        LOG.info(
            f"Selected {len(molecules)} molecules for batch {batch_id} "
            f"(strategy={strategy.value}, policy={failure_policy.value})"
        )

        return Batch(
            batch_id=batch_id,
            molecules=molecules,
            crest_timeout_minutes=crest_timeout_minutes,
            mopac_timeout_minutes=mopac_timeout_minutes,
            sorting_strategy=strategy,
            failure_policy=failure_policy,
        )

    def _apply_sorting_strategy(
        self, df: pd.DataFrame, strategy: BatchSortingStrategy
    ) -> pd.DataFrame:
        """Apply sorting strategy to available molecules.

        Respects mol_id order within each priority group.

        Args:
            df: DataFrame of available molecules
            strategy: Sorting strategy to apply

        Returns:
            Sorted DataFrame
        """
        if strategy == BatchSortingStrategy.RERUN_FIRST_THEN_EASY:
            # RERUN first (sorted by mol_id), then PENDING (sorted by mol_id)
            df = df.copy()
            df["_sort_priority"] = df["status"].apply(
                lambda x: 0 if x == MoleculeStatus.RERUN.value else 1
            )
            return df.sort_values(["_sort_priority", "mol_id"]).drop(
                columns=["_sort_priority"]
            )

        elif strategy == BatchSortingStrategy.RANDOM:
            return df.sample(frac=1)

        elif strategy == BatchSortingStrategy.BY_NHEAVY:
            return df.sort_values("nheavy")

        elif strategy == BatchSortingStrategy.BY_NROTBONDS:
            if "nrotbonds" in df.columns:
                return df.sort_values("nrotbonds")
            else:
                LOG.warning(
                    "Column 'nrotbonds' not found, skipping BY_NROTBONDS sorting"
                )
                return df

        else:
            LOG.warning(f"Unknown strategy {strategy}, using default order")
            return df

    def get_status(self, mol_id: str) -> str:
        """Get current status for a molecule.

        Args:
            mol_id: Molecule identifier

        Returns:
            MoleculeStatus value (string)
        """
        self._ensure_loaded()
        assert self.df is not None  # Guaranteed by _ensure_loaded()
        idx = self._get_row_index(mol_id)
        return str(self.df.at[idx, "status"])

    def mark_running(self, mol_id: str) -> None:
        """Mark molecule as currently running.

        Transition: SELECTED -> RUNNING

        Also resets crest_status and mopac_status to NOT_ATTEMPTED
        to ensure progress tracking starts cleanly for the run.

        Args:
            mol_id: Molecule identifier
        """
        # FUTURE PARALLELIZATION: Wrap with self._lock
        idx = self._get_row_index(mol_id)
        df = self._ensure_loaded()

        current_status = df.at[idx, "status"]
        if current_status != MoleculeStatus.SELECTED.value:
            LOG.warning(
                f"mark_running({mol_id}): expected SELECTED, got {current_status}"
            )

        df.at[idx, "status"] = MoleculeStatus.RUNNING.value
        if "crest_status" in df.columns:
            df.at[idx, "crest_status"] = CREST_STATUS_NOT_ATTEMPTED
        if "mopac_status" in df.columns:
            df.at[idx, "mopac_status"] = MOPAC_STATUS_NOT_ATTEMPTED
        self.save_csv()
        LOG.debug(f"Marked {mol_id} as RUNNING")

    def _update_progress_column(
        self,
        mol_id: str,
        column: str,
        value: str,
        *,
        allowed_current: set[str] | None = None,
    ) -> bool:
        """Update a progress column if the current value allows it.

        Args:
            mol_id: Molecule identifier
            column: CSV column name
            value: Value to write
            allowed_current: Optional set of allowed current values

        Returns:
            True if the column was updated, False otherwise
        """
        df = self._ensure_loaded()
        idx = self._get_row_index(mol_id)

        if column not in df.columns:
            return False

        current_raw = df.at[idx, column]
        current = "" if pd.isna(current_raw) else str(current_raw).strip()

        if (
            allowed_current is not None
            and current not in allowed_current
            and current != ""
        ):
            return False

        if current == value:
            return False

        df.at[idx, column] = value
        return True

    def mark_crest_preopt(self, mol_id: str) -> None:
        """Mark molecule as entering xTB pre-optimization stage.

        Args:
            mol_id: Molecule identifier
        """
        updated = self._update_progress_column(
            mol_id,
            "crest_status",
            CREST_STATUS_PREOPT,
            allowed_current={CREST_STATUS_NOT_ATTEMPTED, CREST_STATUS_PREOPT},
        )
        if updated:
            self.save_csv()

    def mark_crest_search(self, mol_id: str) -> None:
        """Mark molecule as entering CREST conformer search.

        Args:
            mol_id: Molecule identifier
        """
        updated = self._update_progress_column(
            mol_id,
            "crest_status",
            CREST_STATUS_SEARCH,
            allowed_current={
                CREST_STATUS_NOT_ATTEMPTED,
                CREST_STATUS_PREOPT,
                CREST_STATUS_SEARCH,
            },
        )
        if updated:
            self.save_csv()

    def mark_mopac_running(self, mol_id: str) -> None:
        """Mark molecule as entering MOPAC calculation stage.

        Args:
            mol_id: Molecule identifier
        """
        updated = self._update_progress_column(
            mol_id,
            "mopac_status",
            MOPAC_STATUS_RUNNING,
            allowed_current={MOPAC_STATUS_NOT_ATTEMPTED, MOPAC_STATUS_RUNNING},
        )
        if updated:
            self.save_csv()

    def mark_success(
        self,
        mol_id: str,
        result_update: dict[str, Any],
    ) -> None:
        """Mark molecule as successfully processed.

        Transition: RUNNING -> OK

        Args:
            mol_id: Molecule identifier
            result_update: Dict with CSV column updates from _pm7result_to_csv_update
        """
        # FUTURE PARALLELIZATION: Wrap with self._lock
        idx = self._get_row_index(mol_id)
        df = self._ensure_loaded()

        current_status = df.at[idx, "status"]
        if current_status != MoleculeStatus.RUNNING.value:
            LOG.warning(
                f"mark_success({mol_id}): expected RUNNING, got {current_status}"
            )

        # Update status
        df.at[idx, "status"] = MoleculeStatus.OK.value

        # Apply result updates
        for col, val in result_update.items():
            if col in df.columns:
                df.at[idx, col] = val

        self.save_csv()
        LOG.debug(f"Marked {mol_id} as OK")

    def mark_rerun(
        self,
        mol_id: str,
        error_message: str,
        result_update: dict[str, Any] | None = None,
    ) -> None:
        """Mark molecule for retry.

        Transition: RUNNING -> RERUN (if retry_count < max_retries)
        Transition: RUNNING -> SKIP (if retry_count >= max_retries)

        Args:
            mol_id: Molecule identifier
            error_message: Error that caused the failure
            result_update: Optional partial results to save
        """
        # FUTURE PARALLELIZATION: Wrap with self._lock
        idx = self._get_row_index(mol_id)
        df = self._ensure_loaded()

        current_status = df.at[idx, "status"]
        if current_status != MoleculeStatus.RUNNING.value:
            LOG.warning(f"mark_rerun({mol_id}): expected RUNNING, got {current_status}")

        # Increment retry count
        retry_count = self._safe_int(df.at[idx, "retry_count"], default=0) + 1

        # Handle max_retries explicitly (0 is valid)
        max_retries_val = df.at[idx, "max_retries"]
        if pd.isna(max_retries_val):
            max_retries = 3
        else:
            max_retries = int(max_retries_val)
            if max_retries < 0:
                raise ValueError(f"max_retries must be >= 0, got {max_retries}")

        df.at[idx, "retry_count"] = retry_count
        df.at[idx, "last_error_message"] = error_message

        # Determine next status
        if retry_count >= max_retries:
            df.at[idx, "status"] = MoleculeStatus.SKIP.value
            LOG.warning(f"Marked {mol_id} as SKIP (retry_count={retry_count} >= max)")
        else:
            df.at[idx, "status"] = MoleculeStatus.RERUN.value
            LOG.warning(
                f"Marked {mol_id} as RERUN "
                f"({retry_count}/{max_retries}): {error_message}"
            )

        # Apply partial results if provided
        if result_update:
            for col, val in result_update.items():
                if col in df.columns:
                    df.at[idx, col] = val

        self.save_csv()

    def mark_skip(
        self,
        mol_id: str,
        error_message: str,
    ) -> None:
        """Mark molecule as permanently skipped.

        Transition: RUNNING -> SKIP

        Args:
            mol_id: Molecule identifier
            error_message: Error that caused the skip
        """
        # FUTURE PARALLELIZATION: Wrap with self._lock
        idx = self._get_row_index(mol_id)
        df = self._ensure_loaded()

        df.at[idx, "status"] = MoleculeStatus.SKIP.value
        df.at[idx, "last_error_message"] = error_message
        df.at[idx, "error_message"] = error_message

        self.save_csv()
        LOG.warning(f"Marked {mol_id} as SKIP: {error_message}")

    def _update_extra_fields(
        self,
        mol_id: str,
        field_updates: dict[str, Any],
    ) -> None:
        """Update extra CSV fields for a molecule (internal helper).

        This method is used by csv_enhancements to add calculated deltas
        and batch settings after mark_success() has been called.

        Designed to be called AFTER status transition (e.g., after mark_success).
        Does NOT change molecule status, only updates data fields.

        Args:
            mol_id: Molecule identifier
            field_updates: Dict of {column_name: value} to update

        Raises:
            KeyError: If mol_id not found in CSV

        Example:
            >>> csv_manager.mark_success(mol_id, result_update)
            >>> csv_manager._update_extra_fields(mol_id, {
            ...     'delta_1': 0.0,
            ...     'delta_2': 0.45,
            ...     'v3': False,
            ...     'c_method': 'gfn2-xtb',
            ... })
        """
        # FUTURE PARALLELIZATION: Wrap with self._lock
        idx = self._get_row_index(mol_id)
        df = self._ensure_loaded()

        updated_count = 0
        skipped_cols = []

        for col, val in field_updates.items():
            if col in df.columns:
                df.at[idx, col] = val
                updated_count += 1
            else:
                skipped_cols.append(col)
                LOG.warning(f"Column '{col}' not in CSV schema, skipping")

        if updated_count == 0:
            LOG.debug(f"No fields updated for {mol_id}, skipping save")
            return

        self.save_csv()
        if skipped_cols:
            LOG.debug(
                f"Updated {updated_count} extra fields for {mol_id} "
                f"(skipped {len(skipped_cols)} unknown columns: {', '.join(skipped_cols)})"
            )
        else:
            LOG.debug(f"Updated {updated_count} extra fields for {mol_id}")

    def reset_batch(self, batch_id: str) -> int:
        """Reset all molecules in a batch for re-processing.

        Only called when:
        1. failure_policy == ALL_OR_NOTHING
        2. Batch completed with at least one failure

        Resets status to PENDING but keeps retry_count for audit.

        Args:
            batch_id: Batch to reset

        Returns:
            Number of molecules reset
        """
        df = self._ensure_loaded()

        mask = df["batch_id"] == batch_id
        batch_rows = df[mask]

        if batch_rows.empty:
            LOG.warning(f"No molecules found for batch {batch_id}")
            return 0

        # Verify policy
        policy = batch_rows["batch_failure_policy"].iloc[0]
        if policy != BatchFailurePolicy.ALL_OR_NOTHING.value:
            LOG.warning(f"reset_batch called but policy is {policy}, skipping")
            return 0

        reset_count = 0
        for idx in batch_rows.index:
            current_status = df.at[idx, "status"]

            # Reset to PENDING
            df.at[idx, "status"] = MoleculeStatus.PENDING.value
            df.at[idx, "batch_id"] = None
            df.at[idx, "batch_order"] = None

            # Clear execution results
            self._clear_execution_results(idx)

            # Keep retry_count and last_error_message for audit
            reset_count += 1

            LOG.debug(f"Reset {df.at[idx, 'mol_id']}: {current_status} -> PENDING")

        self.save_csv()
        LOG.info(f"Reset {reset_count} molecules from batch {batch_id}")
        return reset_count

    def _clear_execution_results(self, idx: int) -> None:
        """Clear execution result columns for a row.

        Args:
            idx: DataFrame row index
        """
        df = self._ensure_loaded()

        for col in self.RESULT_COLUMNS:
            if col in df.columns:
                df.at[idx, col] = None

    def pm7result_to_csv_update(
        self,
        mol_id: str,  # noqa: ARG002 - kept for API consistency
        result: Any,  # PM7Result - use Any to avoid circular import
        batch_id: str,
        batch_order: int,
        crest_timeout_used: float,
        mopac_timeout_used: float,
    ) -> dict[str, Any]:
        """Convert PM7Result to CSV update dict.

        This method extracts relevant fields from PM7Result for CSV update.
        Status management (OK/Rerun/Skip) is handled separately.

        Args:
            mol_id: Molecule identifier
            result: PM7Result from processing
            batch_id: Current batch ID
            batch_order: Position in batch (1-indexed)
            crest_timeout_used: Actual CREST timeout used (minutes)
            mopac_timeout_used: Actual MOPAC timeout used (minutes)

        Returns:
            Dict with CSV column updates
        """
        # Calculate metrics if H298_cbs and H298_pm7 are available
        h298_pm7 = (
            round(result.most_stable_hof, 2)
            if getattr(result, "most_stable_hof", None) is not None
            else None
        )
        h298_cbs = result.H298_cbs if hasattr(result, "H298_cbs") else None

        abs_diff = None
        abs_diff_pct = None
        if h298_pm7 is not None and h298_cbs is not None:
            abs_diff = round(abs(h298_cbs - h298_pm7), 4)
            if h298_cbs != 0:
                abs_diff_pct = round((abs_diff / abs(h298_cbs)) * 100, 2)

        # Calculate mopac_time (aggregate from all conformers)
        mopac_time = None
        if hasattr(result, "conformers") and result.conformers:
            times = [
                c.mopac_execution_time
                for c in result.conformers
                if hasattr(c, "mopac_execution_time")
                and c.mopac_execution_time is not None
            ]
            if times:
                mopac_time = round(sum(times), 1)

        # Calculate mopac_status (overall status based on conformers)
        mopac_status = None
        if hasattr(result, "conformers") and result.conformers:
            successful = [
                c
                for c in result.conformers
                if hasattr(c, "is_successful") and c.is_successful
            ]
            if successful:
                mopac_status = "OK"
            elif result.conformers:
                mopac_status = "FAILED"

        # NOTE: delta_1/2/3 and conformer_selected are computed later using
        # H298_cbs and conformer HOF values in CSV enhancements.
        conformer_selected = None

        return {
            # RDKit Descriptors (from PM7Result)
            "nrotbonds": result.nrotbonds if hasattr(result, "nrotbonds") else None,
            "tpsa": (
                round(result.tpsa, 2)
                if hasattr(result, "tpsa") and result.tpsa is not None
                else None
            ),
            "aromatic_rings": (
                result.aromatic_rings if hasattr(result, "aromatic_rings") else None
            ),
            # CREST Execution
            "crest_status": (
                result.crest_status.value if result.crest_status is not None else None
            ),
            "crest_conformers_generated": result.crest_conformers_generated,
            "crest_time": (
                round(result.crest_time, 1) if result.crest_time is not None else None
            ),
            "crest_error": result.crest_error,
            # MOPAC Execution
            "mopac_status": mopac_status,  # NEW: Computed from conformers
            "num_conformers_selected": result.num_conformers_selected,
            "mopac_time": mopac_time,  # NEW: Aggregated from all conformers
            "H298_pm7": h298_pm7,  # Renamed from most_stable_hof
            "abs_diff": abs_diff,  # NEW: |H298_cbs - H298_pm7|
            "abs_diff_%": abs_diff_pct,  # NEW: Percentage difference
            "quality_grade": (
                result.quality_grade.value if result.quality_grade is not None else None
            ),
            "success": result.success,
            "error_message": result.error_message,
            "total_execution_time": (
                round(result.total_execution_time, 1)
                if result.total_execution_time is not None
                else None
            ),
            "assigned_crest_timeout": round(crest_timeout_used, 1),  # Renamed
            "assigned_mopac_timeout": round(mopac_timeout_used, 1),  # Renamed
            # Delta-E (renamed for clarity)
            "delta_1": None,
            "delta_2": None,
            "delta_3": None,
            "conformer_selected": conformer_selected,
            # Batch Tracking
            "batch_id": batch_id,
            "batch_order": batch_order,
            "reruns": result.reruns if hasattr(result, "reruns") else 0,  # NEW
            # Timestamp
            "timestamp": (
                result.timestamp.strftime("%d/%m-%H:%M")
                if result.timestamp is not None
                else None
            ),
        }

    def get_status_counts(self) -> dict[str, int]:
        """Get count of molecules by status.

        Returns:
            Dict mapping status to count
        """
        df = self._ensure_loaded()
        counts = df["status"].value_counts().to_dict()
        return {str(k): int(v) for k, v in counts.items()}

    def get_batch_molecules(self, batch_id: str) -> pd.DataFrame:
        """Get all molecules for a batch.

        Args:
            batch_id: Batch identifier

        Returns:
            DataFrame filtered to batch
        """
        df = self._ensure_loaded()
        return df[df["batch_id"] == batch_id].copy()
