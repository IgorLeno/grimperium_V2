"""CSV manager for batch processing status tracking.

This module provides BatchCSVManager for:
- Loading and saving CSV tracking files
- Selecting molecules for batches
- Managing status transitions
- Mapping PM7Result to CSV columns

FUTURE PARALLELIZATION: Add threading.Lock() to protect DataFrame operations.
"""

import logging
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from grimperium.crest_pm7.batch.enums import (
    BatchFailurePolicy,
    BatchSortingStrategy,
    MoleculeStatus,
)
from grimperium.crest_pm7.batch.models import Batch, BatchMolecule

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

    def __init__(self, csv_path: Path) -> None:
        """Initialize CSV manager.

        Args:
            csv_path: Path to CSV tracking file
        """
        self.csv_path = Path(csv_path)
        self.df: Optional[pd.DataFrame] = None
        # FUTURE PARALLELIZATION: Add threading.Lock() here
        # self._lock = threading.Lock()

    def load_csv(self) -> pd.DataFrame:
        """Load CSV file into DataFrame.

        Returns:
            Loaded DataFrame

        Raises:
            FileNotFoundError: If CSV file doesn't exist
        """
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
            },
        )

        # Validate required columns
        required_cols = {"mol_id", "smiles", "nheavy", "status"}
        missing = required_cols - set(self.df.columns)
        if missing:
            raise ValueError(f"CSV missing required columns: {missing}")

        # Validate unique mol_ids
        if self.df["mol_id"].duplicated().any():
            duplicates = self.df[self.df["mol_id"].duplicated()]["mol_id"].tolist()
            raise ValueError(f"Duplicate mol_ids in CSV: {duplicates[:5]}...")

        LOG.info(f"Loaded {len(self.df)} molecules from {self.csv_path}")
        return self.df

    def save_csv(self) -> None:
        """Save DataFrame to CSV file."""
        if self.df is None:
            raise RuntimeError("No DataFrame loaded - call load_csv() first")

        self.df.to_csv(self.csv_path, index=False)
        LOG.debug(f"Saved {len(self.df)} molecules to {self.csv_path}")

    def _ensure_loaded(self) -> pd.DataFrame:
        """Ensure DataFrame is loaded."""
        if self.df is None:
            self.load_csv()
        return self.df  # type: ignore[return-value]

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

    def generate_batch_id(self) -> str:
        """Generate unique batch ID.

        Returns:
            Batch ID in format 'batch_YYYYMMDD_HHMMSS_xxxx'
        """
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        short_uuid = uuid.uuid4().hex[:4]
        return f"batch_{timestamp}_{short_uuid}"

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
        available_mask = df["status"].isin([
            MoleculeStatus.PENDING.value,
            MoleculeStatus.RERUN.value,
        ])
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
                nheavy=int(row["nheavy"]),
                nrotbonds=int(row.get("nrotbonds", 0)),
                retry_count=int(row.get("retry_count", 0)),
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

        Args:
            df: DataFrame of available molecules
            strategy: Sorting strategy to apply

        Returns:
            Sorted DataFrame
        """
        if strategy == BatchSortingStrategy.RERUN_FIRST_THEN_EASY:
            # RERUN first, then PENDING sorted by nheavy (ascending)
            df = df.copy()
            df["_sort_priority"] = df["status"].apply(
                lambda x: 0 if x == MoleculeStatus.RERUN.value else 1
            )
            return df.sort_values(["_sort_priority", "nheavy"]).drop(
                columns=["_sort_priority"]
            )

        elif strategy == BatchSortingStrategy.RANDOM:
            return df.sample(frac=1)

        elif strategy == BatchSortingStrategy.BY_NHEAVY:
            return df.sort_values("nheavy")

        elif strategy == BatchSortingStrategy.BY_NROTBONDS:
            return df.sort_values("nrotbonds")

        else:
            LOG.warning(f"Unknown strategy {strategy}, using default order")
            return df

    def mark_running(self, mol_id: str) -> None:
        """Mark molecule as currently running.

        Transition: SELECTED -> RUNNING

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
        self.save_csv()
        LOG.debug(f"Marked {mol_id} as RUNNING")

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
        result_update: Optional[dict[str, Any]] = None,
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
            LOG.warning(
                f"mark_rerun({mol_id}): expected RUNNING, got {current_status}"
            )

        # Increment retry count
        retry_count = int(df.at[idx, "retry_count"] or 0) + 1
        max_retries = int(df.at[idx, "max_retries"] or 3)
        df.at[idx, "retry_count"] = retry_count
        df.at[idx, "last_error_message"] = error_message

        # Determine next status
        if retry_count >= max_retries:
            df.at[idx, "status"] = MoleculeStatus.SKIP.value
            LOG.warning(f"Marked {mol_id} as SKIP (retry_count={retry_count} >= max)")
        else:
            df.at[idx, "status"] = MoleculeStatus.RERUN.value
            LOG.warning(
                f"Marked {mol_id} as RERUN ({retry_count}/{max_retries}): {error_message}"
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

            LOG.debug(
                f"Reset {df.at[idx, 'mol_id']}: {current_status} -> PENDING"
            )

        self.save_csv()
        LOG.info(f"Reset {reset_count} molecules from batch {batch_id}")
        return reset_count

    def _clear_execution_results(self, idx: int) -> None:
        """Clear execution result columns for a row.

        Args:
            idx: DataFrame row index
        """
        df = self._ensure_loaded()

        result_columns = [
            "crest_status",
            "crest_conformers_generated",
            "crest_time",
            "crest_error",
            "num_conformers_selected",
            "most_stable_hof",
            "quality_grade",
            "delta_e_12",
            "delta_e_13",
            "delta_e_15",
            "success",
            "error_message",
            "total_execution_time",
            "actual_crest_timeout_used",
            "actual_mopac_timeout_used",
            "timestamp",
        ]

        for col in result_columns:
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
        return {
            # CREST Execution
            "crest_status": result.crest_status.value,
            "crest_conformers_generated": result.crest_conformers_generated,
            "crest_time": (
                round(result.crest_time, 1) if result.crest_time else None
            ),
            "crest_error": result.crest_error,
            # MOPAC Execution
            "num_conformers_selected": result.num_conformers_selected,
            "most_stable_hof": (
                round(result.most_stable_hof, 2)
                if result.most_stable_hof
                else None
            ),
            "quality_grade": result.quality_grade.value,
            "success": result.success,
            "error_message": result.error_message,
            "total_execution_time": (
                round(result.total_execution_time, 1)
                if result.total_execution_time
                else None
            ),
            "actual_crest_timeout_used": round(crest_timeout_used, 1),
            "actual_mopac_timeout_used": round(mopac_timeout_used, 1),
            # Delta-E
            "delta_e_12": (
                round(result.delta_e_12, 4) if result.delta_e_12 else None
            ),
            "delta_e_13": (
                round(result.delta_e_13, 4) if result.delta_e_13 else None
            ),
            "delta_e_15": (
                round(result.delta_e_15, 4) if result.delta_e_15 else None
            ),
            # Batch Tracking
            "batch_id": batch_id,
            "batch_order": batch_order,
            # Timestamp
            "timestamp": result.timestamp.isoformat(),
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
