"""CSV manager for batch processing status tracking.

This module provides BatchCSVManager for:
- Loading and saving CSV tracking files
- Selecting molecules for batches
- Managing status transitions
- Mapping PM7Result to CSV columns
"""

import logging
import math
import shutil
from pathlib import Path
from typing import Any, cast

import pandas as pd

from grimperium.core.csv_utils import atomic_to_csv
from grimperium.crest_pm7.batch.enums import (
    BatchFailurePolicy,
    BatchSortingStrategy,
    MoleculeStatus,
)
from grimperium.crest_pm7.batch.models import Batch, BatchMolecule
from grimperium.crest_pm7.batch.output_contracts import BATCH_STATE_COLUMNS
from grimperium.crest_pm7.progress import (
    CREST_STATUS_NOT_ATTEMPTED,
    CREST_STATUS_PREOPT,
    CREST_STATUS_SEARCH,
    MOPAC_STATUS_NOT_ATTEMPTED,
    MOPAC_STATUS_RUNNING,
)

LOG = logging.getLogger("grimperium.crest_pm7.batch.csv_manager")


class CSVCorruptedError(Exception):
    """Raised when both CSV and backup are corrupted or missing."""


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

    # ── Column group definitions (61 columns total) ──

    DEFAULT_MAX_RERUNS: int = 3

    IDENTITY_COLUMNS = ["mol_id", "status", "smiles"]

    MOLECULAR_PROPERTIES_COLUMNS = [
        "multiplicity",
        "charge",
        "nheavy",
        "H298_cbs",
        "H298_pm7",
        "target_delta_kcalmol",
        "abs_diff_%",
        "cbs_quality_flag",
    ]

    BATCH_INFO_COLUMNS = ["batch_id", "timestamp", "total_time", "reruns"]

    RDKIT_COLUMNS = [
        "rdkit_nrotbonds",
        "rdkit_tpsa",
        "rdkit_num_rings",
        "rdkit_fsp3",
        "rdkit_mol_weight",
        "rdkit_hbond_donors",
        "rdkit_hbond_acceptors",
        "rdkit_nC",
        "rdkit_nH",
        "rdkit_nO",
        "rdkit_nN",
        "rdkit_bonds_single",
        "rdkit_bonds_double",
        "rdkit_bonds_triple",
        "rdkit_bonds_aromatic",
    ]

    CREST_COLUMNS = [
        "crest_status",
        "xtb",
        "v3",
        "qm",
        "nci",
        "c_method",
        "energy_window",
        "rmsd_threshold",
        "opt_lvl",
        "crest_conformers_generated",
        "num_conformers_selected",
        "crest_time_s",
    ]

    MOPAC_COLUMNS = [
        "mopac_status",
        "k_selected_pm7",
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

    BATCH_CONFIG_COLUMNS = [
        "batch_order",
        "batch_failure_policy",
        "assigned_crest_timeout",
        "assigned_mopac_timeout",
    ]

    PREDICTION_OUTPUT_COLUMNS = ["H298_predicted", "delta_correction"]

    # Result columns that are cleared when resetting a batch
    RESULT_COLUMNS = [
        "crest_status",
        "crest_conformers_generated",
        "crest_time_s",
        "mopac_status",
        "num_conformers_selected",
        "H298_pm7",
        "abs_diff_%",
        "assigned_crest_timeout",
        "assigned_mopac_timeout",
        "target_delta_kcalmol",
        "k_selected_pm7",
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
        "timestamp",
        "total_time",
        "reruns",
    ]

    def __init__(
        self,
        csv_path: Path | None,
        max_reruns: int = DEFAULT_MAX_RERUNS,
    ) -> None:
        """Initialize CSV manager.

        Args:
            csv_path: Path to CSV tracking file (can be None for schema-only use)
            max_reruns: Maximum number of rerun attempts before a molecule is
                marked SKIP (default: DEFAULT_MAX_RERUNS = 3).
        """
        self.csv_path = Path(csv_path) if csv_path is not None else None
        self.df: pd.DataFrame | None = None
        self.max_reruns = max_reruns

    def get_schema(self) -> list[str]:
        """Get the full CSV schema with all column names.

        Returns:
            List of 61 column names in order:
            - Identity (3): mol_id, status, smiles
            - Molecular properties (8): multiplicity, charge, nheavy, H298_cbs, etc.
            - Batch info (4): batch_id, timestamp, total_time, reruns
            - RDKit descriptors (15): rdkit_nrotbonds, rdkit_tpsa, etc.
            - CREST (12): crest_status, xtb, v3, qm, etc.
            - MOPAC (13): mopac_status, k_selected_pm7, etc.
            - Batch config (4): batch_order, batch_failure_policy, etc.
            - Prediction outputs (2): H298_predicted, delta_correction
        """
        return (
            self.IDENTITY_COLUMNS
            + self.MOLECULAR_PROPERTIES_COLUMNS
            + self.BATCH_INFO_COLUMNS
            + self.RDKIT_COLUMNS
            + self.CREST_COLUMNS
            + self.MOPAC_COLUMNS
            + self.BATCH_CONFIG_COLUMNS
            + self.PREDICTION_OUTPUT_COLUMNS
        )

    def create_empty_csv(self) -> pd.DataFrame:
        """Create empty CSV with full schema on first installation.

        Creates the parent directory if needed, builds a zero-row DataFrame
        with all columns from get_schema(), and saves atomically.

        Returns:
            Empty DataFrame with full schema

        Raises:
            RuntimeError: If csv_path is None
        """
        if self.csv_path is None:
            raise RuntimeError("csv_path is None - cannot create CSV")
        self.csv_path.parent.mkdir(parents=True, exist_ok=True)
        self.df = pd.DataFrame(columns=self.get_schema())
        atomic_to_csv(self.csv_path, self.df)
        LOG.info(
            "Created empty CSV with %d columns at %s",
            len(self.get_schema()),
            self.csv_path,
        )
        return self.df

    _CSV_DTYPE: dict[str, type[Any]] = {
        "mol_id": str,
        "smiles": str,
        "status": str,
        "batch_id": str,
        "batch_failure_policy": str,
        "crest_status": str,
        "mopac_status": str,
        # Force string type for timestamp to avoid float64 issues
        "timestamp": str,
    }

    def load_csv(self) -> pd.DataFrame:
        """Load CSV file into DataFrame with automatic recovery.

        Recovery behaviour:
            1. Orphan ``.csv.tmp`` files from interrupted writes are removed.
            2. If the CSV is 0-byte or unparseable, the ``.bak`` backup is
               restored transparently.
            3. If the CSV and backup are missing, an empty CSV is created with
               the full schema.
            4. If both CSV and backup are corrupt,
               :class:`CSVCorruptedError` is raised.

        Returns:
            Loaded DataFrame

        Raises:
            RuntimeError: If csv_path is None.
            CSVCorruptedError: If CSV and backup are both unreadable.
        """
        if self.csv_path is None:
            raise RuntimeError("csv_path is None - cannot load CSV")

        # --- Clean orphan temp files from interrupted atomic writes -----------
        # Intentionally broad: atomic temp names do not encode target CSV name,
        # so we remove all orphan *.csv.tmp files in this directory.
        for tmp in self.csv_path.parent.glob("*.csv.tmp"):
            try:
                tmp.unlink()
                LOG.warning("Removed orphan temp file: %s", tmp)
            except OSError as exc:
                LOG.debug("Could not remove orphan temp %s: %s", tmp, exc)

        # --- Attempt to read, falling back to .bak on corruption -------------
        backup_path = Path(str(self.csv_path) + ".bak")

        if not self.csv_path.exists():
            # Primary missing — try backup before giving up
            if backup_path.exists() and backup_path.stat().st_size > 0:
                LOG.error("CSV missing, restoring from backup: %s", backup_path)
                shutil.copy2(backup_path, self.csv_path)
            else:
                LOG.warning(
                    "CSV not found at %s. Creating empty database with full schema.",
                    self.csv_path,
                )
                return self.create_empty_csv()

        try:
            if self.csv_path.stat().st_size == 0:
                raise pd.errors.EmptyDataError("CSV file is 0 bytes")
            self.df = pd.read_csv(self.csv_path, dtype=self._CSV_DTYPE)
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
            LOG.error("CSV unreadable (%s), checking backup...", exc)
            if backup_path.exists() and backup_path.stat().st_size > 0:
                LOG.error(
                    "Restoring CSV from backup: %s -> %s",
                    backup_path,
                    self.csv_path,
                )
                shutil.copy2(backup_path, self.csv_path)
                try:
                    self.df = pd.read_csv(self.csv_path, dtype=self._CSV_DTYPE)
                except (pd.errors.ParserError, pd.errors.EmptyDataError) as bak_exc:
                    raise CSVCorruptedError(
                        "Backup restore failed: parsed backup "
                        f"{backup_path} into {self.csv_path}, but resulting CSV "
                        "is still unreadable."
                    ) from bak_exc
                LOG.warning("Recovered %d molecules from backup", len(self.df))
            else:
                raise CSVCorruptedError(
                    f"CSV and backup both corrupted/missing. "
                    f"Restore manually: {self.csv_path}"
                ) from exc

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
        """Save DataFrame to CSV file atomically with backup.

        Uses temp-file-then-rename to prevent corruption on Ctrl-C.
        A ``.bak`` backup is created before overwriting.
        """
        if self.csv_path is None:
            raise RuntimeError("csv_path is None - cannot save CSV")
        if self.df is None:
            raise RuntimeError("No DataFrame loaded - call load_csv() first")

        state_only_columns = set(BATCH_STATE_COLUMNS) - set(self.get_schema())
        self.df = self.df.drop(
            columns=[col for col in self.df.columns if col in state_only_columns]
        )
        atomic_to_csv(self.csv_path, self.df)
        LOG.info("Atomically saved %d molecules to %s", len(self.df), self.csv_path)

    def _ensure_loaded(self) -> pd.DataFrame:
        """Ensure DataFrame is loaded."""
        if self.df is None:
            self.load_csv()
        if self.df is None:
            raise RuntimeError("load_csv() did not populate self.df")
        return self.df

    def _safe_int(self, val: Any, default: int = 0) -> int:
        """Safely convert value to int, handling NaN and invalid types.

        WARNING: Truncates decimal values (e.g., 3.9 -> 3).
        Use with caution on CSV fields that may contain floats.
        A LOG.warning is emitted when truncation occurs.
        Uses math.isclose() to avoid spurious warnings from floating-point precision.

        Args:
            val: Value to convert (Any type)
            default: Default to return if conversion fails

        Returns:
            int: Converted value (truncated if float), or default
        """
        if pd.isna(val):
            return default

        try:
            return int(val)
        except (ValueError, TypeError):
            try:
                if isinstance(val, str):
                    val = val.strip()

                float_val = float(val)
                int_val = int(float_val)

                if not math.isclose(float_val, int_val, rel_tol=1e-9, abs_tol=1e-9):
                    LOG.warning(
                        f"Truncating float {float_val} to int {int_val}. "
                        f"Original value: {val}. This may cause data loss. "
                        f"Consider validating or rounding before CSV import."
                    )

                return int_val

            except (ValueError, TypeError):
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
            if "H298_cbs" in df.columns:
                val = df.at[idx, "H298_cbs"]
            else:
                LOG.debug(f"[{mol_id}] CSV missing H298_cbs column")
                return None
            if pd.isna(val):
                return None
            return float(cast(Any, val))
        except (KeyError, ValueError, TypeError) as e:
            LOG.debug(f"[{mol_id}] Could not get H298_cbs: {e}")
            return None

    def generate_batch_id(self) -> str:
        """Generate unique batch ID with sequential numbering.

        Returns:
            Batch ID in format 'batch_NNNN' where NNNN is 4-digit counter
            (e.g., 'batch_0001', 'batch_0002', ...)
        """
        df = self._ensure_loaded()

        existing_batches = df["batch_id"].dropna().unique()
        batch_numbers = []

        for batch_id in existing_batches:
            if isinstance(batch_id, str) and batch_id.startswith("batch_"):
                try:
                    num = int(batch_id.replace("batch_", ""))
                    batch_numbers.append(num)
                except ValueError:
                    pass

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

        # ── Recover orphaned SELECTED from interrupted batches ──
        if "batch_id" not in df.columns:
            df["batch_id"] = pd.NA
        orphan_mask = (df["status"] == MoleculeStatus.SELECTED.value) & (
            df["batch_id"] != batch_id
        )
        orphan_df = df[orphan_mask].copy()
        orphan_molecules: list[BatchMolecule] = []

        if not orphan_df.empty:
            LOG.warning(
                "Found %d orphaned SELECTED molecule(s) from previous batch(es). "
                "Recovering with maximum priority.",
                len(orphan_df),
            )
            orphan_df = orphan_df.sort_values("mol_id")
            orphan_count = min(batch_size, len(orphan_df))
            orphan_df = orphan_df.head(orphan_count)

            for batch_order, (idx, row) in enumerate(orphan_df.iterrows(), start=1):
                mol = BatchMolecule(
                    mol_id=row["mol_id"],
                    smiles=row["smiles"],
                    batch_order=batch_order,
                    nheavy=self._safe_int(row["nheavy"], 0),
                    nrotbonds=self._safe_int(row.get("rdkit_nrotbonds", 0), 0),
                    reruns=self._safe_int(row.get("reruns", 0), 0),
                )
                orphan_molecules.append(mol)
                df.at[idx, "batch_id"] = batch_id
                df.at[idx, "batch_order"] = batch_order
                df.at[idx, "batch_failure_policy"] = failure_policy.value
                df.at[idx, "assigned_crest_timeout"] = crest_timeout_minutes
                df.at[idx, "assigned_mopac_timeout"] = mopac_timeout_minutes

            self.save_csv()

        remaining_capacity = batch_size - len(orphan_molecules)

        # Find available molecules (PENDING or RERUN)
        new_molecules: list[BatchMolecule] = []
        if remaining_capacity > 0:
            available_mask = df["status"].isin(
                [
                    MoleculeStatus.PENDING.value,
                    MoleculeStatus.RERUN.value,
                ]
            )
            available_df = df[available_mask].copy()

            if not available_df.empty:
                # Apply sorting strategy
                available_df = self._apply_sorting_strategy(available_df, strategy)

                actual_size = min(remaining_capacity, len(available_df))
                if actual_size < remaining_capacity and not orphan_molecules:
                    LOG.info(
                        f"Requested {batch_size} molecules, "
                        f"only {len(available_df)} available"
                    )

                selected_df = available_df.head(actual_size)
                order_offset = len(orphan_molecules)

                for batch_order, (idx, row) in enumerate(
                    selected_df.iterrows(), start=order_offset + 1
                ):
                    mol = BatchMolecule(
                        mol_id=row["mol_id"],
                        smiles=row["smiles"],
                        batch_order=batch_order,
                        nheavy=self._safe_int(row["nheavy"], 0),
                        nrotbonds=self._safe_int(row.get("rdkit_nrotbonds", 0), 0),
                        reruns=self._safe_int(row.get("reruns", 0), 0),
                    )
                    new_molecules.append(mol)

                    # Update DataFrame: mark as SELECTED
                    df.at[idx, "status"] = MoleculeStatus.SELECTED.value
                    df.at[idx, "batch_id"] = batch_id
                    df.at[idx, "batch_order"] = batch_order
                    df.at[idx, "batch_failure_policy"] = failure_policy.value
                    df.at[idx, "assigned_crest_timeout"] = crest_timeout_minutes
                    df.at[idx, "assigned_mopac_timeout"] = mopac_timeout_minutes

                self.save_csv()

        all_molecules = orphan_molecules + new_molecules

        if not all_molecules:
            LOG.warning("No pending/rerun molecules available for batch")
            return Batch(
                batch_id=batch_id,
                molecules=[],
                crest_timeout_minutes=crest_timeout_minutes,
                mopac_timeout_minutes=mopac_timeout_minutes,
                sorting_strategy=strategy,
                failure_policy=failure_policy,
            )

        LOG.info(
            f"Selected {len(all_molecules)} molecules for batch {batch_id} "
            f"(strategy={strategy.value}, policy={failure_policy.value})"
        )

        return Batch(
            batch_id=batch_id,
            molecules=all_molecules,
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
            df = df.copy()
            df["_sort_priority"] = df["status"].apply(
                lambda x: 0 if x == MoleculeStatus.RERUN.value else 1
            )
            return cast(
                pd.DataFrame,
                df.sort_values(["_sort_priority", "mol_id"]).drop(
                    columns=["_sort_priority"]
                ),
            )

        elif strategy == BatchSortingStrategy.RANDOM:
            return cast(pd.DataFrame, df.sample(frac=1))

        elif strategy == BatchSortingStrategy.BY_NHEAVY:
            return cast(pd.DataFrame, df.sort_values("nheavy"))

        elif strategy == BatchSortingStrategy.BY_NROTBONDS:
            if "rdkit_nrotbonds" in df.columns:
                return cast(pd.DataFrame, df.sort_values("rdkit_nrotbonds"))
            else:
                LOG.warning(
                    "Column 'rdkit_nrotbonds' not found, skipping BY_NROTBONDS sorting"
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
        df = self._ensure_loaded()
        idx = self._get_row_index(mol_id)
        return str(df.at[idx, "status"])

    def mark_running(self, mol_id: str) -> None:
        """Mark molecule as currently running.

        Transition: SELECTED -> RUNNING

        Also resets crest_status and mopac_status to NOT_ATTEMPTED
        to ensure progress tracking starts cleanly for the run.

        Args:
            mol_id: Molecule identifier
        """
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
            result_update: Dict with CSV column updates from pm7result_to_csv_update
        """
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

        Transition: RUNNING -> RERUN (if reruns < max allowed)
        Transition: RUNNING -> SKIP (if reruns >= max allowed)

        Args:
            mol_id: Molecule identifier
            error_message: Error that caused the failure
            result_update: Optional partial results to save
        """
        idx = self._get_row_index(mol_id)
        df = self._ensure_loaded()

        current_status = df.at[idx, "status"]
        if current_status != MoleculeStatus.RUNNING.value:
            LOG.warning(f"mark_rerun({mol_id}): expected RUNNING, got {current_status}")

        # Increment reruns count
        reruns = self._safe_int(df.at[idx, "reruns"], default=0) + 1
        max_reruns = self.max_reruns

        df.at[idx, "reruns"] = reruns

        error_base = (error_message or "").strip()
        if not error_base:
            error_base = "Unknown error"

        # Determine next status
        if reruns >= max_reruns:
            df.at[idx, "status"] = MoleculeStatus.SKIP.value
            LOG.warning(f"Marked {mol_id} as SKIP (reruns={reruns} >= max)")
        else:
            df.at[idx, "status"] = MoleculeStatus.RERUN.value
            LOG.warning(
                f"Marked {mol_id} as RERUN " f"({reruns}/{max_reruns}): {error_message}"
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
        idx = self._get_row_index(mol_id)
        df = self._ensure_loaded()

        reruns = self._safe_int(df.at[idx, "reruns"], default=0)

        error_base = (error_message or "").strip()
        if not error_base:
            error_base = "Unknown error"

        df.at[idx, "status"] = MoleculeStatus.SKIP.value

        self.save_csv()
        LOG.warning(f"Marked {mol_id} as SKIP (reruns={reruns}): {error_base}")

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
        """
        idx = self._get_row_index(mol_id)
        df = self._ensure_loaded()

        updated_count = 0
        skipped_cols = []

        for col, val in field_updates.items():
            if col in df.columns:
                df[col] = df[col].astype(object)
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

        Resets status to PENDING but keeps reruns for audit.

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

    def reset_stuck_running(self) -> int:
        """Reset RUNNING molecules on startup recovery.

        When a batch is interrupted (crash, Ctrl+C, etc.), molecules may be
        left in RUNNING status. This method recovers them by incrementing
        reruns and moving to PENDING or SKIP according to max_reruns.

        Returns:
            Number of RUNNING molecules recovered
        """
        df = self._ensure_loaded()

        running_mask = df["status"] == MoleculeStatus.RUNNING.value
        if not running_mask.any():
            return 0

        count = int(running_mask.sum())
        for idx in df[running_mask].index:
            mol_id = df.at[idx, "mol_id"]
            new_reruns = self._safe_int(df.at[idx, "reruns"], default=0) + 1
            df.at[idx, "reruns"] = new_reruns

            if new_reruns >= self.max_reruns:
                df.at[idx, "status"] = MoleculeStatus.SKIP.value
                LOG.warning(
                    f"Reset stuck RUNNING molecule {mol_id} to SKIP "
                    f"(reruns={new_reruns} >= max)"
                )
            else:
                df.at[idx, "status"] = MoleculeStatus.PENDING.value
                LOG.warning(
                    f"Reset stuck RUNNING molecule {mol_id} to PENDING "
                    f"(reruns={new_reruns}/{self.max_reruns})"
                )

        self.save_csv()
        LOG.info(f"Recovered {count} stuck RUNNING molecules")
        return count

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
        # Calculate H298_pm7
        h298_pm7 = (
            round(result.most_stable_hof, 2)
            if getattr(result, "most_stable_hof", None) is not None
            else None
        )

        # Calculate abs_diff_%
        h298_cbs = result.H298_cbs if hasattr(result, "H298_cbs") else None
        abs_diff_pct = None
        if h298_pm7 is not None and h298_cbs is not None:
            abs_diff = abs(h298_cbs - h298_pm7)
            if h298_cbs != 0:
                abs_diff_pct = round((abs_diff / abs(h298_cbs)) * 100, 2)

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
        if mopac_status is None:
            mopac_status = MOPAC_STATUS_NOT_ATTEMPTED

        # Map RDKit descriptors from result
        rdkit_updates: dict[str, Any] = {}
        rdkit_desc = getattr(result, "rdkit_descriptors", None) or {}
        for key in self.RDKIT_COLUMNS:
            rdkit_updates[key] = rdkit_desc.get(key)

        update: dict[str, Any] = {
            # RDKit descriptors
            **rdkit_updates,
            # CREST Execution
            "crest_status": (
                result.crest_status.value if result.crest_status is not None else None
            ),
            "crest_conformers_generated": result.crest_conformers_generated,
            "crest_time_s": (
                round(result.crest_time, 1) if result.crest_time is not None else None
            ),
            # MOPAC Execution
            "mopac_status": mopac_status,
            "num_conformers_selected": result.num_conformers_selected,
            "H298_pm7": h298_pm7,
            "abs_diff_%": abs_diff_pct,
            "total_time": (
                round(result.total_execution_time / 60.0, 2)
                if result.total_execution_time is not None
                else None
            ),
            "assigned_crest_timeout": round(crest_timeout_used, 1),
            "assigned_mopac_timeout": round(mopac_timeout_used, 1),
            # PM7 selection & descriptors (populated by csv_enhancements)
            "target_delta_kcalmol": None,
            "k_selected_pm7": None,
            "mopac_dipole_debye": None,
            "mopac_ionization_potential_ev": None,
            "mopac_homo_ev": None,
            "mopac_lumo_ev": None,
            "mopac_gap_ev": None,
            "mopac_cosmo_area_a2": None,
            "mopac_cosmo_volume_a3": None,
            "mopac_gradient_norm": None,
            "mopac_num_scf_cycles": None,
            "mopac_point_group": None,
            "mopac_time_s": None,
            # Batch Tracking
            "batch_id": batch_id,
            "batch_order": batch_order,
            # Timestamp
            "timestamp": (
                result.timestamp.strftime("%d/%m-%H:%M")
                if result.timestamp is not None
                else None
            ),
        }

        return update

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
        return cast(pd.DataFrame, df[df["batch_id"] == batch_id].copy())

    def claim_single_molecule(self) -> tuple[str, str] | None:
        """Atomically claim one PENDING/RERUN molecule for server-side dispatch.

        Unlike select_batch(), does NOT raise when other molecules are RUNNING.
        Designed for the distributed server where multiple workers process
        different molecules concurrently.

        Returns:
            (mol_id, smiles) if a molecule was available, else None.
        """
        df = self._ensure_loaded()
        pool_mask = df["status"].isin(
            [MoleculeStatus.PENDING.value, MoleculeStatus.RERUN.value]
        )
        available = df[pool_mask]
        if available.empty:
            return None
        row = available.iloc[0]
        idx = int(row.name)
        mol_id = str(row["mol_id"])
        smiles = str(row["smiles"])
        df.at[idx, "status"] = MoleculeStatus.RUNNING.value
        self.save_csv()
        LOG.debug("Claimed %s for distributed processing", mol_id)
        return (mol_id, smiles)
