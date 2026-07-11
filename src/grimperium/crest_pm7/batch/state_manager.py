"""Operational state manager for split batch outputs."""

from __future__ import annotations

import logging
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, cast

import pandas as pd

from grimperium.core.csv_utils import atomic_to_csv
from grimperium.crest_pm7.batch.enums import (
    BatchFailurePolicy,
    MoleculeStatus,
    WorkerStatus,
)
from grimperium.crest_pm7.batch.output_contracts import BATCH_STATE_COLUMNS
from grimperium.crest_pm7.config import PM7Config
from grimperium.crest_pm7.progress import (
    CREST_STATUS_NOT_ATTEMPTED,
    CREST_STATUS_PREOPT,
    CREST_STATUS_SEARCH,
    MOPAC_STATUS_NOT_ATTEMPTED,
    MOPAC_STATUS_RUNNING,
)

LOG = logging.getLogger("grimperium.crest_pm7.batch.state_manager")


class BatchStateManager:
    """Manage operational batch state separately from scientific output.

    This manager operates exclusively on the PR6A ``batch_state.csv`` schema. It
    does not read or write the legacy ``thermo_pm7.csv`` file owned by
    ``BatchCSVManager`` during the transition.
    """

    _CSV_DTYPE: dict[str, type[Any]] = {
        "mol_id": str,
        "status": str,
        "smiles": str,
        "batch_id": str,
        "batch_failure_policy": str,
        "timestamp": str,
        "crest_status": str,
        "xtb_status": str,
        "mopac_status": str,
        "assigned_worker": str,
        "worker_status": str,
        "assigned_at": str,
        "method_id": str,
        "method_version": str,
        "method_definition_snapshot": str,
    }

    DEFAULT_MAX_RERUNS: int = 3

    def __init__(
        self,
        state_csv_path: Path,
        config: PM7Config,
        max_reruns: int = DEFAULT_MAX_RERUNS,
    ) -> None:
        """Initialize the state manager for one ``batch_state.csv`` path.

        Args:
            state_csv_path: Path to the split operational state CSV.
            config: PM7 runtime configuration retained for future PR6 callers.
            max_reruns: Maximum retry attempts before a molecule becomes Skip.
        """
        self.state_csv_path = Path(state_csv_path)
        self.config = config
        self.max_reruns = max_reruns
        self.df: pd.DataFrame | None = None
        self._claim_lock = threading.Lock()

    def claim_single_molecule(self, worker_id: str) -> tuple[str, str] | None:
        """Atomically mark one pending or rerun molecule as assigned to a worker.

        Args:
            worker_id: Identifier of the worker claiming the molecule.

        Returns:
            ``(mol_id, smiles)`` if a molecule was available, else ``None``.
        """
        with self._claim_lock:
            df = self._ensure_loaded()
            claimable = df[
                df["status"].isin(
                    [MoleculeStatus.PENDING.value, MoleculeStatus.RERUN.value]
                )
            ]
            if claimable.empty:
                return None

            row = claimable.iloc[0]
            idx = row.name
            mol_id = str(row["mol_id"])
            smiles = str(row["smiles"])
            self._assign_index(cast(int, idx), worker_id)
            self._save_csv()
            LOG.debug("Assigned %s to worker %s", mol_id, worker_id)
            return (mol_id, smiles)

    def distribute_molecules(
        self, molecule_names: list[str], worker_ids: list[str]
    ) -> dict[str, str]:
        """Distribute molecules round-robin across workers.

        Args:
            molecule_names: Molecule identifiers to assign.
            worker_ids: Worker identifiers that can receive assignments.

        Returns:
            Mapping of molecule identifier to assigned worker identifier.

        Raises:
            ValueError: If ``molecule_names`` is non-empty but no worker is
                available.
        """
        if not molecule_names:
            return {}
        if not worker_ids:
            raise ValueError("worker_ids must not be empty when molecules are given")

        df = self._ensure_loaded()
        assignments: dict[str, str] = {}
        for offset, molecule_name in enumerate(molecule_names):
            idx = self._find_molecule_index(df, molecule_name)
            worker_id = worker_ids[offset % len(worker_ids)]
            self._assign_index(idx, worker_id)
            assignments[molecule_name] = worker_id

        self._save_csv()
        LOG.info(
            "Distributed %d molecules across %d workers",
            len(assignments),
            len(worker_ids),
        )
        return assignments

    def reassign_offline_molecules(
        self, active_worker_ids: list[str], timeout_minutes: int
    ) -> list[str]:
        """Return stale assigned molecules from inactive workers to pending.

        Args:
            active_worker_ids: Worker identifiers currently considered active.
            timeout_minutes: Minimum assignment age before inactive-worker work
                is reclaimed.

        Returns:
            Molecule identifiers returned to pending state.
        """
        df = self._ensure_loaded()
        active_workers = set(active_worker_ids)
        now = datetime.now(timezone.utc)
        threshold = timedelta(minutes=timeout_minutes)
        reassigned: list[str] = []

        assigned_mask = df["status"] == MoleculeStatus.ASSIGNED.value
        for idx in df[assigned_mask].index:
            worker_id = str(df.at[idx, "assigned_worker"])
            if worker_id in active_workers:
                continue
            age = self._parse_assignment_age(df.at[idx, "assigned_at"], now)
            if age is not None and age <= threshold:
                continue
            mol_id = str(df.at[idx, "mol_id"])
            self._clear_assignment(cast(int, idx))
            reassigned.append(mol_id)

        if reassigned:
            self._save_csv()
        return reassigned

    def reset_stuck_assigned(self, timeout_minutes: int) -> list[str]:
        """Reset molecules stuck in assigned state beyond the timeout.

        Args:
            timeout_minutes: Assignment age threshold in minutes.

        Returns:
            Molecule identifiers reset to pending state.
        """
        df = self._ensure_loaded()
        now = datetime.now(timezone.utc)
        threshold = timedelta(minutes=timeout_minutes)
        reset: list[str] = []

        assigned_mask = df["status"] == MoleculeStatus.ASSIGNED.value
        for idx in df[assigned_mask].index:
            age = self._parse_assignment_age(df.at[idx, "assigned_at"], now)
            if age is not None and age <= threshold:
                continue
            mol_id = str(df.at[idx, "mol_id"])
            self._clear_assignment(cast(int, idx))
            reset.append(mol_id)

        if reset:
            self._save_csv()
        return reset

    def reset_stuck_running(self) -> int:
        """Reset molecules stuck in RUNNING state back to PENDING.

        Used during startup recovery to reclaim molecules from interrupted
        processing runs in ``batch_state.csv``.

        Returns:
            Number of molecules reset to Pending.
        """
        reset = self.reset_stuck_running_molecules()
        count = len(reset)
        LOG.info("Startup recovery: reset %d RUNNING molecule(s) to Pending", count)
        return count

    def reset_stuck_running_molecules(self) -> list[str]:
        """Reset RUNNING molecules and return the affected molecule IDs."""
        df = self._ensure_loaded()
        running_mask = df["status"] == MoleculeStatus.RUNNING.value
        running_rows = df[running_mask]
        if running_rows.empty:
            return []

        reset: list[str] = []
        for idx in running_rows.index:
            row_idx = cast(int, idx)
            reset.append(str(df.at[row_idx, "mol_id"]))
            self._clear_assignment(row_idx)

        self._save_csv()
        return reset

    def mark_worker_offline(self, worker_id: str) -> int:
        """Return assigned or running molecules from an offline worker to PENDING.

        Called by the watchdog when a worker's heartbeat expires.

        Args:
            worker_id: Worker identifier whose molecules should be reclaimed.

        Returns:
            Number of molecules returned to Pending.
        """
        df = self._ensure_loaded()
        worker_mask = df["assigned_worker"] == worker_id
        active_mask = df["status"].isin(
            [MoleculeStatus.ASSIGNED.value, MoleculeStatus.RUNNING.value]
        )
        rows_to_reset = df[worker_mask & active_mask]
        if rows_to_reset.empty:
            return 0

        for idx in rows_to_reset.index:
            self._clear_assignment(cast(int, idx))

        self._save_csv()
        count = len(rows_to_reset)
        LOG.warning(
            "Worker %r offline: reset %d molecule(s) to Pending", worker_id, count
        )
        return count

    def seed_from_mol_list(self, molecules: list[dict[str, Any]]) -> int:
        """Populate ``batch_state.csv`` with molecules when the file is empty.

        Only seeds when the state file does not exist or contains no rows, so
        an existing batch is never silently overwritten.

        Args:
            molecules: Sequence of dicts with at least ``"mol_id"`` and
                ``"smiles"`` keys.  Additional keys matching
                ``BATCH_STATE_COLUMNS`` are also written.

        Returns:
            Number of molecules seeded (0 when file already had data).
        """
        if self.state_csv_path.exists():
            df_check = pd.read_csv(
                self.state_csv_path,
                dtype=self._CSV_DTYPE,
                keep_default_na=False,
            )
            if not df_check.empty:
                LOG.info(
                    "batch_state.csv already has %d row(s); skipping seed",
                    len(df_check),
                )
                return 0

        rows = [self._build_state_row(mol) for mol in molecules]

        self.df = pd.DataFrame(rows, columns=BATCH_STATE_COLUMNS)
        self.state_csv_path.parent.mkdir(parents=True, exist_ok=True)
        atomic_to_csv(self.state_csv_path, self.df)
        LOG.info("Seeded batch_state.csv with %d molecule(s)", len(rows))
        return len(rows)

    def reconcile_molecules(self, molecules: list[dict[str, Any]]) -> int:
        """Add missing molecules to ``batch_state.csv`` without touching state.

        Unlike :meth:`seed_from_mol_list`, this works whether or not the file
        already has rows. It is the operational-state entry point for bootstrap
        and recovery: it creates the file when absent, appends every molecule
        whose ``mol_id`` is not yet present, and for already-present molecules
        only fills a missing (empty) SMILES identity. Existing operational state
        — ``status``, worker assignment, timestamps, ``reruns`` and method
        metadata — is never modified, so ``Running``/``Assigned``/``OK``/
        ``Rerun``/``Skip`` rows are preserved. The operation is idempotent and
        persists atomically.

        Args:
            molecules: Sequence of dicts with at least ``"mol_id"``. Keys that
                match ``BATCH_STATE_COLUMNS`` are used when creating new rows.

        Returns:
            Number of molecules newly added.
        """
        df = self._ensure_loaded()
        existing_ids = {str(value) for value in df["mol_id"]}
        added = 0
        changed = False
        for mol in molecules:
            mol_id = str(mol.get("mol_id", ""))
            if not mol_id:
                continue
            if mol_id in existing_ids:
                if self._fill_missing_identity(df, mol_id, mol):
                    changed = True
                continue
            df.loc[len(df)] = self._build_state_row(mol)
            existing_ids.add(mol_id)
            added += 1
            changed = True

        if changed:
            self._save_csv()
        LOG.info("Reconciled batch_state.csv: added %d molecule(s)", added)
        return added

    def mark_selected(
        self,
        mol_ids: list[str],
        *,
        batch_id: str = "",
        extra_fields_by_mol: dict[str, dict[str, Any]] | None = None,
    ) -> None:
        """Mark molecules as Selected without assigning worker fields."""
        if not mol_ids:
            return
        df = self._ensure_loaded()
        fields_by_mol = extra_fields_by_mol or {}
        for mol_id in mol_ids:
            fields = fields_by_mol.get(mol_id, {})
            idx = self._ensure_molecule_index(df, mol_id, fields)
            df.at[idx, "status"] = MoleculeStatus.SELECTED.value
            if batch_id:
                df.at[idx, "batch_id"] = batch_id
            for column, value in fields.items():
                if column not in BATCH_STATE_COLUMNS:
                    raise ValueError(f"Unknown batch_state.csv column: {column}")
                df.at[idx, column] = value
            self._clear_assignment_fields(idx)
        self._save_csv()

    def mark_selected_from_batch(self, batch: Any) -> None:
        """Mirror a selected ``Batch`` into operational state."""
        fields_by_mol: dict[str, dict[str, Any]] = {}
        mol_ids: list[str] = []
        for mol in batch.molecules:
            mol_ids.append(str(mol.mol_id))
            fields_by_mol[str(mol.mol_id)] = {
                "smiles": mol.smiles,
                "nheavy": mol.nheavy,
                "charge": mol.charge,
                "multiplicity": mol.multiplicity,
                "batch_id": batch.batch_id,
                "batch_order": mol.batch_order,
                "batch_failure_policy": batch.failure_policy.value,
                "assigned_crest_timeout": batch.crest_timeout_minutes,
                "assigned_mopac_timeout": batch.mopac_timeout_minutes,
                "method_id": batch.method_id,
                "method_version": batch.method_version,
            }
        self.mark_selected(
            mol_ids,
            batch_id=batch.batch_id,
            extra_fields_by_mol=fields_by_mol,
        )

    def get_molecules_by_batch_id(self, batch_id: str) -> list[str]:
        """Return molecule identifiers currently associated with a batch."""
        df = self._ensure_loaded()
        values = df.loc[df["batch_id"] == batch_id, "mol_id"].tolist()
        return [str(value) for value in values]

    def get_active_worker_assignments(self) -> dict[str, str]:
        """Return mol_id -> worker_id for assigned or running molecules."""
        df = self._ensure_loaded()
        active_mask = df["status"].isin(
            [MoleculeStatus.ASSIGNED.value, MoleculeStatus.RUNNING.value]
        )
        assignments: dict[str, str] = {}
        for _, row in df[active_mask].iterrows():
            worker_id = str(row.get("assigned_worker", ""))
            if not worker_id:
                continue
            assignments[str(row["mol_id"])] = worker_id
        return assignments

    def get_pending_molecules(self) -> list[str]:
        """Return molecule identifiers with pending status."""
        return self.get_molecules_by_status(MoleculeStatus.PENDING.value)

    def get_molecules_by_status(self, status: str) -> list[str]:
        """Return molecule identifiers matching the given status."""
        df = self._ensure_loaded()
        status_value = self._normalize_status(status)
        values = df.loc[df["status"] == status_value, "mol_id"].tolist()
        return [str(value) for value in values]

    def update_molecule_status(
        self,
        name: str,
        status: str,
        worker_id: str | None = None,
        extra_fields: dict[str, Any] | None = None,
    ) -> None:
        """Update status and operational fields for one molecule.

        Args:
            name: Molecule identifier to update.
            status: New status value.
            worker_id: Optional worker assignment to write.
            extra_fields: Additional ``batch_state.csv`` fields to update.
        """
        df = self._ensure_loaded()
        idx = self._ensure_molecule_index(df, name, extra_fields or {})
        df.at[idx, "status"] = self._normalize_status(status)
        if worker_id is not None:
            df.at[idx, "assigned_worker"] = worker_id
        for column, value in (extra_fields or {}).items():
            if column not in BATCH_STATE_COLUMNS:
                raise ValueError(f"Unknown batch_state.csv column: {column}")
            df.at[idx, column] = value
        if df.at[idx, "status"] in self._completion_statuses():
            self._clear_assignment_fields(idx)
        self._save_csv()

    def mark_running(
        self, mol_id: str, extra_fields: dict[str, Any] | None = None
    ) -> None:
        """Mark a molecule as running in ``batch_state.csv``."""
        fields: dict[str, Any] = {
            "crest_status": CREST_STATUS_NOT_ATTEMPTED,
            "mopac_status": MOPAC_STATUS_NOT_ATTEMPTED,
        }
        fields.update(extra_fields or {})
        self.update_molecule_status(
            mol_id,
            MoleculeStatus.RUNNING.value,
            extra_fields=fields,
        )

    def mark_crest_preopt(self, mol_id: str) -> None:
        """Record that a molecule entered xTB pre-optimization."""
        self._update_operational_fields(
            mol_id,
            {"crest_status": CREST_STATUS_PREOPT},
        )

    def mark_crest_search(self, mol_id: str) -> None:
        """Record that a molecule entered CREST conformer search."""
        self._update_operational_fields(
            mol_id,
            {"crest_status": CREST_STATUS_SEARCH},
        )

    def mark_mopac_running(self, mol_id: str) -> None:
        """Record that a molecule entered the MOPAC calculation stage."""
        self._update_operational_fields(
            mol_id,
            {"mopac_status": MOPAC_STATUS_RUNNING},
        )

    def mark_success(self, mol_id: str) -> None:
        """Mark a molecule as successfully completed in operational state."""
        self.update_molecule_status(mol_id, MoleculeStatus.OK.value)

    def mark_rerun(self, mol_id: str, error_message: str) -> str:
        """Mark a failed molecule for retry or permanent skip.

        Args:
            mol_id: Molecule identifier.
            error_message: Failure message retained in logs.

        Returns:
            Final molecule status after retry accounting.
        """
        return self.record_rerun_or_skip(mol_id, error_message)

    def record_rerun_or_skip(self, mol_id: str, error_message: str) -> str:
        """Record failure and return either Rerun or Skip status."""
        df = self._ensure_loaded()
        idx = self._ensure_molecule_index(df, mol_id, {})
        reruns = self._safe_int(df.at[idx, "reruns"]) + 1
        df.at[idx, "reruns"] = reruns
        if reruns >= self.max_reruns:
            new_status = MoleculeStatus.SKIP.value
        else:
            new_status = MoleculeStatus.RERUN.value
        df.at[idx, "status"] = new_status
        self._clear_assignment_fields(idx)
        self._save_csv()
        LOG.warning("%s marked %s after failure: %s", mol_id, new_status, error_message)
        return new_status

    def get_status(self, mol_id: str) -> str:
        """Return the current operational status for one molecule."""
        df = self._ensure_loaded()
        idx = self._find_molecule_index(df, mol_id)
        return str(df.at[idx, "status"])

    def get_reruns(self, mol_id: str) -> int:
        """Return the current rerun count for one molecule."""
        df = self._ensure_loaded()
        idx = self._find_molecule_index(df, mol_id)
        return self._safe_int(df.at[idx, "reruns"])

    def snapshot_row(self, mol_id: str) -> dict[str, Any]:
        """Return a copy of one operational row for compensating rollback."""
        df = self._ensure_loaded()
        idx = self._find_molecule_index(df, mol_id)
        return {col: df.at[idx, col] for col in BATCH_STATE_COLUMNS}

    def restore_row(self, mol_id: str, snapshot: dict[str, Any]) -> None:
        """Restore a previously captured row snapshot (rollback)."""
        df = self._ensure_loaded()
        idx = self._find_molecule_index(df, mol_id)
        for column, value in snapshot.items():
            if column in BATCH_STATE_COLUMNS:
                df.at[idx, column] = value
        self._save_csv()

    def get_status_counts(self) -> dict[str, int]:
        """Return counts by operational status."""
        df = self._ensure_loaded()
        counts = df["status"].value_counts().to_dict()
        return {str(status): int(count) for status, count in counts.items()}

    def status_counts(self) -> dict[str, int]:
        """Return counts by operational status without legacy method naming."""
        return self.get_status_counts()

    def reset_batch(self, batch_id: str) -> int:
        """Reset all molecules in an all-or-nothing batch to pending state."""
        return self.reset_all_or_nothing(batch_id)

    def reset_all_or_nothing(self, batch_id: str) -> int:
        """Reset an all-or-nothing batch in ``batch_state.csv``.

        Args:
            batch_id: Batch identifier to reset.

        Returns:
            Number of rows reset to Pending.
        """
        df = self._ensure_loaded()
        mask = df["batch_id"] == batch_id
        batch_rows = df[mask]
        if batch_rows.empty:
            LOG.warning("No molecules found for batch %s", batch_id)
            return 0

        policy = str(batch_rows["batch_failure_policy"].iloc[0])
        if policy and policy != BatchFailurePolicy.ALL_OR_NOTHING.value:
            LOG.warning("Batch %s has policy %s; reset skipped", batch_id, policy)
            return 0

        reset_count = 0
        for idx in batch_rows.index:
            row_idx = cast(int, idx)
            df.at[row_idx, "status"] = MoleculeStatus.PENDING.value
            df.at[row_idx, "batch_id"] = ""
            df.at[row_idx, "batch_order"] = ""
            df.at[row_idx, "timestamp"] = ""
            df.at[row_idx, "total_time"] = ""
            df.at[row_idx, "crest_status"] = CREST_STATUS_NOT_ATTEMPTED
            df.at[row_idx, "xtb_status"] = ""
            df.at[row_idx, "mopac_status"] = MOPAC_STATUS_NOT_ATTEMPTED
            df.at[row_idx, "assigned_crest_timeout"] = ""
            df.at[row_idx, "assigned_mopac_timeout"] = ""
            df.at[row_idx, "assigned_worker"] = ""
            df.at[row_idx, "worker_status"] = WorkerStatus.UNASSIGNED.value
            df.at[row_idx, "assigned_at"] = ""
            reset_count += 1

        self._save_csv()
        LOG.info("Reset %d molecules from batch %s", reset_count, batch_id)
        return reset_count

    def _ensure_loaded(self) -> pd.DataFrame:
        """Load ``batch_state.csv`` and normalize it to the PR6A schema."""
        if self.df is not None:
            return self.df

        if not self.state_csv_path.exists():
            self.state_csv_path.parent.mkdir(parents=True, exist_ok=True)
            self.df = pd.DataFrame(columns=BATCH_STATE_COLUMNS)
            atomic_to_csv(self.state_csv_path, self.df)
            return self.df

        self.df = pd.read_csv(
            self.state_csv_path,
            dtype=self._CSV_DTYPE,
            keep_default_na=False,
        )
        for column in BATCH_STATE_COLUMNS:
            if column not in self.df.columns:
                self.df[column] = ""
        self.df = self.df.reindex(columns=BATCH_STATE_COLUMNS)
        if "status" in self.df.columns:
            self.df["status"] = self.df["status"].map(self._normalize_status)
        return self.df

    def _save_csv(self) -> None:
        """Persist the current DataFrame atomically."""
        if self.df is None:
            raise RuntimeError("No batch state has been loaded")
        atomic_to_csv(self.state_csv_path, self.df.reindex(columns=BATCH_STATE_COLUMNS))

    def _assign_index(self, idx: int, worker_id: str) -> None:
        """Assign a loaded DataFrame row to a worker."""
        df = self._ensure_loaded()
        df.at[idx, "status"] = MoleculeStatus.ASSIGNED.value
        df.at[idx, "assigned_worker"] = worker_id
        df.at[idx, "worker_status"] = WorkerStatus.ONLINE.value
        df.at[idx, "assigned_at"] = datetime.now(timezone.utc).isoformat()

    def _clear_assignment(self, idx: int) -> None:
        """Clear assignment fields and return a row to pending state."""
        df = self._ensure_loaded()
        df.at[idx, "status"] = MoleculeStatus.PENDING.value
        self._clear_assignment_fields(idx)

    def _clear_assignment_fields(self, idx: int) -> None:
        """Clear worker assignment fields without changing molecule status."""
        df = self._ensure_loaded()
        df.at[idx, "assigned_worker"] = ""
        df.at[idx, "worker_status"] = WorkerStatus.UNASSIGNED.value
        df.at[idx, "assigned_at"] = ""

    def _completion_statuses(self) -> set[str]:
        """Statuses that should not retain worker assignment fields."""
        return {
            MoleculeStatus.OK.value,
            MoleculeStatus.RERUN.value,
            MoleculeStatus.SKIP.value,
        }

    def _find_molecule_index(self, df: pd.DataFrame, name: str) -> int:
        """Return the integer DataFrame index for a molecule identifier."""
        matches = df.index[df["mol_id"] == name].tolist()
        if not matches:
            raise KeyError(f"Molecule not found in batch state: {name}")
        return int(matches[0])

    def _ensure_molecule_index(
        self, df: pd.DataFrame, name: str, initial_fields: dict[str, Any]
    ) -> int:
        """Return a molecule index, creating a row when execution first sees it."""
        matches = df.index[df["mol_id"] == name].tolist()
        if matches:
            return int(matches[0])

        row = dict.fromkeys(BATCH_STATE_COLUMNS, "")
        row["mol_id"] = name
        row["worker_status"] = WorkerStatus.UNASSIGNED.value
        row.update(
            {
                column: value
                for column, value in initial_fields.items()
                if column in BATCH_STATE_COLUMNS
            }
        )
        df.loc[len(df)] = row
        return int(df.index[-1])

    _EXPLICIT_SEED_FIELDS = frozenset(
        {"mol_id", "smiles", "status", "worker_status", "reruns"}
    )

    def _build_state_row(self, mol: dict[str, Any]) -> dict[str, Any]:
        """Build one ``batch_state.csv`` row from a seed/reconcile dict.

        Sets the core identity/defaults, then copies every remaining provided
        key that is a real ``BATCH_STATE_COLUMNS`` column (e.g. ``charge``,
        ``multiplicity``, ``method_id``). Unknown keys are dropped; absent
        fields keep their defaults.
        """
        row: dict[str, Any] = dict.fromkeys(BATCH_STATE_COLUMNS, "")
        row["mol_id"] = mol.get("mol_id", "")
        row["smiles"] = mol.get("smiles", "")
        row["status"] = mol.get("status", MoleculeStatus.PENDING.value)
        row["worker_status"] = WorkerStatus.UNASSIGNED.value
        row["reruns"] = str(mol.get("reruns", 0))
        row.update(
            {
                column: value
                for column, value in mol.items()
                if column in BATCH_STATE_COLUMNS
                and column not in self._EXPLICIT_SEED_FIELDS
            }
        )
        return row

    def _fill_missing_identity(
        self, df: pd.DataFrame, mol_id: str, mol: dict[str, Any]
    ) -> bool:
        """Fill only an empty SMILES identity on an existing row; never state.

        Returns ``True`` when a value was written.
        """
        idx = self._find_molecule_index(df, mol_id)
        smiles = mol.get("smiles", "")
        if smiles and not str(df.at[idx, "smiles"]):
            df.at[idx, "smiles"] = smiles
            return True
        return False

    def _update_operational_fields(
        self, mol_id: str, field_updates: dict[str, Any]
    ) -> None:
        """Update non-scientific operational fields for one molecule."""
        df = self._ensure_loaded()
        idx = self._ensure_molecule_index(df, mol_id, {})
        for column, value in field_updates.items():
            if column not in BATCH_STATE_COLUMNS:
                raise ValueError(f"Unknown batch_state.csv column: {column}")
            df.at[idx, column] = value
        self._save_csv()

    def _parse_assignment_age(
        self, assigned_at_raw: Any, now: datetime
    ) -> timedelta | None:
        """Return assignment age, or ``None`` when the timestamp is absent."""
        if pd.isna(assigned_at_raw) or not assigned_at_raw:
            return None
        try:
            assigned_at = datetime.fromisoformat(str(assigned_at_raw).rstrip("Z"))
        except ValueError:
            return None
        if assigned_at.tzinfo is None:
            assigned_at = assigned_at.replace(tzinfo=timezone.utc)
        return now - assigned_at

    def _normalize_status(self, status: Any) -> str:
        """Normalize status strings to known enum values when possible."""
        if pd.isna(status) or not status:
            return MoleculeStatus.PENDING.value
        status_text = str(status).strip()
        status_map = {item.value.lower(): item.value for item in MoleculeStatus}
        return status_map.get(status_text.lower(), status_text)

    def _safe_int(self, value: Any) -> int:
        """Convert a CSV value to integer, treating blanks as zero."""
        if pd.isna(value) or value == "":
            return 0
        try:
            return int(value)
        except (TypeError, ValueError):
            return 0
