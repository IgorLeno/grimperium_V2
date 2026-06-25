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

    def claim_single_molecule(self, worker_id: str) -> str | None:
        """Atomically mark one pending molecule as assigned to a worker.

        Args:
            worker_id: Identifier of the worker claiming the molecule.

        Returns:
            The claimed molecule ``mol_id`` or ``None`` if no pending molecule
            is available.
        """
        with self._claim_lock:
            df = self._ensure_loaded()
            pending = df[df["status"] == MoleculeStatus.PENDING.value]
            if pending.empty:
                return None

            row = pending.iloc[0]
            idx = row.name
            mol_id = str(row["mol_id"])
            self._assign_index(cast(int, idx), worker_id)
            self._save_csv()
            LOG.debug("Assigned %s to worker %s", mol_id, worker_id)
            return mol_id

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
        self._save_csv()
        LOG.warning("%s marked %s after failure: %s", mol_id, new_status, error_message)
        return new_status

    def get_status(self, mol_id: str) -> str:
        """Return the current operational status for one molecule."""
        df = self._ensure_loaded()
        idx = self._find_molecule_index(df, mol_id)
        return str(df.at[idx, "status"])

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
        df.at[idx, "assigned_worker"] = ""
        df.at[idx, "worker_status"] = WorkerStatus.UNASSIGNED.value
        df.at[idx, "assigned_at"] = ""

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
