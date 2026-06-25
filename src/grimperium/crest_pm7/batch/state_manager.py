"""Operational state manager for split batch outputs."""

from __future__ import annotations

import logging
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, cast

import pandas as pd

from grimperium.core.csv_utils import atomic_to_csv
from grimperium.crest_pm7.batch.enums import MoleculeStatus, WorkerStatus
from grimperium.crest_pm7.batch.output_contracts import BATCH_STATE_COLUMNS
from grimperium.crest_pm7.config import PM7Config

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

    def __init__(self, state_csv_path: Path, config: PM7Config) -> None:
        """Initialize the state manager for one ``batch_state.csv`` path.

        Args:
            state_csv_path: Path to the split operational state CSV.
            config: PM7 runtime configuration retained for future PR6 callers.
        """
        self.state_csv_path = Path(state_csv_path)
        self.config = config
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
        idx = self._find_molecule_index(df, name)
        df.at[idx, "status"] = self._normalize_status(status)
        if worker_id is not None:
            df.at[idx, "assigned_worker"] = worker_id
        for column, value in (extra_fields or {}).items():
            if column not in BATCH_STATE_COLUMNS:
                raise ValueError(f"Unknown batch_state.csv column: {column}")
            df.at[idx, column] = value
        self._save_csv()

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
