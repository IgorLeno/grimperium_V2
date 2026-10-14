"""The CSV side of the training routine.

``pm7_fila_pipeline.csv`` stays the dataset (the ML loaders read it), so the routine
writes its results there, through ``BatchCSVManager`` and from a single
process: the coordinator. Workers never touch the file.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from semi_imperium.training.migration import PROTOCOL_COLUMN
from semi_imperium.training.protocol import CANONICAL_PROTOCOL

#: Rows the routine still has to compute. ``Running`` is an orphan of an
#: interrupted run (the coordinator waits for in-flight work before exiting).
QUEUE_STATUSES = (
    MoleculeStatus.PENDING.value,
    MoleculeStatus.SELECTED.value,
    MoleculeStatus.RERUN.value,
    MoleculeStatus.RUNNING.value,
)


class NotMigratedError(RuntimeError):
    """The CSV has no ``protocol`` column or has unmarked ``OK`` rows."""


@dataclass(frozen=True)
class QueuedMolecule:
    """One molecule waiting for protocol 2."""

    mol_id: str
    smiles: str
    charge: int
    multiplicity: int
    nheavy: int
    h298_cbs: float | None


@dataclass(frozen=True)
class TrainingProgress:
    """Row counts per stage, with the denominator (``total``) explicit."""

    total: int
    protocol2: int
    protocol1: int
    ok_unmarked: int
    pending: int
    selected: int
    rerun: int
    running: int
    skipped: int
    migrated: bool


def read_progress(csv_path: Path) -> TrainingProgress:
    """Count rows by stage without loading the full CSV or touching it."""
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    header = pd.read_csv(csv_path, nrows=0).columns
    has_protocol = PROTOCOL_COLUMN in header
    columns = ["status", *([PROTOCOL_COLUMN] if has_protocol else [])]
    frame = pd.read_csv(csv_path, usecols=columns, dtype=str, keep_default_na=False)
    status = frame["status"].str.strip()
    if has_protocol:
        protocol = (
            frame[PROTOCOL_COLUMN].str.strip().str.replace(r"\.0$", "", regex=True)
        )
    else:
        protocol = pd.Series("", index=frame.index)
    ok = status == MoleculeStatus.OK.value
    unmarked = int((ok & ~protocol.isin(["1", "2"])).sum())
    return TrainingProgress(
        total=len(frame),
        protocol2=int((ok & (protocol == "2")).sum()),
        protocol1=int((ok & (protocol == "1")).sum()),
        ok_unmarked=unmarked,
        pending=int((status == MoleculeStatus.PENDING.value).sum()),
        selected=int((status == MoleculeStatus.SELECTED.value).sum()),
        rerun=int((status == MoleculeStatus.RERUN.value).sum()),
        running=int((status == MoleculeStatus.RUNNING.value).sum()),
        skipped=int((status == MoleculeStatus.SKIP.value).sum()),
        migrated=has_protocol and unmarked == 0,
    )


class TrainingLedger:
    """Queue, state transitions and result columns of the training CSV."""

    def __init__(self, csv_path: Path, *, max_reruns: int = 3) -> None:
        self._manager = BatchCSVManager(csv_path, max_reruns=max_reruns)

    def load(self) -> None:
        """Load the CSV and insist that it was migrated."""
        frame = self._manager.load_csv()
        if PROTOCOL_COLUMN not in frame.columns:
            raise NotMigratedError(
                "pm7_fila_pipeline.csv has no 'protocol' column; run "
                "'python -m semi_imperium.training migrate --apply' first"
            )
        frame[PROTOCOL_COLUMN] = frame[PROTOCOL_COLUMN].fillna("").astype(str)
        unmarked = (frame["status"] == MoleculeStatus.OK.value) & (
            frame[PROTOCOL_COLUMN] == ""
        )
        if unmarked.any():
            raise NotMigratedError(
                f"{int(unmarked.sum())} 'OK' rows have no protocol mark; run "
                "'python -m semi_imperium.training migrate --apply' first"
            )

    def _frame(self) -> pd.DataFrame:
        if self._manager.df is None:
            raise RuntimeError("TrainingLedger.load() must be called first")
        return self._manager.df

    def queue(self, limit: int | None = None) -> list[QueuedMolecule]:
        """Queued molecules, smallest first (ties by ``mol_id``)."""
        frame = self._frame()
        rows = frame.loc[frame["status"].isin(QUEUE_STATUSES)].copy()
        rows["_size"] = pd.to_numeric(rows["nheavy"], errors="coerce").fillna(1e9)
        rows = rows.sort_values(["_size", "mol_id"], kind="stable")
        queue = [
            QueuedMolecule(
                mol_id=str(item["mol_id"]),
                smiles=str(item["smiles"]),
                charge=_as_int(item.get("charge"), 0),
                multiplicity=_as_int(item.get("multiplicity"), 1),
                nheavy=_as_int(item.get("nheavy"), 0),
                h298_cbs=_as_float(item.get("H298_cbs")),
            )
            for _, item in rows.iterrows()
        ]
        return queue if limit is None else queue[:limit]

    def mark_running(self, mol_id: str) -> None:
        """Flag ``mol_id`` as in flight (persisted immediately)."""
        self._manager.apply_operational_status(mol_id, MoleculeStatus.RUNNING.value)

    def record_success(
        self,
        mol_id: str,
        *,
        hof_kcal_mol: float,
        conformers_evaluated: int,
        elapsed_s: float,
        batch_id: str,
    ) -> None:
        """Write a protocol-2 result: value, delta to CBS, provenance columns."""
        cbs = self._cell(mol_id, "H298_cbs")
        delta = None if cbs is None else cbs - hof_kcal_mol
        pct = (
            None
            if cbs is None or cbs == 0
            else round(abs(cbs - hof_kcal_mol) / abs(cbs) * 100, 2)
        )
        update: dict[str, Any] = {
            PROTOCOL_COLUMN: CANONICAL_PROTOCOL,
            "H298_pm7": hof_kcal_mol,
            "target_delta_kcalmol": delta,
            "abs_diff_%": pct,
            "crest_status": "SUCCESS",
            "mopac_status": "OK",
            "num_conformers_selected": conformers_evaluated,
            "total_time": round(elapsed_s / 60.0, 2),
            "timestamp": datetime.now().strftime("%m/%d-%H:%M"),
            "batch_id": batch_id,
        }
        self._manager.apply_operational_status(
            mol_id, MoleculeStatus.OK.value, result_update=update
        )

    def record_failure(self, mol_id: str, reason: str) -> str:
        """Count a failed attempt: ``Rerun`` below the limit, else ``Skip``."""
        del reason  # the reason is persisted by JobState.record_failure
        reruns = _as_int(self._cell(mol_id, "reruns"), 0) + 1
        status = (
            MoleculeStatus.SKIP.value
            if reruns >= self._manager.max_reruns
            else MoleculeStatus.RERUN.value
        )
        self._manager.apply_operational_status(mol_id, status, reruns=reruns)
        return status

    def _cell(self, mol_id: str, column: str) -> float | None:
        frame = self._frame()
        match = frame.loc[frame["mol_id"] == mol_id, column]
        return None if match.empty else _as_float(match.iloc[0])


def _as_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def _as_int(value: Any, default: int) -> int:
    number = _as_float(value)
    return default if number is None else int(number)
