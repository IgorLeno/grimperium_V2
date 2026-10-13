"""One-time CSV migration to the two-protocol layout.

* every ``OK`` row without a protocol is preliminary protocol-1 work and
  gets ``protocol = "1"`` (its values stay exactly as they are);
* every ``Selected`` row goes back to ``Pending`` so protocol 2 recomputes it.

The default is a dry run. Applying first copies the CSV next to itself, then
saves through ``BatchCSVManager`` (atomic write).
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from semi_imperium.training.protocol import LEGACY_PROTOCOL

PROTOCOL_COLUMN = "protocol"


@dataclass(frozen=True)
class MigrationReport:
    """What the migration found and, when applied, did."""

    rows: int
    legacy_marked: int
    selected_reset: int
    already_marked: int
    applied: bool
    backup_path: Path | None


def migrate(csv_path: Path, *, apply: bool = False) -> MigrationReport:
    """Report (and with ``apply=True`` perform) the migration of ``csv_path``."""
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    manager = BatchCSVManager(csv_path)
    frame = manager.load_csv()
    if PROTOCOL_COLUMN not in frame.columns:
        frame[PROTOCOL_COLUMN] = ""
    frame[PROTOCOL_COLUMN] = frame[PROTOCOL_COLUMN].fillna("").astype(str)

    ok = frame["status"] == MoleculeStatus.OK.value
    unmarked = ok & (frame[PROTOCOL_COLUMN] == "")
    selected = frame["status"] == MoleculeStatus.SELECTED.value
    report = MigrationReport(
        rows=len(frame),
        legacy_marked=int(unmarked.sum()),
        selected_reset=int(selected.sum()),
        already_marked=int((ok & (frame[PROTOCOL_COLUMN] != "")).sum()),
        applied=False,
        backup_path=None,
    )
    if not apply:
        return report
    if report.legacy_marked == 0 and report.selected_reset == 0:
        return MigrationReport(**{**report.__dict__, "applied": True})

    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    backup = csv_path.with_name(f"{csv_path.name}.pre-protocol-{stamp}.bak")
    shutil.copy2(csv_path, backup)
    frame.loc[unmarked, PROTOCOL_COLUMN] = LEGACY_PROTOCOL
    frame.loc[selected, "status"] = MoleculeStatus.PENDING.value
    manager.save_csv()
    return MigrationReport(
        **{**report.__dict__, "applied": True, "backup_path": backup}
    )
