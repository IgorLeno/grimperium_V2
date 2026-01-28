"""Shared progress-tracking constants for batch CSV and CLI progress views."""

from __future__ import annotations

from enum import Enum
from typing import Callable, Final

from grimperium.crest_pm7.batch.enums import MoleculeStatus

STATUS_PENDING: Final[str] = MoleculeStatus.PENDING.value
STATUS_SELECTED: Final[str] = MoleculeStatus.SELECTED.value
STATUS_RUNNING: Final[str] = MoleculeStatus.RUNNING.value
STATUS_OK: Final[str] = MoleculeStatus.OK.value
STATUS_RERUN: Final[str] = MoleculeStatus.RERUN.value
STATUS_SKIP: Final[str] = MoleculeStatus.SKIP.value

CREST_STATUS_NOT_ATTEMPTED: Final[str] = "NOT_ATTEMPTED"
CREST_STATUS_PREOPT: Final[str] = "XTB_PREOPT"
CREST_STATUS_SEARCH: Final[str] = "CREST_SEARCH"
CREST_STATUS_SUCCESS: Final[str] = "SUCCESS"
CREST_STATUS_FAILED: Final[str] = "FAILED"

MOPAC_STATUS_NOT_ATTEMPTED: Final[str] = "NOT_ATTEMPTED"
MOPAC_STATUS_RUNNING: Final[str] = "RUNNING"
MOPAC_STATUS_OK: Final[str] = "OK"
MOPAC_STATUS_FAILED: Final[str] = "FAILED"

COMPLETION_STATUSES: Final[set[str]] = {
    STATUS_OK,
    STATUS_RERUN,
    STATUS_SKIP,
}


class BatchProgressStage(str, Enum):
    """Internal batch processing stages for CSV progress updates."""

    XTB_PREOPT = "XTB_PREOPT"
    CREST_SEARCH = "CREST_SEARCH"
    MOPAC_CALC = "MOPAC_CALC"


BatchProgressCallback = Callable[[BatchProgressStage], None]
