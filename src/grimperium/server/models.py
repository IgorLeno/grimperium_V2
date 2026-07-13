"""Pydantic schemas for all server request/response bodies."""

from enum import Enum
from typing import Any

from pydantic import BaseModel


class RegisterRequest(BaseModel):
    worker_id: str
    hostname: str


class MoleculeRef(BaseModel):
    mol_id: str
    smiles: str


class AssignPayload(BaseModel):
    worker_id: str
    server_url: str
    crest_timeout_minutes: int
    mopac_timeout_minutes: int
    molecules: list[MoleculeRef]


class HeartbeatRequest(BaseModel):
    worker_id: str


class ReportSuccessRequest(BaseModel):
    worker_id: str
    mol_id: str
    result_update: dict[str, Any]
    result_id: str | None = None
    completed_at: str | None = None
    attempt_id: str | None = None


class ReportFailureRequest(BaseModel):
    worker_id: str
    mol_id: str
    error: str
    force_skip: bool = False
    result_id: str | None = None
    completed_at: str | None = None
    attempt_id: str | None = None


class SyncResult(BaseModel):
    result_id: str | None = None
    mol_id: str
    success: bool
    result_update: dict[str, Any] | None
    error: str | None
    completed_at: str
    attempt_id: str | None = None


class SyncResultsRequest(BaseModel):
    worker_id: str
    results: list[SyncResult]


class ClaimResponse(BaseModel):
    mol_id: str | None
    smiles: str | None
    attempt_id: str | None = None


class WorkerInfo(BaseModel):
    worker_id: str
    hostname: str
    last_seen: str
    mol_id_current: str | None


class StatusResponse(BaseModel):
    counts: dict[str, int]
    workers: list[WorkerInfo]


class SyncItemOutcome(str, Enum):
    """Per-item outcome for POST /sync_results."""

    APPLIED = "applied"
    DUPLICATE = "duplicate"
    REJECTED = "rejected"
    CONFLICT = "conflict"
    STALE_ATTEMPT = "stale_attempt"


class SyncItemResult(BaseModel):
    result_id: str
    mol_id: str
    status: SyncItemOutcome
    detail: str | None = None


class SyncResponse(BaseModel):
    accepted: int
    rejected: int
    duplicate: bool = False
    items: list[SyncItemResult] = []


class RegisterResponse(BaseModel):
    worker_id: str
    hostname: str
    crest_timeout_minutes: int
    mopac_timeout_minutes: int
    batch_size: int
    profile_name: str


class WorkerInfoExtended(BaseModel):
    worker_id: str
    hostname: str
    last_seen: str
    registered_at: str
    current_mol_id: str | None
    processed: int
    successful: int
    failed: int
    skipped: int
    shutdown_requested: bool


class ConfigurePayload(BaseModel):
    batch_size: int = 10
    crest_timeout_minutes: int = 60
    mopac_timeout_minutes: int = 30
    profile_name: str = "default"


class ShutdownResponse(BaseModel):
    status: str
    worker_id: str | None = None
    workers_signalled: list[str] | None = None
