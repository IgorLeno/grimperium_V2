"""Pydantic schemas for all server request/response bodies."""

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


class ReportFailureRequest(BaseModel):
    worker_id: str
    mol_id: str
    error: str
    force_skip: bool = False


class SyncResult(BaseModel):
    mol_id: str
    success: bool
    result_update: dict[str, Any] | None
    error: str | None
    completed_at: str


class SyncResultsRequest(BaseModel):
    worker_id: str
    results: list[SyncResult]


class ClaimResponse(BaseModel):
    mol_id: str | None
    smiles: str | None


class WorkerInfo(BaseModel):
    worker_id: str
    hostname: str
    last_seen: str
    mol_id_current: str | None


class StatusResponse(BaseModel):
    counts: dict[str, int]
    workers: list[WorkerInfo]


class SyncResponse(BaseModel):
    accepted: int
    rejected: int


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
