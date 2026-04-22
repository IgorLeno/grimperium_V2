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
