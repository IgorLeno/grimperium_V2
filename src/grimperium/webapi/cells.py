"""API the cells call: register, claim, heartbeat and result.

It is a separate app from the dashboard (``create_cell_app``), served on its
own port bound to the Tailscale address, so a cell can reach these four
routes and nothing else; the dashboard stays on 127.0.0.1. The central never
opens a connection to a cell: every exchange starts here, from the cell.

Each request carries ``Authorization: Bearer <token>``; the token names the
cell, so no body ever declares a ``cell_id``. Work is handed out only to a
cell whose last registration matched the central's programs and protocol
and whose admission passed.

Serve with (address from ``tailscale ip -4``)::

    uvicorn grimperium.webapi.cells:create_cell_app --factory \\
        --host <tailscale-ip> --port 8100
"""

from __future__ import annotations

import threading
from collections.abc import Callable
from typing import TYPE_CHECKING

from fastapi import Depends, FastAPI, Header, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.models import (
    CellClaimRequest,
    CellClaimResponse,
    CellHeartbeatRequest,
    CellHeartbeatResponse,
    CellMolecule,
    CellRegisterRequest,
    CellRegisterResponse,
    CellResultRequest,
    CellResultResponse,
)

if TYPE_CHECKING:
    from semi_imperium.training.cells import Cell

ProgramsProbe = Callable[[], dict[str, str]]

#: Well under the 90 s "online" limit, so one lost heartbeat changes nothing.
HEARTBEAT_S = 30.0


def default_central_programs() -> dict[str, str]:
    """The central's CREST/xtb/MOPAC, resolved exactly as ``run`` does."""
    from semi_imperium.training.__main__ import _executables
    from semi_imperium.training.program_info import probe_programs

    return probe_programs(_executables(None))


def create_cell_app(
    paths: DashboardPaths | None = None,
    central_programs: ProgramsProbe | None = None,
) -> FastAPI:
    """Build the cell API over the same training database as the dashboard."""
    from semi_imperium.training.cells import ADMISSION_PASSED, CellRegistry
    from semi_imperium.training.db import database_path
    from semi_imperium.training.geometry_store import GeometryRecord
    from semi_imperium.training.job_state import JobState
    from semi_imperium.training.molecule_queue import (
        DEFAULT_LEASE_S,
        RERUN,
        SKIP,
        TrainingQueue,
    )
    from semi_imperium.training.program_info import program_mismatches
    from semi_imperium.training.protocol import (
        CANONICAL_PROTOCOL,
        HAMILTONIAN,
        canonical_settings,
        canonical_signature,
    )

    files = paths or DashboardPaths.default()
    db = database_path(files.semi_imperium_store)
    registry = CellRegistry(db)
    queue = TrainingQueue(db, files.cbs_csv)
    signature = canonical_signature(canonical_settings())
    probe = central_programs or default_central_programs
    probed: dict[str, str] = {}
    probe_lock = threading.Lock()

    def reference_programs() -> dict[str, str]:
        # Hashing the binaries takes a moment: once per process. A program
        # change on the central means restarting this app (and the queue).
        with probe_lock:
            if not probed:
                probed.update(probe())
            return dict(probed)

    app = FastAPI(title="Grimperium cell API", version="0.1.0")

    @app.exception_handler(RequestValidationError)
    async def invalid_body(
        _request: Request, exc: RequestValidationError
    ) -> JSONResponse:
        # The default 422 echoes each bad input back; a NaN there cannot be
        # encoded as JSON and the refusal itself would become a 500.
        return JSONResponse(
            status_code=422,
            content={
                "detail": [
                    {"loc": list(e["loc"]), "msg": e["msg"], "type": e["type"]}
                    for e in exc.errors()
                ]
            },
        )

    def authenticated(authorization: str | None = Header(None)) -> Cell:
        scheme, _, token = (authorization or "").partition(" ")
        cell = (
            registry.authenticate(token.strip()) if scheme.lower() == "bearer" else None
        )
        if cell is None:
            raise HTTPException(
                401,
                "Token de célula ausente, inválido ou revogado.",
                headers={"WWW-Authenticate": "Bearer"},
            )
        registry.touch(cell.cell_id)
        return cell

    def registered(cell: Cell = Depends(authenticated)) -> Cell:
        if cell.refusal is not None:
            raise HTTPException(409, f"Registro recusado: {cell.refusal}")
        if cell.registered_at is None:
            raise HTTPException(409, "Célula não registrada; registre antes.")
        return cell

    def admitted(cell: Cell = Depends(registered)) -> Cell:
        if cell.admission != ADMISSION_PASSED:
            state = "bloqueada" if cell.admission == "blocked" else "em admissão"
            raise HTTPException(403, f"Célula {state}; não recebe trabalho.")
        return cell

    @app.post("/api/cell/register")
    def register(
        body: CellRegisterRequest, cell: Cell = Depends(authenticated)
    ) -> CellRegisterResponse:
        mismatches = program_mismatches(reference_programs(), body.program_versions)
        if body.protocol_signature != signature:
            mismatches.append(
                f"protocol_signature: {body.protocol_signature} (central {signature})"
            )
        refusal = "; ".join(mismatches) or None
        stored = registry.record_registration(
            cell.cell_id,
            host=body.host,
            cores=body.cores,
            workers=body.workers,
            crest_threads=body.crest_threads,
            program_versions=body.program_versions,
            protocol_signature=body.protocol_signature,
            refusal=refusal,
        )
        # A cell registers when it (re)starts: whatever it held was
        # interrupted and starts over, whether or not it is admitted now.
        released = queue.release(cell.cell_id)
        if refusal is not None:
            raise HTTPException(409, f"Programas diferentes da central: {refusal}")
        return CellRegisterResponse(
            cell_id=stored.cell_id,
            name=stored.name,
            admission=stored.admission,
            released=released,
            lease_s=DEFAULT_LEASE_S,
            heartbeat_s=HEARTBEAT_S,
        )

    @app.post("/api/cell/claim")
    def claim(
        body: CellClaimRequest, cell: Cell = Depends(admitted)
    ) -> CellClaimResponse:
        # Never more than the free workers it registered: a confused cell
        # cannot drain the queue into reservations it will not compute.
        free = (cell.workers or 0) - len(queue.held_attempts(cell.cell_id))
        molecules = queue.claim(
            cell.cell_id, min(body.count, max(free, 0)), avoid_own_failures=True
        )
        return CellClaimResponse(
            molecules=[
                CellMolecule(
                    attempt_id=m.attempt_id,
                    mol_id=m.mol_id,
                    smiles=m.smiles,
                    charge=m.charge,
                    multiplicity=m.multiplicity,
                    nheavy=m.nheavy,
                )
                for m in molecules
            ],
            lease_s=DEFAULT_LEASE_S,
        )

    @app.post("/api/cell/heartbeat")
    def heartbeat(
        body: CellHeartbeatRequest, cell: Cell = Depends(registered)
    ) -> CellHeartbeatResponse:
        renewed = queue.renew(cell.cell_id, body.attempt_ids)
        held = queue.held_attempts(cell.cell_id)
        return CellHeartbeatResponse(
            renewed=renewed,
            not_held=sorted(set(body.attempt_ids) - held),
        )

    @app.post("/api/cell/result")
    def result(
        body: CellResultRequest, cell: Cell = Depends(admitted)
    ) -> CellResultResponse:
        try:
            mol_id, owner = queue.attempt_owner(body.attempt_id)
        except KeyError:
            raise HTTPException(404, "Tentativa desconhecida.") from None
        if owner != cell.cell_id:
            raise HTTPException(403, "Tentativa de outra célula.")
        computed = body.result
        if computed is None:  # a failure: the request holds exactly one of them
            reason = str(body.failure)
            outcome = queue.accept_failure(body.attempt_id, reason)
            if outcome in (RERUN, SKIP):  # a late failure stays in the attempt log
                JobState(db).record_failure(mol_id, f"[{cell.cell_id}] {reason}")
            return CellResultResponse(outcome=outcome)
        if computed.signature != signature:
            raise HTTPException(
                422, "Assinatura do cálculo diferente do protocolo 2 da central."
            )
        try:
            outcome = queue.accept_success(
                body.attempt_id,
                geometry=GeometryRecord(
                    mol_id=computed.mol_id,
                    smiles=computed.smiles,
                    hamiltonian=HAMILTONIAN,
                    protocol=CANONICAL_PROTOCOL,
                    signature=computed.signature,
                    calculation_id=computed.calculation_id,
                    conformer_index=computed.conformer_index,
                    elements=tuple(computed.elements),
                    coordinates=tuple(computed.coordinates),
                    hof_kcal_mol=computed.hof_kcal_mol,
                    lowest_frequency_cm1=computed.lowest_frequency_cm1,
                    minimum_confirmed=True,
                ),
                conformers_evaluated=computed.conformers_evaluated,
                elapsed_s=computed.elapsed_s,
                batch_id=f"protocol{CANONICAL_PROTOCOL}-{cell.cell_id}",
            )
        except ValueError as exc:
            raise HTTPException(422, f"Resultado inconsistente: {exc}") from None
        return CellResultResponse(outcome=outcome)

    return app
