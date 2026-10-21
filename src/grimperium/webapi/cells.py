"""API the cells call (register, claim, heartbeat, result) and the
dashboard's cell actions (list, create, rename, release, revoke).

It is a separate app from the dashboard (``create_cell_app``), served on its
own port bound to the Tailscale address, so a cell can reach these four
routes and nothing else; the dashboard stays on 127.0.0.1. The central never
opens a connection to a cell: every exchange starts here, from the cell.

Each request carries ``Authorization: Bearer <token>``; the token names the
cell, so no body ever declares a ``cell_id``. Work is handed out only to a
cell whose last registration matched the central's programs and protocol
and whose admission passed. A cell still ``pending`` claims the same way but
gets the admission test's reference molecules (attempt ids ``adm-…``); their
results become the admission verdict, never queue results.

Serve with (address from ``tailscale ip -4``)::

    uvicorn grimperium.webapi.cells:create_cell_app --factory \\
        --host <tailscale-ip> --port 8100
"""

from __future__ import annotations

import tempfile
import threading
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from fastapi import Depends, FastAPI, Header, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.models import (
    CellClaimRequest,
    CellClaimResponse,
    CellCreated,
    CellHeartbeatRequest,
    CellHeartbeatResponse,
    CellMolecule,
    CellRegisterRequest,
    CellRegisterResponse,
    CellReleased,
    CellResultRequest,
    CellResultResponse,
    CellSummary,
)

if TYPE_CHECKING:
    from grimperium.webapi.cell_package import WheelBuilder
    from semi_imperium.training.cells import AdmissionMolecule, Cell, CellRegistry
    from semi_imperium.training.molecule_queue import ClaimedMolecule

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
    from semi_imperium.training.cells import (
        ADMISSION_BLOCKED,
        ADMISSION_PASSED,
        ADMISSION_PENDING,
        ADMISSION_PREFIX,
        AdmissionTest,
        AdmissionUnavailable,
        CellRegistry,
    )
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
    admission_test = AdmissionTest(db)
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

    def require_admitted(cell: Cell) -> None:
        if cell.admission != ADMISSION_PASSED:
            state = (
                "bloqueada" if cell.admission == ADMISSION_BLOCKED else "em admissão"
            )
            raise HTTPException(403, f"Célula {state}; não recebe trabalho.")

    def held(cell_id: str) -> set[str]:
        return queue.held_attempts(cell_id) | admission_test.held(cell_id)

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
        admission_test.release(cell.cell_id)
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
        body: CellClaimRequest, cell: Cell = Depends(registered)
    ) -> CellClaimResponse:
        if cell.admission != ADMISSION_PENDING:
            require_admitted(cell)
        # Never more than the free workers it registered: a confused cell
        # cannot drain the queue into reservations it will not compute.
        free = (cell.workers or 0) - len(held(cell.cell_id))
        count = min(body.count, max(free, 0))
        molecules: Sequence[ClaimedMolecule | AdmissionMolecule]
        if cell.admission == ADMISSION_PENDING:
            try:
                molecules = admission_test.claim(cell.cell_id, count)
            except AdmissionUnavailable as exc:
                raise HTTPException(
                    503, f"Teste de admissão indisponível: {exc}"
                ) from None
        else:
            molecules = queue.claim(cell.cell_id, count, avoid_own_failures=True)
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
        # Admission attempts have no lease: holding them is enough.
        holding = held(cell.cell_id)
        renewed += len(admission_test.held(cell.cell_id) & set(body.attempt_ids))
        return CellHeartbeatResponse(
            renewed=renewed,
            not_held=sorted(set(body.attempt_ids) - holding),
        )

    def admission_result(cell: Cell, body: CellResultRequest) -> CellResultResponse:
        try:
            mol_id, owner = admission_test.owner(body.attempt_id)
        except KeyError:
            raise HTTPException(404, "Tentativa desconhecida.") from None
        if owner != cell.cell_id:
            raise HTTPException(403, "Tentativa de outra célula.")
        computed = body.result
        if computed is not None:
            if computed.signature != signature:
                raise HTTPException(
                    422, "Assinatura do cálculo diferente do protocolo 2 da central."
                )
            if computed.mol_id != mol_id:
                raise HTTPException(
                    422, f"Resultado de {computed.mol_id} para a tentativa de {mol_id}."
                )
        outcome = admission_test.record(
            body.attempt_id,
            hof_kcal_mol=None if computed is None else computed.hof_kcal_mol,
            failure=None if computed is not None else str(body.failure),
        )
        return CellResultResponse(outcome=outcome)

    @app.post("/api/cell/result")
    def result(
        body: CellResultRequest, cell: Cell = Depends(registered)
    ) -> CellResultResponse:
        if body.attempt_id.startswith(ADMISSION_PREFIX):
            return admission_result(cell, body)
        require_admitted(cell)
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


# ------------------------------------------------------------ dashboard side
# Called by ``grimperium.webapi.app`` (127.0.0.1 only), never by a cell.


def _summary(cell: Cell, held: int, admission_note: str | None) -> CellSummary:
    from datetime import datetime, timezone

    from semi_imperium.training.cells import cell_state

    return CellSummary(
        cell_id=cell.cell_id,
        name=cell.name,
        state=cell_state(cell, datetime.now(timezone.utc)),
        admission=cell.admission,
        admission_note=admission_note,
        host=cell.host,
        cores=cell.cores,
        workers=cell.workers,
        crest_threads=cell.crest_threads,
        program_versions=cell.program_versions,
        refusal=cell.refusal,
        created_at=cell.created_at,
        registered_at=cell.registered_at,
        last_seen_at=cell.last_seen_at,
        revoked_at=cell.revoked_at,
        held=held,
    )


def _held(paths: DashboardPaths, cell_id: str) -> int:
    from semi_imperium.training.cells import AdmissionTest
    from semi_imperium.training.db import database_path
    from semi_imperium.training.molecule_queue import TrainingQueue

    db = database_path(paths.semi_imperium_store)
    queue = TrainingQueue(db, paths.cbs_csv)
    return len(queue.held_attempts(cell_id)) + len(AdmissionTest(db).held(cell_id))


def _row(paths: DashboardPaths, cell: Cell) -> CellSummary:
    """``cell`` as the dashboard lists it, with its holdings and admission test."""
    from semi_imperium.training.cells import AdmissionTest
    from semi_imperium.training.db import database_path

    test = AdmissionTest(database_path(paths.semi_imperium_store))
    return _summary(
        cell, _held(paths, cell.cell_id), test.note(cell.cell_id, cell.admission)
    )


def _registry(paths: DashboardPaths) -> CellRegistry:
    from semi_imperium.training.cells import CellRegistry
    from semi_imperium.training.db import database_path

    return CellRegistry(database_path(paths.semi_imperium_store))


def _known(paths: DashboardPaths, cell_id: str) -> Cell:
    try:
        return _registry(paths).get(cell_id)
    except KeyError:
        raise HTTPException(404, f"Célula desconhecida: {cell_id}") from None


def list_cells(paths: DashboardPaths) -> list[CellSummary]:
    """Every cell with its state by the central's clock; never creates the DB."""
    return [_row(paths, c) for c in _registry(paths).list()]


def create_cell(paths: DashboardPaths, name: str) -> CellCreated:
    cell, token = _registry(paths).create(name)
    return CellCreated(cell=_summary(cell, 0, None), token=token)


def rename_cell(paths: DashboardPaths, cell_id: str, name: str) -> CellSummary:
    _known(paths, cell_id)
    cell = _registry(paths).rename(cell_id, name)
    return _row(paths, cell)


def release_cell(paths: DashboardPaths, cell_id: str) -> CellReleased:
    """ "Devolver à fila": the only way a silent cell's molecules come back."""
    from semi_imperium.training.db import database_path
    from semi_imperium.training.molecule_queue import TrainingQueue

    _known(paths, cell_id)
    queue = TrainingQueue(database_path(paths.semi_imperium_store), paths.cbs_csv)
    return CellReleased(released=queue.release(cell_id))


def revoke_cell(paths: DashboardPaths, cell_id: str) -> CellSummary:
    """Invalidate the token. Held molecules stay reserved until released."""
    _known(paths, cell_id)
    cell = _registry(paths).revoke(cell_id)
    return _row(paths, cell)


def set_admission(paths: DashboardPaths, cell_id: str, admission: str) -> CellSummary:
    """Admit or block a cell by hand, or start its admission test over (pending)."""
    _known(paths, cell_id)
    cell = _registry(paths).set_admission(cell_id, admission)
    return _row(paths, cell)


def package_cell(
    paths: DashboardPaths,
    cell_id: str,
    central_url: str,
    *,
    build_wheel: WheelBuilder,
    central_programs: ProgramsProbe,
) -> tuple[str, bytes]:
    """Zip name and bytes of a fresh package for ``cell_id``, with a new token.

    The wheel and the central's programs come first: only once the package
    can be built does the token rotate (the installed cell's token then
    stops working, so the new package must replace it).
    """
    from grimperium.webapi.cell_package import (
        CellPackage,
        PackageError,
        package_folder,
    )

    cell = _known(paths, cell_id)
    if cell.revoked_at is not None:
        raise HTTPException(409, "Célula revogada; crie outra célula.")
    with tempfile.TemporaryDirectory(prefix="grimperium-celula-") as scratch:
        try:
            wheel = build_wheel(Path(scratch))
        except PackageError as exc:
            raise HTTPException(503, f"Pacote não gerado: {exc}") from None
        programs = central_programs()
        try:
            cell, token = _registry(paths).rotate_token(cell_id)
        except ValueError:
            raise HTTPException(409, "Célula revogada; crie outra célula.") from None
        data = CellPackage(
            cell_id=cell.cell_id,
            cell_name=cell.name,
            central_url=central_url.rstrip("/"),
            token=token,
            wheel=wheel,
            central_programs=programs,
        ).build()
    return f"{package_folder(cell.cell_id, cell.name)}.zip", data
