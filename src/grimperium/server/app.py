"""FastAPI application for the Grimperium distributed processing server.

Design rules:
- ONE asyncio.Lock guards ALL BatchCSVManager access.
- CSV methods run via asyncio.to_thread() — never called directly in async handlers.
- Watchdog launched as asyncio.create_task() in lifespan.
- Authentication via X-Token header (disabled when api_token == "").
"""

import asyncio
import logging
from collections.abc import AsyncGenerator
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated, Any

from fastapi import Depends, FastAPI, Header, HTTPException, Request

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.server.config import ServerConfig
from grimperium.server.models import (
    AssignPayload,
    ClaimResponse,
    HeartbeatRequest,
    RegisterRequest,
    ReportFailureRequest,
    ReportSuccessRequest,
    StatusResponse,
    SyncResponse,
    SyncResultsRequest,
    WorkerInfo,
)
from grimperium.server.watchdog import make_heartbeat_registry, run_watchdog

LOG = logging.getLogger(__name__)

# ── Application state ─────────────────────────────────────────────────────────

# Populated during lifespan startup; accessed via request.app.state
_STATE_CSV_MANAGER = "csv_manager"
_STATE_LOCK = "lock"
_STATE_HEARTBEAT_REGISTRY = "heartbeat_registry"
_STATE_RUNNING_MOLECULES = "running_molecules"
_STATE_CONFIG = "config"


# ── Lifespan ──────────────────────────────────────────────────────────────────


@asynccontextmanager
async def _lifespan(app: FastAPI) -> AsyncGenerator[None, None]:
    cfg: ServerConfig = app.state.config
    task = asyncio.create_task(
        run_watchdog(
            app.state.csv_manager,
            app.state.heartbeat_registry,
            cfg,
            app.state.lock,
        )
    )
    try:
        yield
    finally:
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass


# ── Authentication ─────────────────────────────────────────────────────────────


async def _verify_token(
    request: Request,
    x_token: Annotated[str | None, Header(alias="X-Token")] = None,
) -> None:
    cfg: ServerConfig = request.app.state.config
    if cfg.api_token == "":
        return
    if x_token != cfg.api_token:
        raise HTTPException(status_code=401, detail="Invalid or missing X-Token header")


AuthDep = Annotated[None, Depends(_verify_token)]


# ── Factory ───────────────────────────────────────────────────────────────────


def create_app(config: ServerConfig) -> FastAPI:
    app = FastAPI(title="Grimperium Server", lifespan=_lifespan)
    app.state.config = config

    # Initialize state synchronously so it's available before lifespan fires
    csv_manager = BatchCSVManager(Path(config.csv_path), max_reruns=config.max_reruns)
    csv_manager.load_csv()
    app.state.csv_manager = csv_manager
    app.state.lock = asyncio.Lock()
    app.state.heartbeat_registry = make_heartbeat_registry()
    app.state.running_molecules = {}  # dict[str, str]: mol_id -> worker_id

    _register_routes(app)
    return app


# ── Routes ────────────────────────────────────────────────────────────────────


def _register_routes(app: FastAPI) -> None:

    @app.post("/register")
    async def register(req: RegisterRequest, request: Request, _: AuthDep) -> dict[str, str]:
        registry: dict[str, tuple[str, datetime]] = request.app.state.heartbeat_registry
        registry[req.worker_id] = (req.hostname, datetime.now(timezone.utc))
        LOG.info("Worker %r registered from %s", req.worker_id, req.hostname)
        return {"status": "registered", "worker_id": req.worker_id}

    @app.post("/assign")
    async def assign(payload: AssignPayload, request: Request, _: AuthDep) -> dict[str, Any]:
        registry: dict[str, tuple[str, datetime]] = request.app.state.heartbeat_registry
        registry[payload.worker_id] = (
            registry.get(payload.worker_id, ("unknown", datetime.now(timezone.utc)))[0],
            datetime.now(timezone.utc),
        )
        LOG.info(
            "Assignment confirmed for worker %r: %d molecules",
            payload.worker_id,
            len(payload.molecules),
        )
        return {"status": "assigned", "count": len(payload.molecules)}

    @app.post("/claim", response_model=ClaimResponse)
    async def claim(req: HeartbeatRequest, request: Request, _: AuthDep) -> ClaimResponse:
        csv_manager: BatchCSVManager = request.app.state.csv_manager
        lock: asyncio.Lock = request.app.state.lock
        running_molecules: dict[str, str] = request.app.state.running_molecules

        async with lock:
            result = await asyncio.to_thread(csv_manager.claim_single_molecule)

        if result is None:
            return ClaimResponse(mol_id=None, smiles=None)

        mol_id, smiles = result
        running_molecules[mol_id] = req.worker_id
        LOG.info("Worker %r claimed %s", req.worker_id, mol_id)
        return ClaimResponse(mol_id=mol_id, smiles=smiles)

    @app.put("/heartbeat/{mol_id}")
    async def heartbeat(
        mol_id: str, req: HeartbeatRequest, request: Request, _: AuthDep
    ) -> dict[str, str]:
        running_molecules: dict[str, str] = request.app.state.running_molecules
        if mol_id not in running_molecules:
            raise HTTPException(status_code=404, detail=f"{mol_id} is not RUNNING")

        registry: dict[str, tuple[str, datetime]] = request.app.state.heartbeat_registry
        if req.worker_id in registry:
            hostname = registry[req.worker_id][0]
            registry[req.worker_id] = (hostname, datetime.now(timezone.utc))

        return {"status": "ok", "mol_id": mol_id}

    @app.post("/report/success")
    async def report_success(
        req: ReportSuccessRequest, request: Request, _: AuthDep
    ) -> dict[str, str]:
        csv_manager: BatchCSVManager = request.app.state.csv_manager
        lock: asyncio.Lock = request.app.state.lock
        running_molecules: dict[str, str] = request.app.state.running_molecules

        if req.mol_id not in running_molecules:
            raise HTTPException(status_code=404, detail=f"{req.mol_id} is not RUNNING")

        try:
            async with lock:
                await asyncio.to_thread(
                    csv_manager.mark_success, req.mol_id, req.result_update
                )
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=f"mol_id not found: {req.mol_id}") from exc

        running_molecules.pop(req.mol_id, None)
        LOG.info("Worker %r reported success for %s", req.worker_id, req.mol_id)
        return {"status": "ok"}

    @app.post("/report/failure")
    async def report_failure(
        req: ReportFailureRequest, request: Request, _: AuthDep
    ) -> dict[str, str]:
        csv_manager: BatchCSVManager = request.app.state.csv_manager
        lock: asyncio.Lock = request.app.state.lock
        running_molecules: dict[str, str] = request.app.state.running_molecules

        if req.mol_id not in running_molecules:
            raise HTTPException(status_code=404, detail=f"{req.mol_id} is not RUNNING")

        try:
            async with lock:
                if req.force_skip:
                    await asyncio.to_thread(
                        csv_manager.mark_skip, req.mol_id, req.error
                    )
                else:
                    await asyncio.to_thread(
                        csv_manager.mark_rerun, req.mol_id, req.error
                    )
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=f"mol_id not found: {req.mol_id}") from exc

        running_molecules.pop(req.mol_id, None)
        LOG.info(
            "Worker %r reported failure for %s (force_skip=%s)",
            req.worker_id,
            req.mol_id,
            req.force_skip,
        )
        return {"status": "ok"}

    @app.post("/sync_results", response_model=SyncResponse)
    async def sync_results(
        req: SyncResultsRequest, request: Request, _: AuthDep
    ) -> SyncResponse:
        csv_manager: BatchCSVManager = request.app.state.csv_manager
        lock: asyncio.Lock = request.app.state.lock
        running_molecules: dict[str, str] = request.app.state.running_molecules

        accepted = 0
        rejected = 0

        for result in req.results:
            try:
                async with lock:
                    if result.success:
                        await asyncio.to_thread(
                            csv_manager.mark_success,
                            result.mol_id,
                            result.result_update or {},
                        )
                    else:
                        await asyncio.to_thread(
                            csv_manager.mark_rerun,
                            result.mol_id,
                            result.error or "sync failure",
                        )
                running_molecules.pop(result.mol_id, None)
                accepted += 1
            except (KeyError, Exception) as exc:
                LOG.warning("sync_results rejected %s: %s", result.mol_id, exc)
                rejected += 1

        return SyncResponse(accepted=accepted, rejected=rejected)

    @app.get("/status", response_model=StatusResponse)
    async def status(request: Request, _: AuthDep) -> StatusResponse:
        csv_manager: BatchCSVManager = request.app.state.csv_manager
        lock: asyncio.Lock = request.app.state.lock
        registry: dict[str, tuple[str, datetime]] = request.app.state.heartbeat_registry
        running_molecules: dict[str, str] = request.app.state.running_molecules

        async with lock:
            counts = await asyncio.to_thread(csv_manager.get_status_counts)

        # Build worker list — invert running_molecules for current mol
        worker_current: dict[str, str | None] = {}
        for mol_id, wid in running_molecules.items():
            worker_current[wid] = mol_id

        workers = [
            WorkerInfo(
                worker_id=wid,
                hostname=hostname,
                last_seen=last_seen.isoformat(),
                mol_id_current=worker_current.get(wid),
            )
            for wid, (hostname, last_seen) in registry.items()
        ]

        return StatusResponse(counts=counts, workers=workers)
