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

from grimperium.cli.settings_manager import DistributedDefaults, SettingsManager
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.server.config import ServerConfig
from grimperium.server.models import (
    AssignPayload,
    ClaimResponse,
    ConfigurePayload,
    HeartbeatRequest,
    RegisterRequest,
    RegisterResponse,
    ReportFailureRequest,
    ReportSuccessRequest,
    ShutdownResponse,
    StatusResponse,
    SyncResponse,
    SyncResultsRequest,
    WorkerInfo,
    WorkerInfoExtended,
)
from grimperium.server.watchdog import make_heartbeat_registry, run_watchdog
from grimperium.server.worker_registry import WorkerRegistry, make_worker_registry

LOG = logging.getLogger(__name__)

# ── Application state ─────────────────────────────────────────────────────────

# Populated during lifespan startup; accessed via request.app.state
_STATE_CSV_MANAGER = "csv_manager"
_STATE_LOCK = "lock"
_STATE_HEARTBEAT_REGISTRY = "heartbeat_registry"
_STATE_RUNNING_MOLECULES = "running_molecules"
_STATE_CONFIG = "config"
_STATE_WORKER_REGISTRY = "worker_registry"
_STATE_DISPATCH_ENABLED = "dispatch_enabled"
_STATE_DISTRIBUTED_DEFAULTS = "distributed_defaults"


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
    app.state.worker_registry = make_worker_registry()
    app.state.dispatch_enabled = False
    try:
        app.state.distributed_defaults = SettingsManager.load_distributed_defaults()
    except RuntimeError:
        app.state.distributed_defaults = DistributedDefaults()

    _register_routes(app)
    return app


# ── Routes ────────────────────────────────────────────────────────────────────


def _register_routes(app: FastAPI) -> None:

    @app.post("/register", response_model=RegisterResponse)
    async def register(
        req: RegisterRequest, request: Request, _: AuthDep
    ) -> RegisterResponse:
        registry: dict[str, tuple[str, datetime]] = request.app.state.heartbeat_registry
        registry[req.worker_id] = (req.hostname, datetime.now(timezone.utc))
        worker_reg: WorkerRegistry = request.app.state.worker_registry
        worker_reg.register(req.worker_id, req.hostname)
        defaults: DistributedDefaults = request.app.state.distributed_defaults
        # Per-worker config override supersedes defaults (set after /configure)
        cfg_override = worker_reg.get_config(req.worker_id) or {}
        LOG.info("Worker %r registered from %s", req.worker_id, req.hostname)
        return RegisterResponse(
            worker_id=req.worker_id,
            hostname=req.hostname,
            crest_timeout_minutes=cfg_override.get(
                "crest_timeout_minutes", defaults.crest_timeout_minutes
            ),
            mopac_timeout_minutes=cfg_override.get(
                "mopac_timeout_minutes", defaults.mopac_timeout_minutes
            ),
            batch_size=cfg_override.get("batch_size", defaults.batch_size),
            profile_name=cfg_override.get("profile_name", defaults.profile_name),
        )

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
        worker_reg: WorkerRegistry = request.app.state.worker_registry

        if not request.app.state.dispatch_enabled:
            return ClaimResponse(mol_id=None, smiles=None)

        async with lock:
            result = await asyncio.to_thread(csv_manager.claim_single_molecule)

        if result is None:
            return ClaimResponse(mol_id=None, smiles=None)

        mol_id, smiles = result
        running_molecules[mol_id] = req.worker_id
        worker_reg.set_current_mol(req.worker_id, mol_id)
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
        worker_reg: WorkerRegistry = request.app.state.worker_registry
        worker_reg.record_success(req.worker_id)
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
        worker_reg: WorkerRegistry = request.app.state.worker_registry
        if req.force_skip:
            worker_reg.record_skip(req.worker_id)
        else:
            worker_reg.record_failure(req.worker_id)
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

    # ── New distributed endpoints ──────────────────────────────────────────────

    @app.get("/workers", response_model=list[WorkerInfo])
    async def list_workers(request: Request, _: AuthDep) -> list[WorkerInfo]:
        worker_reg: WorkerRegistry = request.app.state.worker_registry
        return [
            WorkerInfo(
                worker_id=e.worker_id,
                hostname=e.hostname,
                last_seen=e.last_seen.isoformat(),
                mol_id_current=e.current_mol_id,
            )
            for e in worker_reg.all_workers()
        ]

    @app.get("/workers/status", response_model=list[WorkerInfoExtended])
    async def workers_status(request: Request, _: AuthDep) -> list[WorkerInfoExtended]:
        worker_reg: WorkerRegistry = request.app.state.worker_registry
        return [
            WorkerInfoExtended(
                worker_id=e.worker_id,
                hostname=e.hostname,
                last_seen=e.last_seen.isoformat(),
                registered_at=e.registered_at.isoformat(),
                current_mol_id=e.current_mol_id,
                processed=e.processed,
                successful=e.successful,
                failed=e.failed,
                skipped=e.skipped,
                shutdown_requested=e.shutdown_requested,
            )
            for e in worker_reg.all_workers()
        ]

    @app.post("/configure/{worker_id}", response_model=dict[str, str])
    async def configure_worker(
        worker_id: str,
        payload: ConfigurePayload,
        request: Request,
        _: AuthDep,
    ) -> dict[str, str]:
        worker_reg: WorkerRegistry = request.app.state.worker_registry
        cfg = {
            "batch_size": payload.batch_size,
            "crest_timeout_minutes": payload.crest_timeout_minutes,
            "mopac_timeout_minutes": payload.mopac_timeout_minutes,
            "profile_name": payload.profile_name,
        }
        if not worker_reg.set_config(worker_id, cfg):
            raise HTTPException(status_code=404, detail=f"Worker {worker_id!r} not registered")
        LOG.info("Configured worker %r: %s", worker_id, cfg)
        return {"status": "configured", "worker_id": worker_id}

    @app.get("/configure/{worker_id}", response_model=dict[str, Any])
    async def get_worker_config(
        worker_id: str, request: Request, _: AuthDep
    ) -> dict[str, Any]:
        worker_reg: WorkerRegistry = request.app.state.worker_registry
        cfg = worker_reg.get_config(worker_id)
        if cfg is None:
            # Return distributed defaults if no override set
            defaults: DistributedDefaults = request.app.state.distributed_defaults
            return {
                "worker_id": worker_id,
                "batch_size": defaults.batch_size,
                "crest_timeout_minutes": defaults.crest_timeout_minutes,
                "mopac_timeout_minutes": defaults.mopac_timeout_minutes,
                "profile_name": defaults.profile_name,
            }
        return {"worker_id": worker_id, **cfg}

    @app.post("/shutdown/all", response_model=ShutdownResponse)
    async def shutdown_all(request: Request, _: AuthDep) -> ShutdownResponse:
        worker_reg: WorkerRegistry = request.app.state.worker_registry
        signalled = worker_reg.request_shutdown_all()
        LOG.info("Shutdown requested for all workers: %s", signalled)
        return ShutdownResponse(
            status="shutdown_requested", workers_signalled=signalled
        )

    @app.post("/shutdown/{worker_id}", response_model=ShutdownResponse)
    async def shutdown_worker(
        worker_id: str, request: Request, _: AuthDep
    ) -> ShutdownResponse:
        worker_reg: WorkerRegistry = request.app.state.worker_registry
        if not worker_reg.request_shutdown(worker_id):
            raise HTTPException(status_code=404, detail=f"Worker {worker_id!r} not registered")
        LOG.info("Shutdown requested for worker %r", worker_id)
        return ShutdownResponse(status="shutdown_requested", worker_id=worker_id)

    @app.post("/dispatch/start", response_model=dict[str, str])
    async def dispatch_start(request: Request, _: AuthDep) -> dict[str, str]:
        request.app.state.dispatch_enabled = True
        LOG.info("Dispatch enabled — workers may now claim molecules")
        return {"status": "dispatch_enabled"}
