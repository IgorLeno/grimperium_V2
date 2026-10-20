"""Local dashboard API: reads the datasets and launches the training routine."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from fastapi import Depends, FastAPI, HTTPException, Query, Request
from fastapi.responses import JSONResponse

from grimperium.webapi import cells, profiles, readers, training
from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.models import (
    CellCreated,
    CellNameRequest,
    CellReleased,
    CellSummary,
    DatasetProfile,
    DatasetProgress,
    ExperimentalComparison,
    ExperimentalSummary,
    Health,
    TablePage,
    TrainingStartRequest,
    TrainingStatus,
)

TableKey = Literal["exp", "cbs"]


def create_app(
    paths: DashboardPaths | None = None,
    launcher: training.Launcher | None = None,
) -> FastAPI:
    """Build the app; ``paths`` defaults to the datasets of this checkout."""
    files = paths or DashboardPaths.default()
    dispatch = launcher or training.default_launcher
    app = FastAPI(title="Grimperium dashboard API", version="0.1.0")

    @app.exception_handler(readers.DataUnavailable)
    async def data_unavailable(
        _request: Request, exc: readers.DataUnavailable
    ) -> JSONResponse:
        # 503, not 500: the server works, one of its input files does not.
        return JSONResponse(
            status_code=503,
            content={"detail": exc.reason, "path": str(exc.path)},
        )

    @app.get("/api/health")
    def health() -> Health:
        return Health(status="ok", can_dispatch=True)

    @app.get("/api/dataset/progress")
    def dataset_progress() -> DatasetProgress:
        return readers.dataset_progress(files.cbs_csv, files.semi_imperium_store)

    @app.get("/api/experimental/summary")
    def experimental_summary() -> ExperimentalSummary:
        return readers.experimental_summary(files.exp_csv)

    @app.get("/api/experimental/comparison")
    def experimental_comparison() -> ExperimentalComparison:
        return readers.experimental_comparison(files.exp_vs_cbs_csv)

    def table_path(key: TableKey) -> Path:
        return {"exp": files.exp_csv, "cbs": files.cbs_csv}[key]

    @app.get("/api/tables/{key}/profile")
    def table_profile(key: TableKey) -> DatasetProfile:
        return profiles.profile(table_path(key), key)

    @app.get("/api/tables/{key}/rows")
    def table_rows(
        key: TableKey,
        q: str = "",
        offset: int = Query(0, ge=0),
        limit: int = Query(50, ge=1, le=profiles.MAX_PAGE),
        sort: str | None = None,
        desc: bool = False,
    ) -> TablePage:
        return profiles.rows(
            table_path(key),
            key,
            query=q,
            offset=offset,
            limit=limit,
            sort=sort,
            descending=desc,
        )

    allowed_origins = {
        f"http://{host}:{port}"
        for host in ("localhost", "127.0.0.1")
        for port in (5173, 4173, 8000)
    }

    def require_local_origin(request: Request) -> None:
        """Refuse browser POSTs from pages that are not the local dashboard."""
        origin = request.headers.get("origin")
        if origin is not None and origin not in allowed_origins:
            raise HTTPException(403, "Origem não permitida.")

    @app.get("/api/training/status")
    def training_status() -> TrainingStatus:
        return training.build_status(files)

    @app.post("/api/training/start", dependencies=[Depends(require_local_origin)])
    def training_start(body: TrainingStartRequest) -> TrainingStatus:
        return training.start(files, body, dispatch)

    @app.post("/api/training/stop", dependencies=[Depends(require_local_origin)])
    def training_stop() -> TrainingStatus:
        return training.stop(files)

    # The cells themselves talk to ``cells.create_cell_app`` on another port;
    # these routes are the dashboard's view of them and stay local.
    @app.get("/api/cells")
    def cells_list() -> list[CellSummary]:
        return cells.list_cells(files)

    @app.post("/api/cells", dependencies=[Depends(require_local_origin)])
    def cells_create(body: CellNameRequest) -> CellCreated:
        return cells.create_cell(files, body.name)

    @app.post(
        "/api/cells/{cell_id}/rename", dependencies=[Depends(require_local_origin)]
    )
    def cells_rename(cell_id: str, body: CellNameRequest) -> CellSummary:
        return cells.rename_cell(files, cell_id, body.name)

    @app.post(
        "/api/cells/{cell_id}/release", dependencies=[Depends(require_local_origin)]
    )
    def cells_release(cell_id: str) -> CellReleased:
        return cells.release_cell(files, cell_id)

    @app.post(
        "/api/cells/{cell_id}/revoke", dependencies=[Depends(require_local_origin)]
    )
    def cells_revoke(cell_id: str) -> CellSummary:
        return cells.revoke_cell(files, cell_id)

    return app
