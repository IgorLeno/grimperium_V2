"""FastAPI app of the local dashboard: read-only, GET endpoints only."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from fastapi import FastAPI, Query, Request
from fastapi.responses import JSONResponse

from grimperium.webapi import profiles, readers
from grimperium.webapi.catalog import routine_catalog
from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.models import (
    DatasetProfile,
    DatasetProgress,
    ExperimentalComparison,
    ExperimentalSummary,
    Health,
    RoutineCatalog,
    TablePage,
)

TableKey = Literal["exp", "cbs", "pm7"]


def create_app(paths: DashboardPaths | None = None) -> FastAPI:
    """Build the app; ``paths`` defaults to the datasets of this checkout."""
    files = paths or DashboardPaths.default()
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
        return Health(status="ok", read_only=True)

    @app.get("/api/dataset/progress")
    def dataset_progress() -> DatasetProgress:
        return readers.dataset_progress(files.pm7_csv, files.semi_imperium_store)

    @app.get("/api/experimental/summary")
    def experimental_summary() -> ExperimentalSummary:
        return readers.experimental_summary(files.exp_csv)

    @app.get("/api/experimental/comparison")
    def experimental_comparison() -> ExperimentalComparison:
        return readers.experimental_comparison(files.exp_vs_cbs_csv)

    def table_path(key: TableKey) -> Path:
        path = {"exp": files.exp_csv, "cbs": files.cbs_csv, "pm7": files.pm7_csv}[key]
        if path is None:
            raise readers.DataUnavailable(Path(key), "table path not configured")
        return path

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

    @app.get("/api/routines/catalog")
    def catalog() -> RoutineCatalog:
        return routine_catalog()

    return app
