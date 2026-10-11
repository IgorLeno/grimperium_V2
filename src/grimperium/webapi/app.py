"""FastAPI app of the local dashboard: read-only, GET endpoints only."""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from grimperium.webapi import readers
from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.models import (
    DatasetProgress,
    ExperimentalComparison,
    ExperimentalSummary,
    Health,
)


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

    return app
