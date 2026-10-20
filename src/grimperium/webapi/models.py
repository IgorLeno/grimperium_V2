"""Response bodies of the dashboard API (kcal/mol for every enthalpy)."""

from __future__ import annotations

from pydantic import BaseModel, Field

BLANK_LABEL = "unset"
"""Count key for empty cells, so a blank flag is visible, not dropped."""


class Health(BaseModel):
    status: str
    can_dispatch: bool


class TrainingProgress(BaseModel):
    """Molecules of the training queue per stage; ``total`` is the denominator."""

    total: int
    ok: int
    pending: int
    rerun: int
    running: int
    skipped: int


class StoreProgress(BaseModel):
    path: str
    exists: bool
    molecules: int
    runs: int
    calculations: int
    states: dict[str, int]


class DatasetProgress(BaseModel):
    queue: TrainingProgress
    semi_imperium: StoreProgress


class ExperimentalSummary(BaseModel):
    path: str
    rows: int
    by_source: dict[str, int]
    by_source_family: dict[str, int]
    flags: dict[str, int]
    uncertainty_kind: dict[str, int]


class ErrorStats(BaseModel):
    n: int
    mae: float | None
    bias: float | None
    rmse: float | None


class ComparisonPoint(BaseModel):
    exp_id: str
    mol_id: str
    smiles: str
    nheavy: int
    source: str
    H298_exp: float
    uncertainty: float | None
    H298_cbs: float
    cbs_minus_exp: float
    conflict_flag: bool
    low_priority_only: bool
    validation_only: bool


class ExperimentalComparison(BaseModel):
    path: str
    units: str
    stats: ErrorStats
    points: list[ComparisonPoint]


# ------------------------------------------------------- dataset profiles


class HistogramBin(BaseModel):
    start: float
    end: float
    count: int


class Histogram(BaseModel):
    """Fixed-width bins over the central 99% of a numeric column.

    Values outside ``[bins[0].start, bins[-1].end]`` are counted in
    ``below``/``above`` instead of stretching the axis, so the chart stays
    readable and no row silently disappears.
    """

    column: str
    label: str
    unit: str
    n: int
    below: int
    above: int
    mean: float | None
    median: float | None
    p05: float | None
    p95: float | None
    minimum: float | None
    maximum: float | None
    bins: list[HistogramBin]


class SizeStats(BaseModel):
    """Error column summarised per number of heavy atoms."""

    nheavy: int
    n: int
    mean: float
    mean_abs: float


class ProfilePair(BaseModel):
    id: str
    x: float
    y: float
    nheavy: int | None


class DatasetProfile(BaseModel):
    key: str
    title: str
    path: str
    rows: int
    unique_smiles: int | None
    nheavy: dict[str, int]
    composition: dict[str, int]
    categories: dict[str, dict[str, int]]
    histograms: list[Histogram]
    by_nheavy: list[SizeStats]
    pairs: list[ProfilePair]
    pairs_total: int
    missing_columns: list[str]
    manifest: dict[str, object] | None


class TablePage(BaseModel):
    key: str
    columns: list[str]
    total: int
    matched: int
    offset: int
    limit: int
    rows: list[dict[str, str | float | int | bool | None]]


# ------------------------------------------------------------------ training


class TrainingJob(BaseModel):
    state: str
    pid: int | None
    started_at: str | None
    heartbeat_at: str | None
    finished_at: str | None
    workers: int
    limit: int | None
    done: int
    failed: int
    message: str
    stop_requested: bool


class TrainingFailure(BaseModel):
    mol_id: str
    reason: str
    at: str


class TrainingInFlight(BaseModel):
    mol_id: str
    smiles: str
    nheavy: int
    started_at: str | None
    cell_id: str | None = None


class TrainingCompleted(BaseModel):
    mol_id: str
    smiles: str
    h298_pm7: float
    delta_kcal_mol: float | None
    elapsed_s: float | None
    conformers_evaluated: int | None
    finished_at: str


class TrainingActivity(BaseModel):
    """What the workers compute now and the latest accepted molecules."""

    in_flight: list[TrainingInFlight]
    recent: list[TrainingCompleted]
    mean_elapsed_s: float | None


class TrainingStatus(BaseModel):
    progress: TrainingProgress
    job: TrainingJob
    geometries: int
    failures: list[TrainingFailure]
    activity: TrainingActivity
    cpu_count: int


class TrainingStartRequest(BaseModel):
    workers: int = Field(3, ge=1, le=8)
    crest_threads: int = Field(4, ge=1, le=64)
    limit: int | None = Field(None, ge=1)
