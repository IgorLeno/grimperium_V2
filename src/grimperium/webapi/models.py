"""Response bodies of the dashboard API (kcal/mol for every enthalpy)."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, FiniteFloat, model_validator

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


# --------------------------------------------------------------------- cells
# Bodies of the cell API (``grimperium.webapi.cells``), called by the cells.


class CellRegisterRequest(BaseModel):
    host: str = Field(min_length=1, max_length=255)
    cores: int = Field(ge=1, le=1024)
    workers: int = Field(ge=1, le=256)
    crest_threads: int = Field(ge=1, le=256)
    #: ``probe_programs`` output on the cell (version, sha256, path per program).
    program_versions: dict[str, str]
    protocol_signature: str


class CellRegisterResponse(BaseModel):
    cell_id: str
    name: str
    admission: str
    #: Molecules the cell held before this registration, back in the queue.
    released: list[str]
    lease_s: float
    #: How often the cell should send a heartbeat (well under 90 s).
    heartbeat_s: float


class CellClaimRequest(BaseModel):
    count: int = Field(ge=1, le=256)


class CellMolecule(BaseModel):
    attempt_id: str
    mol_id: str
    smiles: str
    charge: int
    multiplicity: int
    nheavy: int


class CellClaimResponse(BaseModel):
    molecules: list[CellMolecule]
    lease_s: float


class CellHeartbeatRequest(BaseModel):
    attempt_ids: list[str] = Field(default_factory=list, max_length=256)


class CellHeartbeatResponse(BaseModel):
    renewed: int
    #: Attempts the cell reported that it no longer holds (released or settled).
    not_held: list[str]


class CellGeometry(BaseModel):
    """A protocol-2 result, as ``WorkerResult`` holds it on the cell."""

    mol_id: str
    smiles: str
    signature: str
    calculation_id: str
    conformer_index: int | None
    elements: list[str] = Field(min_length=1)
    coordinates: list[tuple[FiniteFloat, FiniteFloat, FiniteFloat]]
    hof_kcal_mol: FiniteFloat
    lowest_frequency_cm1: FiniteFloat | None
    conformers_evaluated: int = Field(ge=0)
    elapsed_s: FiniteFloat = Field(ge=0)

    @model_validator(mode="after")
    def _one_position_per_atom(self) -> CellGeometry:
        if len(self.elements) != len(self.coordinates):
            raise ValueError(
                f"{len(self.elements)} elementos e {len(self.coordinates)} posições"
            )
        return self


class CellResultRequest(BaseModel):
    """Exactly one of ``result`` (success) or ``failure`` (reason)."""

    attempt_id: str
    result: CellGeometry | None = None
    failure: str | None = Field(None, min_length=1, max_length=2000)

    @model_validator(mode="after")
    def _exactly_one(self) -> CellResultRequest:
        if (self.result is None) == (self.failure is None):
            raise ValueError("envie result ou failure, um dos dois")
        return self


class CellResultResponse(BaseModel):
    outcome: str


# Dashboard side of the cells (``grimperium.webapi.app``, local only).


class CellSummary(BaseModel):
    cell_id: str
    name: str
    #: online · late · silent · never_seen · version_mismatch · blocked · revoked
    state: str
    admission: str
    #: Admission test progress while pending, or why it blocked the cell.
    admission_note: str | None
    host: str | None
    cores: int | None
    workers: int | None
    crest_threads: int | None
    program_versions: dict[str, str]
    refusal: str | None
    created_at: str
    registered_at: str | None
    last_seen_at: str | None
    revoked_at: str | None
    #: Molecules reserved to this cell right now.
    held: int


class CellNameRequest(BaseModel):
    name: str = Field(min_length=1, max_length=80, pattern=r"\S")


class CellCreated(BaseModel):
    cell: CellSummary
    #: Shown only in this response; the central keeps just its hash.
    token: str


class CellReleased(BaseModel):
    released: list[str]


class CellPackageRequest(BaseModel):
    #: What the cell calls: the cell API on the central's Tailscale address.
    central_url: str = Field(
        min_length=8, max_length=200, pattern=r"^https?://[^\s/]+(:\d+)?/?$"
    )


class CellPackageDefaults(BaseModel):
    #: Suggested ``central_url`` (``""`` when the central's address is unknown).
    central_url: str


class CellAdmissionRequest(BaseModel):
    admission: Literal["pending", "passed", "blocked"]
