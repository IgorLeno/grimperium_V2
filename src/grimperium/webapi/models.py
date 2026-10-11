"""Response bodies of the dashboard API (kcal/mol for every enthalpy)."""

from __future__ import annotations

from pydantic import BaseModel

LEGACY_PROTOCOL = "crest_pm7"
"""Protocol that produced ``thermo_pm7.csv``; the canonical one is
Semi-Imperium, so these rows are shown as legacy, never as training v2."""

BLANK_LABEL = "unset"
"""Count key for empty cells, so a blank flag is visible, not dropped."""


class Health(BaseModel):
    status: str
    read_only: bool


class LegacyDatasetProgress(BaseModel):
    protocol: str
    path: str
    total: int
    status: dict[str, int]
    cbs_quality_flag: dict[str, int]


class StoreProgress(BaseModel):
    path: str
    exists: bool
    molecules: int
    runs: int
    calculations: int
    states: dict[str, int]


class DatasetProgress(BaseModel):
    legacy: LegacyDatasetProgress
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
