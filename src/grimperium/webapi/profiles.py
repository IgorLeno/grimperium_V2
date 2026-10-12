"""Descriptive profiles and paginated rows of the three reference tables.

A profile is what the dashboard needs to *explain* a table to a reader:
how many molecules, how large they are, which elements they contain, how
the enthalpies are spread and, for PM7, how far the cheap method lands from
the CBS reference. Nothing here writes to disk.

Unlike the strict readers in :mod:`grimperium.webapi.readers`, optional
columns that are absent are reported in ``missing_columns`` instead of
failing, because older PM7 exports do not carry every column.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from grimperium.webapi.models import (
    BLANK_LABEL,
    DatasetProfile,
    Histogram,
    HistogramBin,
    ProfilePair,
    SizeStats,
    TablePage,
)
from grimperium.webapi.readers import DataUnavailable

#: Upper bound of rows per page and of points in a parity cloud.
MAX_PAGE = 200
MAX_PAIRS = 4000
TARGET_BINS = 32

_HAS_N = re.compile(r"[Nn]")
_HAS_O = re.compile(r"[Oo]")


@dataclass(frozen=True)
class HistogramSpec:
    column: str
    label: str
    unit: str


@dataclass(frozen=True)
class TableSpec:
    key: str
    title: str
    id_column: str
    required: tuple[str, ...]
    optional: tuple[str, ...]
    table_columns: tuple[str, ...]
    search_columns: tuple[str, ...]
    numeric: tuple[str, ...]
    categories: dict[str, Callable[[pd.DataFrame], pd.Series]] = field(
        default_factory=dict
    )
    histograms: tuple[HistogramSpec, ...] = ()
    derived: Callable[[pd.DataFrame], pd.DataFrame] | None = None
    error_column: str | None = None
    """Signed error summarised per heavy-atom count, if any."""
    pair_columns: tuple[str, str] | None = None
    """(x, y) of a parity cloud, if any."""


def _column(name: str) -> Callable[[pd.DataFrame], pd.Series]:
    return lambda frame: frame[name]


def _source_family(frame: pd.DataFrame) -> pd.Series:
    return frame["source"].astype(str).str.split(":").str[0]


def _pm7_delta(frame: pd.DataFrame) -> pd.DataFrame:
    # Delta-learning target, same sign convention as the training code.
    frame["delta_cbs_pm7"] = frame["H298_cbs"] - frame["H298_pm7"]
    return frame


EXP_SPEC = TableSpec(
    key="exp",
    title="Experimental ΔHf",
    id_column="exp_id",
    required=("exp_id", "smiles", "nheavy", "H298_exp", "source"),
    optional=(
        "name",
        "formula",
        "phase",
        "uncertainty",
        "uncertainty_kind",
        "multiplicity",
        "charge",
        "n_sources",
        "validation_only",
        "conflict_flag",
        "low_priority_only",
    ),
    table_columns=(
        "exp_id",
        "name",
        "formula",
        "smiles",
        "nheavy",
        "H298_exp",
        "uncertainty",
        "source",
        "n_sources",
        "validation_only",
        "conflict_flag",
    ),
    search_columns=("exp_id", "name", "formula", "smiles", "source"),
    numeric=("nheavy", "H298_exp", "uncertainty", "n_sources"),
    categories={
        "source_family": _source_family,
        "phase": _column("phase"),
        "uncertainty_kind": _column("uncertainty_kind"),
        "multiplicity": _column("multiplicity"),
    },
    histograms=(
        HistogramSpec("H298_exp", "Entalpia de formação experimental", "kcal/mol"),
        HistogramSpec("uncertainty", "Incerteza informada", "kcal/mol"),
    ),
)

CBS_SPEC = TableSpec(
    key="cbs",
    title="Referência CBS (CHON)",
    id_column="mol_id",
    required=("mol_id", "smiles", "nheavy", "H298_cbs"),
    optional=(
        "H298_b3",
        "cbs_b3",
        "charge",
        "multiplicity",
        "smiles_count",
        "conformer_rank_by_h298",
    ),
    table_columns=(
        "mol_id",
        "smiles",
        "nheavy",
        "charge",
        "multiplicity",
        "H298_cbs",
        "H298_b3",
        "cbs_b3",
        "smiles_count",
    ),
    search_columns=("mol_id", "smiles"),
    numeric=(
        "nheavy",
        "H298_cbs",
        "H298_b3",
        "cbs_b3",
        "smiles_count",
        "conformer_rank_by_h298",
    ),
    categories={
        "multiplicity": _column("multiplicity"),
        "charge": _column("charge"),
        "smiles_count": _column("smiles_count"),
    },
    histograms=(
        HistogramSpec("H298_cbs", "Entalpia de formação CBS", "kcal/mol"),
        HistogramSpec("cbs_b3", "Correção CBS − B3LYP", "kcal/mol"),
    ),
)

PM7_SPEC = TableSpec(
    key="pm7",
    title="CREST + PM7 (legado)",
    id_column="mol_id",
    required=("mol_id", "status", "smiles"),
    optional=(
        "nheavy",
        "H298_cbs",
        "H298_pm7",
        "cbs_quality_flag",
        "crest_status",
        "mopac_status",
        "crest_conformers_generated",
        "num_conformers_selected",
        "crest_time_s",
        "mopac_time_s",
        "reruns",
    ),
    table_columns=(
        "mol_id",
        "smiles",
        "status",
        "nheavy",
        "H298_cbs",
        "H298_pm7",
        "delta_cbs_pm7",
        "cbs_quality_flag",
        "crest_status",
        "mopac_status",
        "crest_conformers_generated",
        "crest_time_s",
        "mopac_time_s",
    ),
    search_columns=("mol_id", "smiles", "status"),
    numeric=(
        "nheavy",
        "H298_cbs",
        "H298_pm7",
        "crest_conformers_generated",
        "num_conformers_selected",
        "crest_time_s",
        "mopac_time_s",
        "reruns",
    ),
    categories={
        "status": _column("status"),
        "cbs_quality_flag": _column("cbs_quality_flag"),
        "crest_status": _column("crest_status"),
        "mopac_status": _column("mopac_status"),
    },
    histograms=(
        HistogramSpec("delta_cbs_pm7", "Erro do PM7 (CBS − PM7)", "kcal/mol"),
        HistogramSpec("H298_pm7", "Entalpia de formação PM7", "kcal/mol"),
        HistogramSpec("crest_time_s", "Tempo do CREST", "s"),
        HistogramSpec("mopac_time_s", "Tempo do MOPAC", "s"),
    ),
    derived=_pm7_delta,
    error_column="delta_cbs_pm7",
    pair_columns=("H298_cbs", "H298_pm7"),
)

SPECS = {spec.key: spec for spec in (EXP_SPEC, CBS_SPEC, PM7_SPEC)}

_cache: dict[tuple[str, str], tuple[int, pd.DataFrame, list[str]]] = {}


def _load(path: Path, spec: TableSpec) -> tuple[pd.DataFrame, list[str]]:
    """Read the spec's columns once per file version (mtime)."""
    if not path.is_file():
        raise DataUnavailable(path, "file not found")
    mtime = path.stat().st_mtime_ns
    cached = _cache.get((str(path), spec.key))
    if cached and cached[0] == mtime:
        return cached[1], cached[2]
    try:
        header = list(pd.read_csv(path, nrows=0).columns)
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as e:
        raise DataUnavailable(path, f"unreadable CSV ({e})") from e
    absent = [name for name in spec.required if name not in header]
    if absent:
        raise DataUnavailable(path, f"missing columns {absent}")
    present = [name for name in (*spec.required, *spec.optional) if name in header]
    missing = [name for name in spec.optional if name not in header]
    frame = pd.read_csv(path, usecols=present, dtype=str, keep_default_na=False)
    for name in spec.numeric:
        if name in frame:
            # Blank and malformed cells both become NaN (descriptive, not strict).
            frame[name] = pd.to_numeric(frame[name].str.strip(), errors="coerce")
    if spec.derived is not None and {"H298_cbs", "H298_pm7"} <= set(frame):
        frame = spec.derived(frame)
    _cache[(str(path), spec.key)] = (mtime, frame, missing)
    return frame, missing


def _nice_step(span: float) -> float:
    """Round a raw bin width to 1, 2 or 5 times a power of ten."""
    if span <= 0 or not math.isfinite(span):
        return 1.0
    raw = span / TARGET_BINS
    power = 10 ** math.floor(math.log10(raw))
    for factor in (1, 2, 5, 10):
        if raw <= factor * power:
            return float(factor * power)
    return float(10 * power)


def _optional(value: float) -> float | None:
    return float(value) if math.isfinite(value) else None


def histogram(values: pd.Series, spec: HistogramSpec) -> Histogram:
    data = values.dropna().to_numpy(dtype=float)
    if data.size == 0:
        return Histogram(
            column=spec.column,
            label=spec.label,
            unit=spec.unit,
            n=0,
            below=0,
            above=0,
            mean=None,
            median=None,
            p05=None,
            p95=None,
            minimum=None,
            maximum=None,
            bins=[],
        )
    lo_q, hi_q = np.quantile(data, [0.005, 0.995])
    step = _nice_step(float(hi_q - lo_q))
    lo = math.floor(lo_q / step) * step
    hi = math.ceil(hi_q / step) * step
    if hi <= lo:
        hi = lo + step
    edges = np.arange(lo, hi + step / 2, step)
    inside = data[(data >= lo) & (data <= hi)]
    counts, _ = np.histogram(inside, bins=edges)
    p05, median, p95 = np.quantile(data, [0.05, 0.5, 0.95])
    return Histogram(
        column=spec.column,
        label=spec.label,
        unit=spec.unit,
        n=int(data.size),
        below=int((data < lo).sum()),
        above=int((data > hi).sum()),
        mean=_optional(float(data.mean())),
        median=_optional(float(median)),
        p05=_optional(float(p05)),
        p95=_optional(float(p95)),
        minimum=_optional(float(data.min())),
        maximum=_optional(float(data.max())),
        bins=[
            HistogramBin(
                start=round(float(edges[i]), 6),
                end=round(float(edges[i + 1]), 6),
                count=int(count),
            )
            for i, count in enumerate(counts)
        ],
    )


def _labelled_counts(values: pd.Series) -> dict[str, int]:
    text = values.astype(str).str.strip()
    # Numeric categories stored as floats ("1.0") read better as "1".
    text = text.str.replace(r"^(-?\d+)\.0$", r"\1", regex=True)
    text = text.replace({"": BLANK_LABEL, "nan": BLANK_LABEL})
    return {str(k): int(v) for k, v in text.value_counts().items()}


def composition(smiles: pd.Series) -> dict[str, int]:
    """Element family of each CHON SMILES: which of N and O it contains."""
    text = smiles.astype(str)
    has_n = text.str.contains(_HAS_N)
    has_o = text.str.contains(_HAS_O)
    return {
        "CH": int((~has_n & ~has_o).sum()),
        "CHO": int((~has_n & has_o).sum()),
        "CHN": int((has_n & ~has_o).sum()),
        "CHNO": int((has_n & has_o).sum()),
    }


def _manifest(path: Path) -> dict[str, object] | None:
    manifest = path.with_suffix(".manifest.json")
    if not manifest.is_file():
        manifest = path.with_name(path.stem + "_v2.manifest.json")
    if not manifest.is_file():
        return None
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    keep = ("created_at", "source_citation", "filter", "rows", "unique_smiles")
    return {key: payload[key] for key in keep if key in payload}


def profile(path: Path, key: str) -> DatasetProfile:
    spec = SPECS[key]
    frame, missing = _load(path, spec)
    nheavy: dict[str, int] = {}
    if "nheavy" in frame:
        sizes = frame["nheavy"].dropna().astype(int).value_counts().sort_index()
        nheavy = {str(k): int(v) for k, v in sizes.items()}
    categories = {
        name: _labelled_counts(read(frame))
        for name, read in spec.categories.items()
        if _has_inputs(frame, name)
    }
    histograms = [
        histogram(frame[item.column], item)
        for item in spec.histograms
        if item.column in frame
    ]
    by_nheavy: list[SizeStats] = []
    if spec.error_column and spec.error_column in frame and "nheavy" in frame:
        errors = frame[[spec.error_column, "nheavy"]].dropna()
        for size, group in errors.groupby(errors["nheavy"].astype(int)):
            values = group[spec.error_column]
            by_nheavy.append(
                SizeStats(
                    nheavy=int(size),
                    n=len(values),
                    mean=float(values.mean()),
                    mean_abs=float(values.abs().mean()),
                )
            )
    pairs: list[ProfilePair] = []
    pairs_total = 0
    if spec.pair_columns and set(spec.pair_columns) <= set(frame):
        x, y = spec.pair_columns
        both = frame.dropna(subset=[x, y])
        pairs_total = len(both)
        if pairs_total > MAX_PAIRS:
            # Deterministic thinning: every k-th row keeps the shape of the cloud.
            both = both.iloc[:: math.ceil(pairs_total / MAX_PAIRS)]
        pairs = [
            ProfilePair(
                id=str(row[spec.id_column]),
                x=float(row[x]),
                y=float(row[y]),
                nheavy=(
                    int(row["nheavy"])
                    if "nheavy" in both and pd.notna(row["nheavy"])
                    else None
                ),
            )
            for _, row in both.iterrows()
        ]
    return DatasetProfile(
        key=spec.key,
        title=spec.title,
        path=str(path),
        rows=len(frame),
        unique_smiles=int(frame["smiles"].nunique()),
        nheavy=nheavy,
        composition=composition(frame["smiles"]),
        categories=categories,
        histograms=histograms,
        by_nheavy=by_nheavy,
        pairs=pairs,
        pairs_total=pairs_total,
        missing_columns=missing,
        manifest=_manifest(path),
    )


def _has_inputs(frame: pd.DataFrame, category: str) -> bool:
    if category == "source_family":
        return "source" in frame
    return category in frame


def _cell(value: object) -> str | float | int | bool | None:
    if value is None:
        return None
    if isinstance(value, bool | np.bool_):
        return bool(value)
    if isinstance(value, int | np.integer):
        return int(value)
    if isinstance(value, float | np.floating):
        return float(value) if math.isfinite(float(value)) else None
    text = str(value)
    if text in ("True", "False"):
        return text == "True"
    return text if text != "" else None


def rows(
    path: Path,
    key: str,
    *,
    query: str = "",
    offset: int = 0,
    limit: int = 50,
    sort: str | None = None,
    descending: bool = False,
) -> TablePage:
    spec = SPECS[key]
    frame, _missing = _load(path, spec)
    columns = [name for name in spec.table_columns if name in frame]
    view = frame
    term = query.strip().lower()
    if term:
        mask = pd.Series(False, index=frame.index)
        for name in spec.search_columns:
            if name in frame:
                mask |= (
                    frame[name].astype(str).str.lower().str.contains(term, regex=False)
                )
        view = frame[mask]
    if sort and sort in columns:
        view = view.sort_values(
            sort, ascending=not descending, na_position="last", kind="stable"
        )
    limit = max(1, min(limit, MAX_PAGE))
    offset = max(0, offset)
    page = view.iloc[offset : offset + limit]
    return TablePage(
        key=key,
        columns=columns,
        total=len(frame),
        matched=len(view),
        offset=offset,
        limit=limit,
        rows=[
            {name: _cell(record[name]) for name in columns}
            for record in page.to_dict("records")
        ],
    )
