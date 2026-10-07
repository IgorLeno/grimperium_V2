"""Read the dashboard inputs; nothing here writes to disk.

Every failure to read an input (missing file, missing column, malformed
value, corrupt store record) becomes :class:`DataUnavailable`, which the app
turns into an explicit HTTP error instead of an empty chart.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from grimperium.webapi.models import (
    BLANK_LABEL,
    LEGACY_PROTOCOL,
    ComparisonPoint,
    DatasetProgress,
    ErrorStats,
    ExperimentalComparison,
    ExperimentalSummary,
    LegacyDatasetProgress,
    StoreProgress,
)

EXP_FLAGS = (
    "conflict_flag",
    "validation_only",
    "low_priority_only",
    "stewart_reference",
)
COMPARISON_UNITS = "kcal/mol; cbs_minus_exp = H298_cbs - H298_exp"


class DataUnavailable(Exception):
    """An input file cannot be read as the API expects."""

    def __init__(self, path: Path, reason: str) -> None:
        super().__init__(f"{path}: {reason}")
        self.path = path
        self.reason = reason


def _read_csv(path: Path, columns: Iterable[str]) -> pd.DataFrame:
    """Read ``columns`` as strings; blanks stay ``""``, never NaN."""
    wanted = list(columns)
    if not path.is_file():
        raise DataUnavailable(path, "file not found")
    try:
        header = pd.read_csv(path, nrows=0).columns
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError) as e:
        raise DataUnavailable(path, f"unreadable CSV ({e})") from e
    missing = [name for name in wanted if name not in header]
    if missing:
        raise DataUnavailable(path, f"missing columns {missing}")
    return pd.read_csv(path, usecols=wanted, dtype=str, keep_default_na=False)


def _counts(values: pd.Series) -> dict[str, int]:
    """Value counts, most frequent first, blanks under ``BLANK_LABEL``."""
    labelled = values.str.strip().replace("", BLANK_LABEL)
    return {str(k): int(v) for k, v in labelled.value_counts().items()}


def _bools(path: Path, frame: pd.DataFrame, column: str) -> pd.Series:
    values = frame[column].str.strip()
    bad = sorted(set(values) - {"True", "False"})
    if bad:
        raise DataUnavailable(path, f"column {column!r} is not boolean: {bad[:5]}")
    return values == "True"


def _floats(path: Path, frame: pd.DataFrame, column: str, *, blank: bool) -> pd.Series:
    """Numeric column; with ``blank=True`` empty cells become NaN, else fail."""
    text = frame[column].str.strip()
    values = pd.to_numeric(text.replace("", math.nan), errors="coerce")
    invalid = values.isna() & (text != "") if blank else values.isna()
    if invalid.any():
        rows = [int(i) for i in values.index[invalid][:5]]
        raise DataUnavailable(path, f"column {column!r} not numeric at rows {rows}")
    return values


# ---------------------------------------------------------------- progress


def legacy_progress(path: Path) -> LegacyDatasetProgress:
    frame = _read_csv(path, ("status", "cbs_quality_flag"))
    return LegacyDatasetProgress(
        protocol=LEGACY_PROTOCOL,
        path=str(path),
        total=len(frame),
        status=_counts(frame["status"]),
        cbs_quality_flag=_counts(frame["cbs_quality_flag"]),
    )


def store_progress(root: Path) -> StoreProgress:
    """Summarise the Semi-Imperium store; an absent store is empty, not an error.

    The store is created by its first run, so a missing directory is a
    legitimate "nothing computed yet" state.
    """
    # Local imports keep the semi_imperium package (and its domain model)
    # out of the import path of the rest of grimperium.webapi.
    from semi_imperium.persistence import SemiImperiumStore, StoreIntegrityError
    from semi_imperium.workflows.database import build_summary

    try:
        summary = build_summary(SemiImperiumStore(root))
    except (StoreIntegrityError, ValueError, KeyError) as e:
        raise DataUnavailable(root, f"unreadable store record ({e})") from e
    return StoreProgress(
        path=str(root),
        exists=root.is_dir(),
        molecules=summary.molecule_count,
        runs=summary.run_count,
        calculations=summary.calculation_count,
        states=summary.state_counts(),
    )


def dataset_progress(pm7_csv: Path, store_root: Path) -> DatasetProgress:
    return DatasetProgress(
        legacy=legacy_progress(pm7_csv), semi_imperium=store_progress(store_root)
    )


# ------------------------------------------------------------ experimental


def experimental_summary(path: Path) -> ExperimentalSummary:
    frame = _read_csv(path, ("source", "uncertainty_kind", *EXP_FLAGS))
    return ExperimentalSummary(
        path=str(path),
        rows=len(frame),
        by_source=_counts(frame["source"]),
        by_source_family=_counts(frame["source"].str.split(":").str[0]),
        flags={flag: int(_bools(path, frame, flag).sum()) for flag in EXP_FLAGS},
        uncertainty_kind=_counts(frame["uncertainty_kind"]),
    )


COMPARISON_COLUMNS = tuple(ComparisonPoint.model_fields)


def _stats(errors: pd.Series) -> ErrorStats:
    if errors.empty:
        return ErrorStats(n=0, mae=None, bias=None, rmse=None)
    return ErrorStats(
        n=len(errors),
        mae=float(errors.abs().mean()),
        bias=float(errors.mean()),
        rmse=float(math.sqrt((errors**2).mean())),
    )


def experimental_comparison(path: Path) -> ExperimentalComparison:
    frame = _read_csv(path, COMPARISON_COLUMNS)
    numbers = {
        name: _floats(path, frame, name, blank=name == "uncertainty")
        for name in ("nheavy", "H298_exp", "uncertainty", "H298_cbs", "cbs_minus_exp")
    }
    flags = {
        name: _bools(path, frame, name)
        for name in ("conflict_flag", "low_priority_only", "validation_only")
    }
    points = [
        ComparisonPoint(
            exp_id=frame.at[i, "exp_id"],
            mol_id=frame.at[i, "mol_id"],
            smiles=frame.at[i, "smiles"],
            nheavy=int(numbers["nheavy"][i]),
            source=frame.at[i, "source"],
            H298_exp=float(numbers["H298_exp"][i]),
            uncertainty=(
                None
                if math.isnan(numbers["uncertainty"][i])
                else float(numbers["uncertainty"][i])
            ),
            H298_cbs=float(numbers["H298_cbs"][i]),
            cbs_minus_exp=float(numbers["cbs_minus_exp"][i]),
            conflict_flag=bool(flags["conflict_flag"][i]),
            low_priority_only=bool(flags["low_priority_only"][i]),
            validation_only=bool(flags["validation_only"][i]),
        )
        for i in frame.index
    ]
    return ExperimentalComparison(
        path=str(path),
        units=COMPARISON_UNITS,
        stats=_stats(numbers["cbs_minus_exp"]),
        points=points,
    )
