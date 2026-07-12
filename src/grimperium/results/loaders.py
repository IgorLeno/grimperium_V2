"""CSV loaders and adapters for results analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from grimperium.calculation.contracts.quantity import KJ_PER_KCAL

LEGACY_REQUIRED_COLUMNS = {"H298_cbs", "H298_predicted"}
CANONICAL_MARKER_COLUMNS = {"run_id", "role", "canonical_value"}


def load_analysis_dataframe(dataset_path: Path | str) -> pd.DataFrame:
    """Load a legacy wide or canonical long-form results CSV."""
    path = Path(dataset_path)
    df = pd.read_csv(path)
    if CANONICAL_MARKER_COLUMNS.issubset(df.columns):
        return canonical_long_form_to_analysis_dataframe(df)
    return df


def canonical_long_form_to_analysis_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Adapt canonical long-form calculation rows to analysis columns."""
    rows: dict[str, dict[str, Any]] = {}
    for _, row in df.iterrows():
        mol_id = _molecule_id(row)
        record = rows.setdefault(
            mol_id,
            {
                "mol_id": mol_id,
                "smiles": _string_or_none(row.get("molecule_smiles")),
                "run_id": _string_or_none(row.get("run_id")),
            },
        )
        value = _value_kcal(row)
        role = str(row.get("role", "")).lower()
        hamiltonian = str(row.get("hamiltonian", "")).upper()
        if role == "baseline":
            if hamiltonian in {"", "PM7"} or "H298_pm7" not in record:
                record["H298_pm7"] = value
        elif role == "correction":
            record["delta_correction"] = value
        elif role == "final":
            record["H298_predicted"] = value
        elif role in {"reference", "cbs"}:
            record["H298_cbs"] = value

    return pd.DataFrame(rows.values())


def _molecule_id(row: pd.Series) -> str:
    name = _string_or_none(row.get("molecule_name"))
    if name:
        return name
    smiles = _string_or_none(row.get("molecule_smiles"))
    if smiles:
        return smiles
    estimate_id = _string_or_none(row.get("estimate_id"))
    if estimate_id:
        return estimate_id
    return "unknown"


def _value_kcal(row: pd.Series) -> float:
    value_kcal = row.get("value_kcal_mol")
    if pd.notna(value_kcal) and str(value_kcal) != "":
        return float(value_kcal)
    value = float(row["canonical_value"])
    unit = str(row.get("canonical_unit", "kcal/mol"))
    if unit == "kJ/mol":
        return value / KJ_PER_KCAL
    return value


def _string_or_none(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value)
    return text if text else None
