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


def load_canonical_long_form(dataset_path: Path | str) -> pd.DataFrame | None:
    """Return the raw canonical long-form CSV, or None when not canonical."""
    path = Path(dataset_path)
    df = pd.read_csv(path)
    if CANONICAL_MARKER_COLUMNS.issubset(df.columns):
        return df
    return None


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


def join_reference_from_dataset(
    analysis_df: pd.DataFrame,
    dataset_path: Path | str,
) -> tuple[pd.DataFrame, list[str]]:
    """Join H298_cbs from an external dataset by mol_id, then smiles.

    Does not modify the source dataset. Returns the joined frame and warnings.
    """
    warnings: list[str] = []
    path = Path(dataset_path)
    if not path.exists():
        warnings.append(f"reference dataset not found: {path}")
        return analysis_df, warnings

    ref = pd.read_csv(path)
    if "H298_cbs" not in ref.columns:
        warnings.append(f"reference dataset missing H298_cbs: {path}")
        return analysis_df, warnings

    ref = ref.copy()
    if "mol_id" not in ref.columns and "molecule_name" in ref.columns:
        ref["mol_id"] = ref["molecule_name"]
    if "smiles" not in ref.columns and "molecule_smiles" in ref.columns:
        ref["smiles"] = ref["molecule_smiles"]

    joined = analysis_df.copy()
    if "H298_cbs" not in joined.columns:
        joined["H298_cbs"] = pd.NA

    by_mol: dict[str, list[float]] = {}
    by_smiles: dict[str, list[float]] = {}
    if "mol_id" in ref.columns:
        for _, row in ref.iterrows():
            mol_id = _string_or_none(row.get("mol_id"))
            value = row.get("H298_cbs")
            if mol_id is None or pd.isna(value):
                continue
            by_mol.setdefault(mol_id, []).append(float(value))
    if "smiles" in ref.columns:
        for _, row in ref.iterrows():
            smiles = _string_or_none(row.get("smiles"))
            value = row.get("H298_cbs")
            if smiles is None or pd.isna(value):
                continue
            by_smiles.setdefault(smiles, []).append(float(value))

    unmatched = 0
    ambiguous = 0
    for idx, row in joined.iterrows():
        if pd.notna(row.get("H298_cbs")):
            continue
        mol_id = _string_or_none(row.get("mol_id"))
        smiles = _string_or_none(row.get("smiles"))
        candidates: list[float] = []
        if mol_id is not None and mol_id in by_mol:
            candidates = by_mol[mol_id]
        elif smiles is not None and smiles in by_smiles:
            candidates = by_smiles[smiles]
        if not candidates:
            unmatched += 1
            continue
        unique_values = {round(v, 8) for v in candidates}
        if len(unique_values) > 1:
            ambiguous += 1
            warnings.append(
                f"ambiguous reference for molecule {mol_id or smiles}: "
                f"{sorted(unique_values)}"
            )
            continue
        joined.at[idx, "H298_cbs"] = float(candidates[0])

    if unmatched:
        warnings.append(f"{unmatched} molecule(s) could not be matched to reference")
    if ambiguous:
        warnings.append(f"{ambiguous} molecule(s) had ambiguous reference keys")
    joined["H298_cbs"] = pd.to_numeric(joined["H298_cbs"], errors="coerce")
    return joined, warnings


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
