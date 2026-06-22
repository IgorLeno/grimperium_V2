"""CSV writer for canonical calculation results."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Literal

from grimperium.calculation.contracts.models import MoleculeCalculationResult
from grimperium.calculation.contracts.quantity import KJ_PER_KCAL

UnitSelection = Literal["kcal/mol", "kJ/mol", "both"]

CANONICAL_CSV_COLUMNS = [
    "estimate_id",
    "run_id",
    "molecule_smiles",
    "molecule_name",
    "method_id",
    "method_version",
    "property_id",
    "role",
    "hamiltonian",
    "canonical_value",
    "canonical_unit",
    "value_kcal_mol",
    "value_kj_mol",
    "conformer_source_id",
    "overall_status",
    "execution_phase",
    "started_at",
    "completed_at",
    "schema_version",
]


def _datetime_to_str(value) -> str:  # type: ignore[no-untyped-def]
    return value.isoformat() if value is not None else ""


def _value_or_blank(value: object | None) -> object:
    return "" if value is None else value


def estimate_rows(
    results: list[MoleculeCalculationResult],
    units: UnitSelection = "kcal/mol",
) -> list[dict[str, object]]:
    """Build canonical long-form estimate rows."""
    rows: list[dict[str, object]] = []
    for result in results:
        for estimate in result.estimates:
            canonical_kcal = estimate.value.to_kcal_mol()
            value_kcal_mol: float | str = ""
            value_kj_mol: float | str = ""
            if units in ("kcal/mol", "both"):
                value_kcal_mol = canonical_kcal
            if units in ("kJ/mol", "both"):
                value_kj_mol = canonical_kcal * KJ_PER_KCAL

            rows.append(
                {
                    "estimate_id": estimate.estimate_id,
                    "run_id": result.run.run_id,
                    "molecule_smiles": result.molecule.smiles,
                    "molecule_name": _value_or_blank(result.molecule.name),
                    "method_id": estimate.method_id,
                    "method_version": estimate.method_version,
                    "property_id": estimate.property_id,
                    "role": estimate.role.value,
                    "hamiltonian": _value_or_blank(estimate.hamiltonian),
                    "canonical_value": estimate.value.value,
                    "canonical_unit": estimate.value.unit,
                    "value_kcal_mol": value_kcal_mol,
                    "value_kj_mol": value_kj_mol,
                    "conformer_source_id": _value_or_blank(
                        estimate.conformer_source_id
                    ),
                    "overall_status": result.overall_status.value,
                    "execution_phase": result.run.execution_phase,
                    "started_at": _datetime_to_str(result.run.started_at),
                    "completed_at": _datetime_to_str(result.run.completed_at),
                    "schema_version": result.schema_version,
                }
            )
    return rows


def write_canonical_csv(
    results: list[MoleculeCalculationResult],
    output_path: Path,
    units: UnitSelection = "kcal/mol",
) -> None:
    """Write canonical long-form estimates CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = estimate_rows(results, units=units)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANONICAL_CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
