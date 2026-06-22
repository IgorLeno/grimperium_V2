"""XLSX writer for canonical calculation results."""

from __future__ import annotations

from pathlib import Path

try:
    import openpyxl
except ImportError as e:
    raise ImportError(
        "openpyxl is required for XLSX export. "
        "Install with: pip install grimperium[output]"
    ) from e

from grimperium.calculation.contracts.models import MoleculeCalculationResult

from .csv_writer import CANONICAL_CSV_COLUMNS, UnitSelection, estimate_rows

STAGE_EXECUTION_COLUMNS = [
    "run_id",
    "molecule_smiles",
    "stage_id",
    "program",
    "role",
    "status",
    "requested",
    "execution_time_s",
    "error_message",
]

ARTIFACT_COLUMNS = [
    "run_id",
    "molecule_smiles",
    "artifact_id",
    "artifact_type",
    "hamiltonian",
    "conformer_index",
    "relative_path",
]


def _value_or_blank(value: object | None) -> object:
    return "" if value is None else value


def _write_sheet(
    workbook: openpyxl.Workbook,
    title: str,
    columns: list[str],
    rows: list[dict[str, object]],
    *,
    use_active: bool = False,
) -> None:
    worksheet = workbook.active if use_active else workbook.create_sheet(title)
    worksheet.title = title
    worksheet.append(columns)
    for row in rows:
        worksheet.append([row.get(column, "") for column in columns])


def _stage_rows(results: list[MoleculeCalculationResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        for stage in result.stage_executions:
            rows.append(
                {
                    "run_id": result.run.run_id,
                    "molecule_smiles": result.molecule.smiles,
                    "stage_id": stage.stage_id,
                    "program": stage.program,
                    "role": stage.role,
                    "status": stage.status.value,
                    "requested": stage.requested,
                    "execution_time_s": _value_or_blank(stage.execution_time_s),
                    "error_message": _value_or_blank(stage.error_message),
                }
            )
    return rows


def _artifact_rows(results: list[MoleculeCalculationResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        for artifact in result.artifacts:
            rows.append(
                {
                    "run_id": result.run.run_id,
                    "molecule_smiles": result.molecule.smiles,
                    "artifact_id": artifact.artifact_id,
                    "artifact_type": artifact.artifact_type.value,
                    "hamiltonian": _value_or_blank(artifact.hamiltonian),
                    "conformer_index": _value_or_blank(artifact.conformer_index),
                    "relative_path": artifact.relative_path,
                }
            )
    return rows


def write_canonical_xlsx(
    results: list[MoleculeCalculationResult],
    output_path: Path,
    units: UnitSelection = "kcal/mol",
) -> None:
    """Write canonical calculation results to a simple XLSX workbook."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    workbook = openpyxl.Workbook()
    _write_sheet(
        workbook,
        "estimates",
        CANONICAL_CSV_COLUMNS,
        estimate_rows(results, units=units),
        use_active=True,
    )
    _write_sheet(
        workbook,
        "stage_executions",
        STAGE_EXECUTION_COLUMNS,
        _stage_rows(results),
    )
    _write_sheet(workbook, "artifacts", ARTIFACT_COLUMNS, _artifact_rows(results))
    workbook.save(output_path)
