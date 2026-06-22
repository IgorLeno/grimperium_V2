import pytest

from grimperium.calculation.output.csv_writer import CANONICAL_CSV_COLUMNS
from tests.calculation.test_contract_serialization import _sample_result


def _xlsx_writer_modules():  # type: ignore[no-untyped-def]
    openpyxl = pytest.importorskip("openpyxl")
    from grimperium.calculation.output.xlsx_writer import (
        ARTIFACT_COLUMNS,
        STAGE_EXECUTION_COLUMNS,
        write_canonical_xlsx,
    )

    return (
        openpyxl.load_workbook,
        write_canonical_xlsx,
        STAGE_EXECUTION_COLUMNS,
        ARTIFACT_COLUMNS,
    )


def _headers(sheet) -> list[str]:  # type: ignore[no-untyped-def]
    return [cell.value for cell in sheet[1]]


def test_write_canonical_xlsx_creates_required_sheets(tmp_path) -> None:  # type: ignore[no-untyped-def]
    output_path = tmp_path / "results.xlsx"
    load_workbook, write_canonical_xlsx, _, _ = _xlsx_writer_modules()

    write_canonical_xlsx([_sample_result()], output_path)

    workbook = load_workbook(output_path)
    assert workbook.sheetnames == ["estimates", "stage_executions", "artifacts"]


def test_write_canonical_xlsx_estimates_headers(tmp_path) -> None:  # type: ignore[no-untyped-def]
    output_path = tmp_path / "results.xlsx"
    load_workbook, write_canonical_xlsx, _, _ = _xlsx_writer_modules()

    write_canonical_xlsx([_sample_result()], output_path)

    workbook = load_workbook(output_path)
    assert _headers(workbook["estimates"]) == CANONICAL_CSV_COLUMNS


def test_write_canonical_xlsx_stage_execution_headers(tmp_path) -> None:  # type: ignore[no-untyped-def]
    output_path = tmp_path / "results.xlsx"
    (
        load_workbook,
        write_canonical_xlsx,
        STAGE_EXECUTION_COLUMNS,
        _,
    ) = _xlsx_writer_modules()

    write_canonical_xlsx([_sample_result()], output_path)

    workbook = load_workbook(output_path)
    assert _headers(workbook["stage_executions"]) == STAGE_EXECUTION_COLUMNS


def test_write_canonical_xlsx_artifact_headers(tmp_path) -> None:  # type: ignore[no-untyped-def]
    output_path = tmp_path / "results.xlsx"
    load_workbook, write_canonical_xlsx, _, ARTIFACT_COLUMNS = _xlsx_writer_modules()

    write_canonical_xlsx([_sample_result()], output_path)

    workbook = load_workbook(output_path)
    assert _headers(workbook["artifacts"]) == ARTIFACT_COLUMNS
