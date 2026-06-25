import csv

import pytest

from grimperium.calculation.contracts.quantity import KJ_PER_KCAL
from grimperium.calculation.output.csv_writer import CANONICAL_CSV_COLUMNS
from grimperium.crest_pm7.batch.output_contracts import (
    BATCH_STATE_COLUMNS,
    SCIENTIFIC_RESULT_COLUMNS,
    BatchOutputLayout,
    BatchResultWriteReport,
    write_batch_calculation_results,
)
from tests.calculation.test_contract_serialization import _sample_result


def _read_rows(path):  # type: ignore[no-untyped-def]
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_batch_state_schema_keeps_scientific_results_out() -> None:
    assert "mol_id" in BATCH_STATE_COLUMNS
    assert "status" in BATCH_STATE_COLUMNS
    assert "batch_id" in BATCH_STATE_COLUMNS
    assert "assigned_worker" in BATCH_STATE_COLUMNS
    assert "method_id" in BATCH_STATE_COLUMNS
    assert "method_version" in BATCH_STATE_COLUMNS
    assert "method_definition_snapshot" in BATCH_STATE_COLUMNS

    forbidden_columns = set(SCIENTIFIC_RESULT_COLUMNS)
    assert forbidden_columns
    assert not forbidden_columns.intersection(BATCH_STATE_COLUMNS)


def test_batch_output_layout_uses_split_default_filenames(tmp_path) -> None:  # type: ignore[no-untyped-def]
    layout = BatchOutputLayout(tmp_path)

    assert layout.batch_state_csv == tmp_path / "batch_state.csv"
    assert layout.calculation_results_csv == tmp_path / "calculation_results.csv"
    assert layout.calculation_results_xlsx == tmp_path / "calculation_results.xlsx"


def test_batch_calculation_result_writer_reuses_canonical_csv(tmp_path) -> None:  # type: ignore[no-untyped-def]
    report = write_batch_calculation_results(
        [_sample_result()],
        BatchOutputLayout(tmp_path),
        units="both",
        include_xlsx=False,
    )

    assert report == BatchResultWriteReport(
        csv_path=tmp_path / "calculation_results.csv",
        xlsx_path=None,
        result_count=1,
        estimate_row_count=1,
    )
    rows = _read_rows(report.csv_path)
    assert list(rows[0].keys()) == CANONICAL_CSV_COLUMNS
    assert rows[0]["estimate_id"] == "estimate-001"
    assert rows[0]["value_kcal_mol"] == "-56.7890123456789"
    assert float(rows[0]["value_kj_mol"]) == pytest.approx(
        -56.7890123456789 * KJ_PER_KCAL
    )


def test_batch_calculation_result_writer_can_emit_xlsx(tmp_path) -> None:  # type: ignore[no-untyped-def]
    openpyxl = pytest.importorskip("openpyxl")

    report = write_batch_calculation_results(
        [_sample_result()],
        BatchOutputLayout(tmp_path),
        include_xlsx=True,
    )

    assert report.xlsx_path == tmp_path / "calculation_results.xlsx"
    assert report.xlsx_path.exists()
    workbook = openpyxl.load_workbook(report.xlsx_path)
    assert workbook.sheetnames == ["estimates", "stage_executions", "artifacts"]
