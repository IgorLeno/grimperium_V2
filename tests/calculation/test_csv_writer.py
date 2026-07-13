import csv
from datetime import datetime, timezone

import pytest

from grimperium.calculation.contracts.enums import OverallStatus, PropertyRole
from grimperium.calculation.contracts.models import (
    CalculationMethodReference,
    MoleculeCalculationResult,
    MoleculeData,
    PropertyEstimate,
    RunMetadata,
)
from grimperium.calculation.contracts.quantity import Quantity
from grimperium.calculation.output.csv_writer import write_canonical_csv


def _result_with_estimates(count: int = 1) -> MoleculeCalculationResult:
    estimates = [
        PropertyEstimate(
            estimate_id=f"estimate-{idx}",
            property_id="standard_enthalpy_of_formation",
            role=PropertyRole.BASELINE,
            method_id="pm7_delta_legacy",
            method_version="0.0.0",
            hamiltonian="PM7",
            value=Quantity(value=-10.0 - idx, unit="kcal/mol"),
            value_kcal_mol=None,
            value_kj_mol=None,
            conformer_source_id=idx + 1,
            uncertainty=None,
            model_path=None,
        )
        for idx in range(count)
    ]
    return MoleculeCalculationResult(
        molecule=MoleculeData(smiles="CCO", name="ethanol"),
        run=RunMetadata(
            run_id="run-001",
            execution_phase="A",
            method_ref=CalculationMethodReference(
                method_id="pm7_delta_legacy",
                method_version="0.0.0",
                property_id="standard_enthalpy_of_formation",
            ),
            started_at=datetime(2026, 6, 22, 10, 0, tzinfo=timezone.utc),
            completed_at=datetime(2026, 6, 22, 10, 1, tzinfo=timezone.utc),
            grimperium_version=None,
        ),
        overall_status=OverallStatus.SUCCESS,
        conformers=[],
        molecular_descriptors=None,
        estimates=estimates,
        artifacts=[],
        stage_executions=[],
    )


def _read_rows(path) -> list[dict[str, str]]:  # type: ignore[no-untyped-def]
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_write_canonical_csv_kcal_units(tmp_path) -> None:  # type: ignore[no-untyped-def]
    output_path = tmp_path / "results.csv"

    write_canonical_csv([_result_with_estimates()], output_path, units="kcal/mol")

    rows = _read_rows(output_path)
    assert rows[0]["canonical_value"] == "-10.0"
    assert rows[0]["canonical_unit"] == "kcal/mol"
    assert rows[0]["value_kcal_mol"] == "-10.0"
    assert rows[0]["value_kj_mol"] == ""


def test_write_canonical_csv_kj_units(tmp_path) -> None:  # type: ignore[no-untyped-def]
    output_path = tmp_path / "results.csv"

    write_canonical_csv([_result_with_estimates()], output_path, units="kJ/mol")

    rows = _read_rows(output_path)
    assert rows[0]["canonical_value"] == "-10.0"
    assert rows[0]["canonical_unit"] == "kcal/mol"
    assert rows[0]["value_kcal_mol"] == ""
    assert rows[0]["value_kj_mol"] == "-41.84"


def test_write_canonical_csv_both_units(tmp_path) -> None:  # type: ignore[no-untyped-def]
    output_path = tmp_path / "results.csv"

    write_canonical_csv([_result_with_estimates()], output_path, units="both")

    rows = _read_rows(output_path)
    assert rows[0]["canonical_value"] == "-10.0"
    assert rows[0]["canonical_unit"] == "kcal/mol"
    assert rows[0]["value_kcal_mol"] == "-10.0"
    assert rows[0]["value_kj_mol"] == "-41.84"


def test_write_canonical_csv_has_one_row_per_estimate(tmp_path) -> None:  # type: ignore[no-untyped-def]
    output_path = tmp_path / "results.csv"

    write_canonical_csv([_result_with_estimates(count=2)], output_path)

    rows = _read_rows(output_path)
    assert [row["estimate_id"] for row in rows] == ["estimate-0", "estimate-1"]


def test_write_canonical_csv_preserves_existing_file_on_failure(
    tmp_path,
    monkeypatch,
) -> None:  # type: ignore[no-untyped-def]
    output_path = tmp_path / "results.csv"
    original = "estimate_id,run_id\nold,run-001\n"
    output_path.write_text(original, encoding="utf-8")

    class FailingWriter:
        def __init__(self, handle, fieldnames):  # type: ignore[no-untyped-def]
            self.handle = handle

        def writeheader(self) -> None:
            self.handle.write("partial\n")

        def writerows(self, rows) -> None:  # type: ignore[no-untyped-def]
            raise RuntimeError("disk full")

    monkeypatch.setattr(
        "grimperium.calculation.output.csv_writer.csv.DictWriter",
        FailingWriter,
    )

    with pytest.raises(RuntimeError, match="disk full"):
        write_canonical_csv([_result_with_estimates()], output_path)

    assert output_path.read_text(encoding="utf-8") == original
