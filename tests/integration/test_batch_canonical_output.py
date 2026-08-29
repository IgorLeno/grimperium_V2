"""Integration tests for canonical batch output wiring (Phase 3)."""

from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from grimperium.calculation.output.csv_writer import CANONICAL_CSV_COLUMNS
from grimperium.crest_pm7.batch.enums import BatchFailurePolicy, MoleculeStatus
from grimperium.crest_pm7.batch.execution_manager import BatchExecutionManager
from grimperium.crest_pm7.batch.models import Batch, BatchMolecule
from grimperium.crest_pm7.batch.output_contracts import BatchOutputLayout
from grimperium.crest_pm7.config import CRESTStatus, MOPACStatus
from grimperium.crest_pm7.csv_enhancements import CSVManagerExtensions
from grimperium.crest_pm7.molecule_processor import ConformerData, PM7Result

# Operational columns that must NEVER leak into the canonical scientific output.
_OPERATIONAL_COLUMNS = {"status", "assigned_worker", "worker_status", "reruns"}


def _success_result(mol_id: str, hof: float) -> PM7Result:
    conformer = ConformerData(index=0, mol_id=mol_id, crest_rank=1)
    conformer.crest_status = CRESTStatus.SUCCESS
    conformer.crest_geometry_file = Path(f"work/{mol_id}/conf_1.xyz")
    conformer.mopac_status = MOPACStatus.SUCCESS
    conformer.mopac_output_file = Path(f"work/{mol_id}/conf_1.out")
    conformer.energy_hof = hof
    conformer.hof_extraction_successful = True
    return PM7Result(
        mol_id=mol_id,
        smiles="CCO",
        phase="A",
        nheavy=3,
        rdkit_descriptors={},
        crest_status=CRESTStatus.SUCCESS,
        crest_conformers_generated=3,
        crest_time=5.0,
        conformers=[conformer],
        num_conformers_selected=1,
        k_selected_pm7=1,
        total_execution_time=20.0,
        success=True,
    )


def _failed_result(mol_id: str) -> PM7Result:
    return PM7Result(
        mol_id=mol_id,
        smiles="CCO",
        conformers=[],
        success=False,
        error_message="MOPAC failed",
    )


class _FakeProcessor:
    def __init__(self, results: dict[str, PM7Result]) -> None:
        self._results = results

    def update_timeouts(self, **_kwargs: object) -> None:
        pass

    def process_with_fixed_timeout(
        self, mol_id: str, smiles: str, progress_callback: object = None
    ) -> PM7Result:
        return self._results[mol_id]


def _batch(mol_ids: list[str]) -> Batch:
    return Batch(
        batch_id="batch_0001",
        molecules=[
            BatchMolecule(mol_id=m, smiles="CCO", batch_order=i, nheavy=3)
            for i, m in enumerate(mol_ids, start=1)
        ],
        crest_timeout_minutes=30,
        mopac_timeout_minutes=60,
        failure_policy=BatchFailurePolicy.PARTIAL_OK,
        method_id="pm7_delta_learning",
        method_version="1.0.0",
        method_snapshot={},
    )


def _make_manager(
    tmp_path: Path,
    results: dict[str, PM7Result],
    *,
    output_dir: Path | None = None,
    result_writer: object = None,
    write_xlsx: bool = True,
) -> tuple[BatchExecutionManager, MagicMock, MagicMock]:
    csv_manager = MagicMock()
    csv_manager.pm7result_to_csv_update.return_value = {"H298_pm7": -55.0}
    csv_manager.get_reference_hof.return_value = -50.0
    state_manager = MagicMock()
    state_manager.record_rerun_or_skip.return_value = MoleculeStatus.RERUN.value
    config = MagicMock()
    config.temp_dir = tmp_path
    layout = BatchOutputLayout(output_dir) if output_dir is not None else None
    extra = {"result_writer": result_writer} if result_writer is not None else {}
    manager = BatchExecutionManager(
        csv_manager=csv_manager,
        state_manager=state_manager,
        detail_manager=MagicMock(),
        pm7_config=config,
        processor_adapter=_FakeProcessor(results),
        output_layout=layout,
        write_xlsx=write_xlsx,
        **extra,  # type: ignore[arg-type]
    )
    return manager, csv_manager, state_manager


def _run(manager: BatchExecutionManager, mol_ids: list[str]) -> None:
    with patch.object(
        CSVManagerExtensions, "update_molecule_with_mopac_results", return_value=True
    ):
        manager.execute_batch(_batch(mol_ids))


def _read_canonical(output_dir: Path) -> list[dict[str, str]]:
    with (output_dir / "calculation_results.csv").open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_batch_writes_canonical_csv_for_one_success(tmp_path: Path) -> None:
    out = tmp_path / "out"
    manager, _, _ = _make_manager(
        tmp_path, {"m1": _success_result("m1", -55.0)}, output_dir=out
    )
    _run(manager, ["m1"])

    assert (out / "calculation_results.csv").exists()
    rows = _read_canonical(out)
    assert len(rows) == 1
    assert rows[0]["role"] == "baseline"
    assert rows[0]["method_id"]
    assert rows[0]["schema_version"] == "1"


def test_batch_writes_authoritative_run_id(tmp_path: Path) -> None:
    out = tmp_path / "out"
    manager, _, _ = _make_manager(
        tmp_path, {"m1": _success_result("m1", -55.0)}, output_dir=out
    )
    manager.canonical_run_id = "run_from_manifest"

    _run(manager, ["m1"])

    rows = _read_canonical(out)
    assert {row["run_id"] for row in rows} == {"run_from_manifest"}


def test_batch_generates_xlsx_when_enabled(tmp_path: Path) -> None:
    pytest.importorskip(
        "openpyxl", reason="XLSX integration requires the optional output extra"
    )
    out = tmp_path / "out"
    manager, _, _ = _make_manager(
        tmp_path, {"m1": _success_result("m1", -55.0)}, output_dir=out, write_xlsx=True
    )
    _run(manager, ["m1"])
    assert (out / "calculation_results.xlsx").exists()


def test_batch_skips_xlsx_when_disabled(tmp_path: Path) -> None:
    out = tmp_path / "out"
    manager, _, _ = _make_manager(
        tmp_path, {"m1": _success_result("m1", -55.0)}, output_dir=out, write_xlsx=False
    )
    _run(manager, ["m1"])
    assert (out / "calculation_results.csv").exists()
    assert not (out / "calculation_results.xlsx").exists()


def test_batch_partial_success_writes_only_successful(tmp_path: Path) -> None:
    out = tmp_path / "out"
    results = {
        "m1": _success_result("m1", -55.0),
        "m2": _failed_result("m2"),
        "m3": _success_result("m3", -60.0),
    }
    manager, _, _ = _make_manager(tmp_path, results, output_dir=out)
    _run(manager, ["m1", "m2", "m3"])

    rows = _read_canonical(out)
    names = {row["molecule_name"] for row in rows}
    assert names == {"m1", "m3"}  # failed m2 absent


def test_batch_total_failure_writes_empty_canonical_csv(tmp_path: Path) -> None:
    out = tmp_path / "out"
    manager, _, _ = _make_manager(
        tmp_path, {"m1": _failed_result("m1")}, output_dir=out
    )
    _run(manager, ["m1"])
    assert (out / "calculation_results.csv").exists()
    assert _read_canonical(out) == []  # header only, no data rows


def test_batch_pm7_only_never_writes_ml_final_or_correction(tmp_path: Path) -> None:
    out = tmp_path / "out"
    manager, _, _ = _make_manager(
        tmp_path, {"m1": _success_result("m1", -55.0)}, output_dir=out
    )
    _run(manager, ["m1"])
    roles = {row["role"] for row in _read_canonical(out)}
    assert roles == {"baseline"}
    assert "correction" not in roles
    assert "final" not in roles


def test_canonical_output_has_no_operational_columns(tmp_path: Path) -> None:
    out = tmp_path / "out"
    manager, _, _ = _make_manager(
        tmp_path, {"m1": _success_result("m1", -55.0)}, output_dir=out
    )
    _run(manager, ["m1"])
    with (out / "calculation_results.csv").open(encoding="utf-8") as handle:
        header = next(csv.reader(handle))
    assert header == CANONICAL_CSV_COLUMNS
    assert _OPERATIONAL_COLUMNS.isdisjoint(header)


def test_writer_called_with_actually_computed_results(tmp_path: Path) -> None:
    captured: dict[str, object] = {}

    def fake_writer(results, layout, *, include_xlsx=True):  # type: ignore[no-untyped-def]
        captured["results"] = list(results)
        captured["include_xlsx"] = include_xlsx
        report = MagicMock()
        report.result_count = len(results)
        report.csv_path = layout.calculation_results_csv
        return report

    out = tmp_path / "out"
    manager, _, _ = _make_manager(
        tmp_path,
        {"m1": _success_result("m1", -55.0), "m2": _failed_result("m2")},
        output_dir=out,
        result_writer=fake_writer,
    )
    _run(manager, ["m1", "m2"])

    results = captured["results"]
    assert len(results) == 1  # only the successful molecule
    assert results[0].estimates[0].role.value == "baseline"


def test_legacy_and_operational_paths_still_invoked(tmp_path: Path) -> None:
    out = tmp_path / "out"
    manager, csv_manager, state_manager = _make_manager(
        tmp_path, {"m1": _success_result("m1", -55.0)}, output_dir=out
    )
    _run(manager, ["m1"])

    # Legacy scientific CSV path preserved.
    csv_manager.mark_success.assert_called_once()
    # Operational state still marked OK via the shared result applier.
    state_manager.mark_success.assert_called_once_with("m1")


def test_no_canonical_files_when_output_disabled(tmp_path: Path) -> None:
    manager, _, _ = _make_manager(
        tmp_path, {"m1": _success_result("m1", -55.0)}, output_dir=None
    )
    _run(manager, ["m1"])
    assert not (tmp_path / "calculation_results.csv").exists()
    assert manager._canonical_results == []
