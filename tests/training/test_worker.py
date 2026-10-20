from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from semi_imperium.conformers import ConformerGeometry
from semi_imperium.domain import (
    CalculationResultData,
    CalculationState,
    VerificationOutcome,
)
from semi_imperium.training.worker import WorkerTask, run_task
from semi_imperium.workflows.calculation import ExecutionOutcome, ExecutionRequest

GEOMETRY = ConformerGeometry(
    elements=("C", "H", "H", "H", "H"),
    coordinates=(
        (0.0, 0.0, 0.0),
        (0.63, 0.63, 0.63),
        (-0.63, -0.63, 0.63),
        (-0.63, 0.63, -0.63),
        (0.63, -0.63, -0.63),
    ),
)


def _with_outputs(
    store: Path, request: ExecutionRequest, outcome: ExecutionOutcome
) -> ExecutionOutcome:
    """Leave the two MOPAC outputs of a verified result where MOPAC would."""
    if outcome.state is not CalculationState.VERIFIED or outcome.result is None:
        return outcome
    folder = store / "artifacts" / request.calculation_id / "mopac"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "opt.out").write_text("OPT")
    (folder / "freq.out").write_text("FREQ")
    base = f"artifacts/{request.calculation_id}/mopac"
    paths = (f"{base}/opt.out", f"{base}/freq.out")
    return replace(outcome, result=replace(outcome.result, artifact_paths=paths))


class FakeExecutor:
    def __init__(self, outcome: ExecutionOutcome, *, store: Path | None = None) -> None:
        self.outcome = outcome
        self.store = store
        self.requests: list[ExecutionRequest] = []

    def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
        self.requests.append(request)
        if self.store is None:
            return self.outcome
        return _with_outputs(self.store, request, self.outcome)


def _task(tmp_path: Path, **overrides: object) -> WorkerTask:
    values: dict[str, object] = {
        "mol_id": "mol_00001",
        "smiles": "C",
        "charge": 0,
        "multiplicity": 1,
        "run_id": "run-test",
        "store_root": str(tmp_path / "store"),
        "work_dir": str(tmp_path / "work"),
        "crest_executable": "crest",
        "xtb_executable": "xtb",
        "mopac_executable": "mopac",
        "crest_threads": 1,
        "crest_timeout_s": 60.0,
        "mopac_timeout_s": 60.0,
        "keep_artifacts": False,
    }
    values.update(overrides)
    return WorkerTask(**values)  # type: ignore[arg-type]


def _verified() -> ExecutionOutcome:
    return ExecutionOutcome(
        state=CalculationState.VERIFIED,
        verification=VerificationOutcome.MINIMUM_CONFIRMED,
        result=CalculationResultData(
            energy_hof_kcal_mol=-17.9,
            conformer_index=2,
            conformers_evaluated=3,
            lowest_frequency_cm1=1306.0,
        ),
        geometry=GEOMETRY,
    )


def test_verified_minimum_is_ok_and_carries_the_geometry(tmp_path: Path) -> None:
    executor = FakeExecutor(_verified(), store=tmp_path / "store")
    result = run_task(_task(tmp_path), executor=executor)
    assert result.ok is True and result.state == "verified"
    assert result.hof_kcal_mol == -17.9 and result.conformers_evaluated == 3
    assert result.elements == GEOMETRY.elements
    assert result.coordinates == GEOMETRY.coordinates
    request = executor.requests[0]
    assert request.hamiltonian == "PM7"
    assert request.configuration.conformer_selection.folding_filter.enabled is True
    assert result.signature == request.configuration.signature().digest


def test_charge_and_multiplicity_reach_the_identity(tmp_path: Path) -> None:
    executor = FakeExecutor(_verified(), store=tmp_path / "store")
    run_task(_task(tmp_path, smiles="[CH3]", multiplicity=2), executor=executor)
    identity = executor.requests[0].identity
    assert (identity.charge, identity.multiplicity) == (0, 2)


def test_a_saddle_is_never_accepted(tmp_path: Path) -> None:
    outcome = ExecutionOutcome(
        state=CalculationState.SADDLE,
        verification=VerificationOutcome.SADDLE_POINT,
        result=CalculationResultData(
            energy_hof_kcal_mol=-10.0, lowest_frequency_cm1=-420.0
        ),
        geometry=GEOMETRY,
    )
    result = run_task(_task(tmp_path), executor=FakeExecutor(outcome))
    assert result.ok is False and result.state == "saddle"
    assert "saddle" in (result.error or "")


def test_executor_exception_becomes_a_failed_result(tmp_path: Path) -> None:
    class Boom:
        def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
            raise RuntimeError("crest exploded")

    result = run_task(_task(tmp_path), executor=Boom())
    assert result.ok is False and result.state == "failed"
    assert "crest exploded" in (result.error or "")


def test_success_removes_artifacts_but_failure_keeps_them(tmp_path: Path) -> None:
    store = tmp_path / "store"

    class Leaves:
        def __init__(self, outcome: ExecutionOutcome) -> None:
            self.outcome = outcome

        def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
            folder = store / "artifacts" / request.calculation_id
            folder.mkdir(parents=True)
            (folder / "trace.out").write_text("x")
            return _with_outputs(store, request, self.outcome)

    ok = run_task(_task(tmp_path), executor=Leaves(_verified()))
    assert not (store / "artifacts" / ok.calculation_id).exists()

    failed_outcome = ExecutionOutcome(
        state=CalculationState.FAILED,
        verification=VerificationOutcome.FAILED,
        error_message="no geometry",
    )
    bad = run_task(_task(tmp_path, run_id="run-2"), executor=Leaves(failed_outcome))
    assert (store / "artifacts" / bad.calculation_id / "trace.out").exists()


def test_the_calculation_record_is_persisted_in_the_store(tmp_path: Path) -> None:
    from semi_imperium.persistence import SemiImperiumStore

    result = run_task(
        _task(tmp_path),
        executor=FakeExecutor(_verified(), store=tmp_path / "store"),
    )
    records = SemiImperiumStore(tmp_path / "store").list_calculations()
    assert [r.calculation_id for r in records] == [result.calculation_id]
    assert records[0].state is CalculationState.VERIFIED
    assert records[0].result is not None and records[0].result.artifact_paths == ()


def test_success_keeps_the_chosen_conformer_outputs_outside_the_workdir(
    tmp_path: Path,
) -> None:
    store = tmp_path / "store"
    result = run_task(_task(tmp_path), executor=FakeExecutor(_verified(), store=store))
    folder = store / "geometrias" / "mol_00001"
    assert sorted(p.name for p in folder.iterdir()) == [
        "mol_00001.freq.out",
        "mol_00001.json",
        "mol_00001.out",
        "mol_00001.xyz",
    ]
    data = json.loads((folder / "mol_00001.json").read_text())
    assert (data["conformer_index"], data["hof_kcal_mol"]) == (2, -17.9)
    assert data["lowest_frequency_cm1"] == 1306.0
    assert data["calculation_id"] == result.calculation_id
    assert not (store / "artifacts" / result.calculation_id).exists()


def test_a_saddle_or_failure_keeps_no_final_outputs(tmp_path: Path) -> None:
    outcome = ExecutionOutcome(
        state=CalculationState.SADDLE,
        verification=VerificationOutcome.SADDLE_POINT,
        result=CalculationResultData(energy_hof_kcal_mol=-10.0),
        geometry=GEOMETRY,
    )
    run_task(_task(tmp_path), executor=FakeExecutor(outcome, store=tmp_path / "store"))
    assert not (tmp_path / "store" / "geometrias").exists()


def test_missing_mopac_outputs_turn_a_verified_result_into_a_failure(
    tmp_path: Path,
) -> None:
    # No store given: the executor leaves no .out files behind.
    result = run_task(_task(tmp_path), executor=FakeExecutor(_verified()))
    assert result.ok is False and result.state == "failed"
    assert "MOPAC" in (result.error or "")
    assert not (tmp_path / "store" / "geometrias" / "mol_00001").exists()


def test_the_calculation_record_carries_program_versions_and_hostname(
    tmp_path: Path,
) -> None:
    from semi_imperium.persistence import SemiImperiumStore

    run_task(
        _task(
            tmp_path,
            program_versions=(("mopac", "23.2.2"), ("mopac_sha256", "def")),
            hostname="central",
        ),
        executor=FakeExecutor(_verified(), store=tmp_path / "store"),
    )
    (record,) = SemiImperiumStore(tmp_path / "store").list_calculations()
    assert record.provenance.program_versions == {
        "mopac": "23.2.2",
        "mopac_sha256": "def",
    }
    assert record.provenance.hostname == "central"
