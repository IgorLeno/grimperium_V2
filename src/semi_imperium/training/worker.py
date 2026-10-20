"""One protocol-2 molecule, run in a worker process.

``run_task`` is a plain module-level function over primitive arguments so a
``ProcessPoolExecutor`` can pickle it. The worker persists the Semi-Imperium
calculation record (the full audit trail) but never touches the CSV or the
geometry table: those belong to the coordinator, the only writer.
"""

from __future__ import annotations

import shutil
import time
from dataclasses import dataclass, replace
from pathlib import Path

from semi_imperium.domain import (
    CalculationRecord,
    CalculationState,
    MolecularIdentity,
    VerificationOutcome,
)
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.settings import RuntimeSettings, SemiImperiumSettings
from semi_imperium.training.final_outputs import FinalOutputError, export_final_outputs
from semi_imperium.training.protocol import (
    CANONICAL_PROTOCOL,
    HAMILTONIAN,
    canonical_configuration,
    canonical_settings,
)
from semi_imperium.workflows.calculation import (
    CalculationExecutor,
    ExecutionOutcome,
    ExecutionRequest,
    build_provenance,
)
from semi_imperium.workflows.execution import ARTIFACTS_DIRNAME


@dataclass(frozen=True)
class WorkerTask:
    """Everything one worker needs, as picklable primitives."""

    mol_id: str
    smiles: str
    charge: int
    multiplicity: int
    run_id: str
    store_root: str
    work_dir: str
    crest_executable: str
    xtb_executable: str
    mopac_executable: str
    crest_threads: int
    crest_timeout_s: float
    mopac_timeout_s: float
    keep_artifacts: bool


@dataclass(frozen=True)
class WorkerResult:
    """What the coordinator needs to record one finished molecule."""

    mol_id: str
    smiles: str
    ok: bool
    state: str
    hof_kcal_mol: float | None
    lowest_frequency_cm1: float | None
    conformers_evaluated: int
    conformer_index: int | None
    elements: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    calculation_id: str
    signature: str
    elapsed_s: float
    error: str | None


def _settings(task: WorkerTask) -> SemiImperiumSettings:
    runtime = RuntimeSettings(
        crest_executable=task.crest_executable,
        mopac_executable=task.mopac_executable,
        crest_threads=task.crest_threads,
        crest_timeout_seconds=task.crest_timeout_s,
        mopac_timeout_seconds=task.mopac_timeout_s,
        work_dir=Path(task.work_dir),
        store_root=Path(task.store_root),
    )
    return canonical_settings(SemiImperiumSettings(runtime=runtime))


def _default_executor(
    settings: SemiImperiumSettings, task: WorkerTask
) -> CalculationExecutor:
    from semi_imperium.conformers.crest_runner import SubprocessCrestRunner
    from semi_imperium.workflows.execution import ScientificCalculationExecutor

    runner = SubprocessCrestRunner(
        work_root=Path(task.work_dir) / "crest",
        executable=task.crest_executable,
        xtb_executable=task.xtb_executable,
        threads=task.crest_threads,
        timeout_seconds=task.crest_timeout_s,
    )
    return ScientificCalculationExecutor(settings, crest_runner=runner)


def run_task(
    task: WorkerTask, *, executor: CalculationExecutor | None = None
) -> WorkerResult:
    """Compute one molecule under protocol 2 and persist its audit record."""
    started = time.perf_counter()
    settings = _settings(task)
    configuration = canonical_configuration(settings)
    signature = configuration.signature()
    identity = MolecularIdentity.from_smiles(
        task.smiles, charge=task.charge, multiplicity=task.multiplicity
    )
    store = SemiImperiumStore(task.store_root)
    record = CalculationRecord(
        run_id=task.run_id,
        molecule=identity,
        signature=signature,
        provenance=build_provenance(settings),
    ).transition_to(CalculationState.RUNNING)
    store.save_calculation(record)

    runner = executor or _default_executor(settings, task)
    try:
        outcome = runner.execute(
            ExecutionRequest(
                run_id=task.run_id,
                calculation_id=record.calculation_id,
                identity=identity,
                configuration=configuration,
                hamiltonian=HAMILTONIAN,
            )
        )
    except Exception as exc:  # a crashing executor is a failed result
        outcome = ExecutionOutcome(
            state=CalculationState.FAILED,
            verification=VerificationOutcome.FAILED,
            error_message=f"{type(exc).__name__}: {exc}",
        )

    accepted = (
        outcome.state is CalculationState.VERIFIED
        and outcome.result is not None
        and outcome.result.energy_hof_kcal_mol is not None
        and outcome.geometry is not None
    )
    artifacts = Path(task.store_root) / ARTIFACTS_DIRNAME / record.calculation_id
    if accepted and outcome.result is not None and outcome.geometry is not None:
        # The working directory is disposable; keep the chosen conformer's
        # evidence first. A molecule whose outputs cannot be kept is a failure,
        # not a silent success without files.
        try:
            export_final_outputs(
                store_root=Path(task.store_root),
                mol_id=task.mol_id,
                smiles=task.smiles,
                hamiltonian=HAMILTONIAN,
                protocol=CANONICAL_PROTOCOL,
                calculation_id=record.calculation_id,
                signature=signature.digest,
                geometry=outcome.geometry,
                result=outcome.result,
            )
        except FinalOutputError as exc:
            accepted = False
            outcome = ExecutionOutcome(
                state=CalculationState.FAILED,
                verification=VerificationOutcome.FAILED,
                result=outcome.result,
                error_message=str(exc),
            )
    result_data = outcome.result
    if accepted and not task.keep_artifacts:
        shutil.rmtree(artifacts, ignore_errors=True)
        shutil.rmtree(
            Path(task.work_dir) / "crest" / identity.molecule_id, ignore_errors=True
        )
        if result_data is not None:
            result_data = replace(result_data, artifact_paths=())

    record = record.transition_to(
        outcome.state,
        verification=outcome.verification,
        result=result_data,
        error_message=outcome.error_message,
    )
    store.save_calculation(record)

    geometry = outcome.geometry
    error = None
    if not accepted:
        error = f"{outcome.state.value}: " + (
            outcome.error_message or outcome.verification.value
        )
    return WorkerResult(
        mol_id=task.mol_id,
        smiles=task.smiles,
        ok=accepted,
        state=outcome.state.value,
        hof_kcal_mol=None if result_data is None else result_data.energy_hof_kcal_mol,
        lowest_frequency_cm1=(
            None if result_data is None else result_data.lowest_frequency_cm1
        ),
        conformers_evaluated=(
            0 if result_data is None else result_data.conformers_evaluated
        ),
        conformer_index=None if result_data is None else result_data.conformer_index,
        elements=() if geometry is None else geometry.elements,
        coordinates=() if geometry is None else geometry.coordinates,
        calculation_id=record.calculation_id,
        signature=signature.digest,
        elapsed_s=time.perf_counter() - started,
        error=error,
    )
