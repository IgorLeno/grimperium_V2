"""Run the three arms over the sample and append one JSON line per arm.

Resumption: ``results.jsonl`` is append-only and written as soon as a
molecule finishes; a (mol_id, arm) pair already present is skipped on the
next run (pairs that ended in an error are retried only with
``retry_errors``). CREST itself is resumable through the runner's run
record, so an interrupted molecule never repeats a finished search.

Parallelism: one molecule per worker process, arms in sequence inside
it, so CREST runs once and every arm reuses its ensemble. Any failure
becomes an error line for that molecule/arm; the batch goes on.
"""

from __future__ import annotations

import json
import os
import time
import traceback
from collections import Counter
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from conformer_validation.arms import (
    ARM_NAMES,
    HAMILTONIAN,
    StoredEnsembleRunner,
    configuration,
    prepare,
)
from conformer_validation.metrics import (
    ProbeTemplate,
    best_rmsd,
    reference_skeleton,
)
from conformer_validation.sampling import sha256_of
from grimperium.crest_pm7.config import PM7Config
from semi_imperium.conformers import (
    ConformerEnsemble,
    ConformerPreparation,
    ConformerRequest,
    CrestConformerSearch,
)
from semi_imperium.conformers.crest_runner import SubprocessCrestRunner
from semi_imperium.mopac import (
    CandidateAttempt,
    CandidateState,
    HamiltonianResult,
    JsonWorkflowJournal,
    MopacExecutableBackend,
    MopacMinimumBackend,
    MopacMinimumWorkflow,
)

ROW_VERSION = 1
RUN_ID = "conformer_validation"

#: Mean CREST wall time of the stored ``runs/cbt_*`` ensembles (larger
#: molecules than the sample, so the estimate is pessimistic).
ESTIMATED_CREST_SECONDS = 536.0
#: Median MOPAC PM7 job in ``runs/``; each candidate costs one
#: optimization and one FORCE job.
ESTIMATED_MOPAC_JOB_SECONDS = 0.75
ESTIMATED_CANDIDATES_PER_ARM = 10

_VERIFICATION_STATES = {
    CandidateState.MINIMUM_VERIFIED,
    CandidateState.SADDLE_DETECTED,
    CandidateState.VERIFICATION_FAILED,
}
_OUTCOME_CLASS = {
    CandidateState.MINIMUM_VERIFIED: "verified",
    CandidateState.SADDLE_DETECTED: "saddle",
    CandidateState.VERIFICATION_FAILED: "unverified",
    CandidateState.OPTIMIZED_UNVERIFIED: "unverified",
    CandidateState.OPTIMIZATION_FAILED: "failed",
}


@dataclass(frozen=True)
class MoleculeTask:
    """Everything one worker needs for one molecule; picklable on purpose."""

    mol_id: str
    smiles: str
    arms: tuple[str, ...]
    reference_xyz: str | None
    h298_cbs: float | None
    work_dir: str
    ensemble_dir: str | None = None
    crest_threads: int = 4
    crest_timeout: float = 7200.0
    crest_executable: str = "crest"
    xtb_executable: str = "xtb"
    mopac_executable: str = "mopac"


MopacBackendFactory = Callable[[MoleculeTask, str, Path], MopacMinimumBackend]


def real_mopac_backend(
    task: MoleculeTask, arm: str, mol_dir: Path
) -> MopacMinimumBackend:
    """Real MOPAC backend writing under ``mol_dir/arm_<arm>/mopac``."""
    settings = configuration(arm)
    config = PM7Config(
        mopac_executable=task.mopac_executable,
        temp_dir=mol_dir,
        mopac_scf_threshold=settings.semiempirical.scf_convergence,
    )
    return MopacExecutableBackend.from_pm7_config(
        config,
        calculation_id=f"{task.mol_id}-{arm}",
        work_dir=mol_dir / f"arm_{arm}" / "mopac",
    )


def now() -> str:
    """UTC timestamp for result lines."""
    return datetime.now(UTC).isoformat()


def _error_row(
    task: MoleculeTask, arm: str, stage: str, exc: BaseException
) -> dict[str, Any]:
    return {
        "row_version": ROW_VERSION,
        "mol_id": task.mol_id,
        "arm": arm,
        "finished_at": now(),
        "error_stage": stage,
        "error_code": getattr(exc, "code", type(exc).__name__),
        "error": str(exc) or type(exc).__name__,
        "traceback": "".join(traceback.format_exception(exc))[-4000:],
    }


def _search(
    task: MoleculeTask, request: ConformerRequest
) -> tuple[ConformerEnsemble, dict[str, Any]]:
    """Return the molecule's CREST ensemble and how it was obtained."""
    settings = configuration("A").conformer_search
    if task.ensemble_dir is not None:
        search = CrestConformerSearch(
            runner=StoredEnsembleRunner(Path(task.ensemble_dir))
        )
        ensemble = search.search(request, settings)
        return ensemble, {"crest_source": "stored", "crest_seconds": None}
    runner = SubprocessCrestRunner(
        work_root=Path(task.work_dir),
        executable=task.crest_executable,
        xtb_executable=task.xtb_executable,
        threads=task.crest_threads,
        timeout_seconds=task.crest_timeout,
    )
    cached = runner.run_record(request) is not None
    ensemble = CrestConformerSearch(runner=runner).search(request, settings)
    record = runner.run_record(request) or {}
    return ensemble, {
        "crest_source": "cache" if cached else "run",
        "crest_seconds": record.get("elapsed_seconds"),
    }


def _rmsd_by_attempt(
    attempts: Iterable[CandidateAttempt],
    template: ProbeTemplate | None,
    reference: Any,
) -> dict[str, float]:
    if template is None or reference is None:
        return {}
    rmsds: dict[str, float] = {}
    for attempt in attempts:
        geometry = attempt.optimized_geometry
        if geometry is None:
            continue
        probe = template.skeleton(geometry.elements, geometry.coordinates)
        rmsds[attempt.attempt_id] = best_rmsd(probe, reference)
    return rmsds


def _arm_row(
    task: MoleculeTask,
    arm: str,
    *,
    ensemble: ConformerEnsemble,
    crest: dict[str, Any],
    prepared: ConformerPreparation,
    selection_seconds: float,
    a_indices: Sequence[int] | None,
    result: HamiltonianResult,
    mopac_seconds: float,
    template: ProbeTemplate | None,
    reference: Any,
) -> dict[str, Any]:
    selected = list(prepared.selection.selected_indices)
    attempts = {attempt.attempt_id: attempt for attempt in result.attempts}
    chosen_id = result.verified_attempt_id or result.provisional_lowest_attempt_id
    chosen = attempts.get(chosen_id or "")
    rmsds = _rmsd_by_attempt(result.attempts, template, reference)
    verified_hof = result.verified_heat_of_formation_kcal_mol
    folding = prepared.folding
    return {
        "row_version": ROW_VERSION,
        "mol_id": task.mol_id,
        "arm": arm,
        "finished_at": now(),
        "config_digest": configuration(arm).signature().digest,
        "ensemble_size": ensemble.size,
        "crest_source": crest["crest_source"],
        "crest_seconds": crest["crest_seconds"],
        "crest_command": list(ensemble.provenance.command),
        "crest_version": ensemble.provenance.program_version,
        "selection_seconds": selection_seconds,
        "considered": prepared.selection.considered,
        "ranking_basis": prepared.selection.ranking_basis,
        "selected_indices": selected,
        "evidence": list(prepared.selection.evidence),
        "folding": (
            None
            if folding is None
            else {
                "bypassed": folding.bypassed,
                "discarded_indices": list(folding.discarded_indices),
                "measurements": [item.to_dict() for item in folding.measurements],
            }
        ),
        "overlap_with_A": (
            None if a_indices is None else len(set(selected) & set(a_indices))
        ),
        "final_state": result.state.value,
        "final_class": _OUTCOME_CLASS[result.state],
        "verified_hof_kcal_mol": verified_hof,
        "chosen_attempt_id": chosen_id,
        "chosen_conformer_index": (
            None if chosen is None else chosen.source_conformer_index
        ),
        "chosen_hof_kcal_mol": (
            None if chosen is None else chosen.provisional_heat_of_formation_kcal_mol
        ),
        "attempt_states": dict(
            Counter(attempt.state.value for attempt in result.attempts)
        ),
        "n_optimizations": len(result.attempts),
        "n_force": sum(
            any(state in _VERIFICATION_STATES for state in attempt.state_history)
            for attempt in result.attempts
        ),
        "mopac_seconds": mopac_seconds,
        "rmsd_by_attempt": rmsds,
        "min_rmsd": min(rmsds.values()) if rmsds else None,
        "final_rmsd": rmsds.get(chosen_id or ""),
        "h298_cbs": task.h298_cbs,
        "delta_kcal_mol": (
            None
            if verified_hof is None or task.h298_cbs is None
            else task.h298_cbs - verified_hof
        ),
        "error": None,
    }


def process_molecule(
    task: MoleculeTask,
    mopac_backend_factory: MopacBackendFactory | None = None,
) -> list[dict[str, Any]]:
    """Run every pending arm of one molecule; never raises."""
    factory = mopac_backend_factory or real_mopac_backend
    request = ConformerRequest(
        molecule_id=task.mol_id, smiles=task.smiles, run_id=RUN_ID
    )
    mol_dir = Path(task.work_dir) / task.mol_id
    try:
        ensemble, crest = _search(task, request)
    except Exception as exc:  # noqa: BLE001 - becomes an error line per arm
        return [_error_row(task, arm, "crest", exc) for arm in task.arms]

    try:
        template = ProbeTemplate(task.smiles)
        reference = (
            None
            if task.reference_xyz is None
            else reference_skeleton(task.reference_xyz)
        )
    except Exception as exc:  # noqa: BLE001
        return [_error_row(task, arm, "reference", exc) for arm in task.arms]

    # Arm A's selection is cheap and every other arm is compared with it
    # before any MOPAC job runs.
    prepared: dict[str, tuple[ConformerPreparation, float]] = {}
    a_indices: tuple[int, ...] | None = None
    a_error: Exception | None = None
    try:
        started = time.perf_counter()
        prepared["A"] = (prepare(ensemble, request, "A"), time.perf_counter() - started)
        a_indices = prepared["A"][0].selection.selected_indices
    except Exception as exc:  # noqa: BLE001 - reported on arm A's line
        a_error = exc

    rows: list[dict[str, Any]] = []
    for arm in task.arms:
        stage = "selection"
        try:
            if arm == "A" and a_error is not None:
                raise a_error
            if arm not in prepared:
                started = time.perf_counter()
                prepared[arm] = (
                    prepare(ensemble, request, arm),
                    time.perf_counter() - started,
                )
            preparation, selection_seconds = prepared[arm]
            stage = "mopac"
            arm_dir = mol_dir / f"arm_{arm}"
            arm_dir.mkdir(parents=True, exist_ok=True)
            started = time.perf_counter()
            minima = MopacMinimumWorkflow(
                factory(task, arm, mol_dir),
                verification=configuration(arm).verification,
                journal=JsonWorkflowJournal(arm_dir / "workflow.json"),
            ).run(preparation.selection, hamiltonians=(HAMILTONIAN,))
            mopac_seconds = time.perf_counter() - started
            stage = "metrics"
            rows.append(
                _arm_row(
                    task,
                    arm,
                    ensemble=ensemble,
                    crest=crest,
                    prepared=preparation,
                    selection_seconds=selection_seconds,
                    a_indices=a_indices,
                    result=minima.for_hamiltonian(HAMILTONIAN),
                    mopac_seconds=mopac_seconds,
                    template=template,
                    reference=reference,
                )
            )
        except Exception as exc:  # noqa: BLE001 - isolated to this arm
            rows.append(_error_row(task, arm, stage, exc))
    return rows


def read_results(path: Path) -> list[dict[str, Any]]:
    """Read every JSON line; a torn last line from a crash is ignored."""
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
    lines = path.read_text(encoding="utf-8").splitlines()
    for number, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            if number == len(lines):
                continue
            raise
    return rows


def done_pairs(
    rows: Iterable[dict[str, Any]], *, retry_errors: bool
) -> set[tuple[str, str]]:
    """Pairs that must not run again; the latest line of a pair decides."""
    latest: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        latest[(row["mol_id"], row["arm"])] = row
    return {
        key for key, row in latest.items() if not (retry_errors and row.get("error"))
    }


def append_rows(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    """Append rows and push them to disk before returning."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def build_tasks(
    sample: pd.DataFrame,
    references: dict[str, tuple[str, float]],
    *,
    arms: Sequence[str],
    done: set[tuple[str, str]],
    limit: int | None,
    template: MoleculeTask,
) -> list[MoleculeTask]:
    """One task per molecule that still has a pending arm, in sample order."""
    tasks: list[MoleculeTask] = []
    for mol_id, smiles in zip(sample["mol_id"], sample["smiles"], strict=True):
        pending = tuple(arm for arm in arms if (mol_id, arm) not in done)
        if not pending:
            continue
        reference = references.get(mol_id)
        fields = asdict(template)
        fields.update(
            mol_id=mol_id,
            smiles=smiles,
            arms=pending,
            reference_xyz=None if reference is None else reference[0],
            h298_cbs=None if reference is None else reference[1],
        )
        tasks.append(MoleculeTask(**fields))
        if limit is not None and len(tasks) >= limit:
            break
    return tasks


def load_references(
    dataset: Path, sample: pd.DataFrame
) -> dict[str, tuple[str, float]]:
    """Reference XYZ and H298_cbs of the sampled molecules.

    Raises:
        ValueError: If the sample records a dataset digest that does not
            match ``dataset``.
    """
    recorded = {
        str(value)
        for value in sample.get("dataset_sha256", pd.Series(dtype=str)).dropna()
        if str(value)
    }
    if recorded:
        actual = sha256_of(dataset)
        if recorded != {actual}:
            raise ValueError(
                f"The sample was drawn from dataset sha256 {sorted(recorded)} but "
                f"{dataset} has {actual}"
            )
    frame = pd.read_csv(dataset, usecols=["mol_id", "xyz", "H298_cbs"])
    frame = frame[frame["mol_id"].isin(set(sample["mol_id"]))]
    return {
        str(mol_id): (str(xyz), float(h298))
        for mol_id, xyz, h298 in zip(
            frame["mol_id"], frame["xyz"], frame["H298_cbs"], strict=True
        )
    }


def plan_text(tasks: Sequence[MoleculeTask], *, workers: int) -> str:
    """Dry-run account of what would run and a rough cost."""
    crest_needed = 0
    for task in tasks:
        if task.ensemble_dir is not None:
            continue
        if not (Path(task.work_dir) / task.mol_id / "crest_run.json").exists():
            crest_needed += 1
    pairs = sum(len(task.arms) for task in tasks)
    crest_hours = crest_needed * ESTIMATED_CREST_SECONDS / 3600
    mopac_hours = (
        pairs * ESTIMATED_CANDIDATES_PER_ARM * 2 * ESTIMATED_MOPAC_JOB_SECONDS / 3600
    )
    by_arm = Counter(arm for task in tasks for arm in task.arms)
    return "\n".join(
        [
            f"molecules with pending arms: {len(tasks)}",
            "pending (molecule, arm) pairs: "
            + ", ".join(f"{arm}={by_arm.get(arm, 0)}" for arm in ARM_NAMES),
            f"CREST searches to run: {crest_needed} "
            f"(~{crest_hours:.1f} h serial, ~{crest_hours / max(workers, 1):.1f} h "
            f"with {workers} workers)",
            f"MOPAC (rough, <= {ESTIMATED_CANDIDATES_PER_ARM} candidates x "
            f"opt+FORCE per arm): ~{mopac_hours:.2f} h serial",
            "CONFPASS on very large ensembles is not in this estimate "
            "(4702 conformers took 415 s and 778 MB).",
        ]
    )


def run_tasks(
    tasks: Sequence[MoleculeTask],
    results: Path,
    *,
    workers: int,
    mopac_backend_factory: MopacBackendFactory | None = None,
    log: Callable[[str], None] = print,
) -> int:
    """Run ``tasks``, appending each molecule's rows as it finishes.

    ``workers == 1`` runs in this process (and is the only mode that can
    take a custom MOPAC backend factory, which need not be picklable).
    Returns the number of rows written.
    """
    written = 0

    def record(task: MoleculeTask, rows: list[dict[str, Any]]) -> None:
        nonlocal written
        append_rows(results, rows)
        written += len(rows)
        states = ", ".join(
            f"{row['arm']}={row.get('final_class') or 'error:' + str(row.get('error_code'))}"
            for row in rows
        )
        log(f"[{written}] {task.mol_id}: {states}")

    if workers == 1:
        for task in tasks:
            record(task, process_molecule(task, mopac_backend_factory))
        return written

    if mopac_backend_factory is not None:
        raise ValueError("A custom MOPAC backend factory needs workers == 1")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures: dict[Future[list[dict[str, Any]]], MoleculeTask] = {
            pool.submit(process_molecule, task): task for task in tasks
        }
        for future in as_completed(futures):
            task = futures[future]
            try:
                rows = future.result()
            except Exception as exc:  # noqa: BLE001 - e.g. a killed worker
                rows = [_error_row(task, arm, "worker", exc) for arm in task.arms]
            record(task, rows)
    return written


__all__ = [
    "MoleculeTask",
    "append_rows",
    "build_tasks",
    "done_pairs",
    "load_references",
    "plan_text",
    "process_molecule",
    "read_results",
    "real_mopac_backend",
    "run_tasks",
]
