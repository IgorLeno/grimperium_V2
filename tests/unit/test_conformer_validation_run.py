"""Arms and runner of the validation experiment, on stored CREST fixtures.

No CREST or MOPAC runs here: ensembles come from
``tests/fixtures/crest_atom_order`` and MOPAC is a deterministic double.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from semi_imperium.conformers import ConformerRequest, CrestConformerSearch
from semi_imperium.domain import ConformerSearchSettings
from tests.unit.conformer_validation_support import FIXTURES, FakeMopac, load

arms = load("arms")
run = load("run")

CASES = {
    case["run"]: case["source_order_smiles"]
    for case in json.loads((FIXTURES / "cases.json").read_text())["cases"]
}


def _ensemble(mol_id: str) -> Any:
    runner = arms.StoredEnsembleRunner(FIXTURES)
    return CrestConformerSearch(runner=runner).search(
        ConformerRequest(molecule_id=mol_id, smiles=CASES[mol_id]),
        ConformerSearchSettings(),
    )


def _task(tmp_path: Path, mol_id: str, **overrides: Any) -> Any:
    fields: dict[str, Any] = {
        "mol_id": mol_id,
        "smiles": CASES[mol_id],
        "arms": ("A", "B", "C"),
        "reference_xyz": (FIXTURES / mol_id / "input.xyz").read_text(),
        "h298_cbs": -100.0,
        "work_dir": str(tmp_path / "work"),
        "ensemble_dir": str(FIXTURES),
    }
    fields.update(overrides)
    return run.MoleculeTask(**fields)


def _fake_factory(*broken: str, **per_arm: FakeMopac) -> Any:
    """MOPAC doubles per arm; arms in ``broken`` fail to build a backend."""

    def factory(task: Any, arm: str, mol_dir: Path) -> FakeMopac:
        del task, mol_dir
        if arm in broken:
            raise RuntimeError(f"no MOPAC backend for arm {arm}")
        return per_arm.get(arm, FakeMopac())

    return factory


def test_arm_configurations_differ_only_in_selection() -> None:
    configs = {arm: arms.configuration(arm) for arm in arms.ARM_NAMES}
    assert len({config.signature().digest for config in configs.values()}) == 3
    assert {config.conformer_search for config in configs.values()} == {
        ConformerSearchSettings()
    }
    assert all(c.verification.requires_minimum for c in configs.values())
    folding = configs["C"].conformer_selection.folding_filter
    assert folding.enabled and folding.max_rg_ratio is None
    assert folding.min_topological_distance == 6
    assert not configs["B"].conformer_selection.folding_filter.enabled
    assert arms.parse_arms("c, a") == ("A", "C")
    with pytest.raises(ValueError):
        arms.parse_arms("A,D")


def test_arms_select_from_the_same_stored_ensemble() -> None:
    mol_id = "cbt_013"
    ensemble = _ensemble(mol_id)
    request = ConformerRequest(molecule_id=mol_id, smiles=CASES[mol_id])

    prepared = {arm: arms.prepare(ensemble, request, arm) for arm in arms.ARM_NAMES}

    assert ensemble.size == 73
    for preparation in prepared.values():
        assert preparation.ensemble is ensemble
        assert len(preparation.selected) == 10
    energy_order = [c.index for c in ensemble.ranked_by_energy()[:10]]
    assert list(prepared["A"].selection.selected_indices) == energy_order
    assert prepared["B"].selection.is_experimental
    assert prepared["A"].folding is None and prepared["B"].folding is None
    measurements = prepared["C"].folding.measurements
    assert len(measurements) == 73
    assert {m.index for m in measurements} == {c.index for c in ensemble.conformers}


def test_rigid_molecule_falls_back_to_crest_order_in_confpass_arms() -> None:
    mol_id = "cbt_020"
    ensemble = _ensemble(mol_id)
    request = ConformerRequest(molecule_id=mol_id, smiles=CASES[mol_id])

    prepared = arms.prepare(ensemble, request, "B")

    assert "confpass_fallback_crest_order" in prepared.selection.evidence


def test_process_molecule_writes_one_complete_row_per_arm(tmp_path: Path) -> None:
    task = _task(tmp_path, "cbt_013")

    rows = run.process_molecule(task, _fake_factory())

    assert [row["arm"] for row in rows] == ["A", "B", "C"]
    for row in rows:
        assert row["error"] is None
        assert row["ensemble_size"] == 73
        assert row["crest_source"] == "stored"
        assert row["crest_command"][0] == "stored"
        assert row["final_class"] == "verified"
        assert row["n_optimizations"] == 10
        assert row["n_force"] == 10
        assert len(row["selected_indices"]) == 10
        chosen = row["chosen_conformer_index"]
        assert chosen == min(row["selected_indices"])
        assert row["verified_hof_kcal_mol"] == pytest.approx(-100.0 + 0.1 * chosen)
        assert row["delta_kcal_mol"] == pytest.approx(
            -100.0 - row["verified_hof_kcal_mol"]
        )
        assert row["min_rmsd"] is not None and row["min_rmsd"] >= 0
        assert row["final_rmsd"] == row["rmsd_by_attempt"][row["chosen_attempt_id"]]
        assert len(row["rmsd_by_attempt"]) == 10
        assert row["config_digest"] == arms.configuration(row["arm"]).signature().digest
    by_arm = {row["arm"]: row for row in rows}
    assert by_arm["A"]["overlap_with_A"] == 10
    assert by_arm["A"]["folding"] is None
    assert len(by_arm["C"]["folding"]["measurements"]) == 73
    assert {"fold_contacts", "rg_ratio"} <= set(
        by_arm["C"]["folding"]["measurements"][0]
    )
    assert (tmp_path / "work" / "cbt_013" / "arm_B" / "workflow.json").exists()


def test_failure_in_one_arm_does_not_touch_the_others(tmp_path: Path) -> None:
    task = _task(tmp_path, "cbt_003")

    rows = run.process_molecule(task, _fake_factory("B"))

    by_arm = {row["arm"]: row for row in rows}
    assert by_arm["B"]["error_stage"] == "mopac"
    assert "no MOPAC backend for arm B" in by_arm["B"]["error"]
    assert by_arm["A"]["error"] is None and by_arm["C"]["error"] is None


def test_failed_optimizations_are_a_failed_outcome_not_an_error(
    tmp_path: Path,
) -> None:
    task = _task(tmp_path, "cbt_003", arms=("A",))

    (row,) = run.process_molecule(task, _fake_factory(A=FakeMopac(fail=True)))

    assert row["error"] is None
    assert row["final_class"] == "failed"
    assert row["verified_hof_kcal_mol"] is None and row["delta_kcal_mol"] is None
    assert row["min_rmsd"] is None


def test_missing_ensemble_becomes_an_error_line_per_arm(tmp_path: Path) -> None:
    task = _task(tmp_path, "cbt_003", ensemble_dir=str(tmp_path / "nowhere"))

    rows = run.process_molecule(task, _fake_factory())

    assert [(row["arm"], row["error_stage"]) for row in rows] == [
        ("A", "crest"),
        ("B", "crest"),
        ("C", "crest"),
    ]


def test_run_is_resumable_and_retries_errors_only_on_request(tmp_path: Path) -> None:
    results = tmp_path / "results.jsonl"
    sample = pd.DataFrame(
        {
            "mol_id": ["cbt_003", "cbt_020"],
            "smiles": [CASES["cbt_003"], CASES["cbt_020"]],
        }
    )
    template = _task(tmp_path, "cbt_003", reference_xyz=None, h298_cbs=None)
    tasks = run.build_tasks(
        sample, {}, arms=("A", "B"), done=set(), limit=None, template=template
    )
    assert [task.mol_id for task in tasks] == ["cbt_003", "cbt_020"]

    written = run.run_tasks(
        tasks,
        results,
        workers=1,
        mopac_backend_factory=_fake_factory("B"),
        log=lambda line: None,
    )

    rows = run.read_results(results)
    assert written == len(rows) == 4
    assert all(row["delta_kcal_mol"] is None for row in rows if not row.get("error"))
    done = run.done_pairs(rows, retry_errors=False)
    assert (
        run.build_tasks(
            sample, {}, arms=("A", "B"), done=done, limit=None, template=template
        )
        == []
    )
    retry = run.done_pairs(rows, retry_errors=True)
    pending = run.build_tasks(
        sample, {}, arms=("A", "B"), done=retry, limit=None, template=template
    )
    assert [(task.mol_id, task.arms) for task in pending] == [
        ("cbt_003", ("B",)),
        ("cbt_020", ("B",)),
    ]

    run.run_tasks(
        pending,
        results,
        workers=1,
        mopac_backend_factory=_fake_factory(),
        log=lambda line: None,
    )

    final = run.done_pairs(run.read_results(results), retry_errors=True)
    assert final == {(m, a) for m in ("cbt_003", "cbt_020") for a in ("A", "B")}


def test_limit_counts_molecules_with_pending_arms(tmp_path: Path) -> None:
    sample = pd.DataFrame({"mol_id": list(CASES), "smiles": list(CASES.values())})
    template = _task(tmp_path, "cbt_003")
    done = {("cbt_003", arm) for arm in arms.ARM_NAMES}

    tasks = run.build_tasks(
        sample, {}, arms=arms.ARM_NAMES, done=done, limit=1, template=template
    )

    assert [task.mol_id for task in tasks] == ["cbt_013"]


def test_torn_last_line_is_ignored_but_earlier_corruption_is_not(
    tmp_path: Path,
) -> None:
    path = tmp_path / "results.jsonl"
    path.write_text('{"mol_id": "m", "arm": "A"}\n{"mol_id": "m", "ar')
    assert run.read_results(path) == [{"mol_id": "m", "arm": "A"}]
    path.write_text('{"mol_id": "m", "ar\n{"mol_id": "m", "arm": "A"}\n')
    with pytest.raises(json.JSONDecodeError):
        run.read_results(path)


def test_references_refuse_a_dataset_the_sample_was_not_drawn_from(
    tmp_path: Path,
) -> None:
    dataset = tmp_path / "dataset.csv"
    pd.DataFrame(
        {"mol_id": ["m1", "m2"], "xyz": ["x1", "x2"], "H298_cbs": [-1.0, -2.0]}
    ).to_csv(dataset, index=False)
    sample = pd.DataFrame({"mol_id": ["m2"], "dataset_sha256": ["0" * 64]})

    with pytest.raises(ValueError, match="sha256"):
        run.load_references(dataset, sample)

    sample["dataset_sha256"] = load("sampling").sha256_of(dataset)
    assert run.load_references(dataset, sample) == {"m2": ("x2", -2.0)}


def test_dry_run_plan_counts_crest_searches_and_pairs(tmp_path: Path) -> None:
    template = _task(tmp_path, "cbt_003", ensemble_dir=None)
    sample = pd.DataFrame({"mol_id": list(CASES), "smiles": list(CASES.values())})
    tasks = run.build_tasks(
        sample,
        {},
        arms=arms.ARM_NAMES,
        done={("cbt_003", "A")},
        limit=None,
        template=template,
    )
    finished = tmp_path / "work" / "cbt_013"
    finished.mkdir(parents=True)
    (finished / "crest_run.json").write_text("{}")

    text = run.plan_text(tasks, workers=4)

    assert "molecules with pending arms: 3" in text
    assert "A=2, B=3, C=3" in text
    assert "CREST searches to run: 2" in text


def test_worker_processes_import_the_script_and_isolate_failures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Workers import ``conformer_validation`` by name, as under ``__main__``.
    monkeypatch.syspath_prepend(str(FIXTURES.parents[2] / "scripts"))
    sample = pd.DataFrame({"mol_id": list(CASES), "smiles": list(CASES.values())})
    template = _task(tmp_path, "cbt_003", ensemble_dir=str(tmp_path / "nowhere"))
    tasks = run.build_tasks(
        sample, {}, arms=("A",), done=set(), limit=None, template=template
    )
    results = tmp_path / "results.jsonl"

    written = run.run_tasks(tasks, results, workers=2, log=lambda line: None)

    rows = run.read_results(results)
    assert written == 3
    assert {row["mol_id"] for row in rows} == set(CASES)
    assert {row["error_stage"] for row in rows} == {"crest"}
