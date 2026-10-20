from __future__ import annotations

import json
from pathlib import Path

import pytest

from semi_imperium.conformers import ConformerGeometry
from semi_imperium.domain import CalculationResultData
from semi_imperium.training.final_outputs import (
    OUTPUT_DIRNAME,
    FinalOutputError,
    export_final_outputs,
)

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


def _artifacts(store: Path) -> tuple[str, str]:
    folder = store / "artifacts" / "calc-1" / "mopac"
    folder.mkdir(parents=True)
    (folder / "opt.out").write_text("OPT OUTPUT")
    (folder / "freq.out").write_text("FREQ OUTPUT")
    return "artifacts/calc-1/mopac/opt.out", "artifacts/calc-1/mopac/freq.out"


def _export(store: Path, paths: tuple[str, ...], mol_id: str = "cbs_00002") -> Path:
    return export_final_outputs(
        store_root=store,
        mol_id=mol_id,
        smiles="C",
        hamiltonian="PM7",
        protocol="2",
        calculation_id="calc-1",
        signature="abc123",
        geometry=GEOMETRY,
        result=CalculationResultData(
            energy_hof_kcal_mol=-17.9,
            conformer_index=2,
            conformers_evaluated=3,
            lowest_frequency_cm1=1306.0,
            artifact_paths=paths,
        ),
    )


def test_keeps_xyz_both_mopac_outputs_and_the_json(tmp_path: Path) -> None:
    folder = _export(tmp_path, _artifacts(tmp_path))
    assert folder == tmp_path / OUTPUT_DIRNAME / "cbs_00002"
    assert sorted(p.name for p in folder.iterdir()) == [
        "cbs_00002.freq.out",
        "cbs_00002.json",
        "cbs_00002.out",
        "cbs_00002.xyz",
    ]
    assert (folder / "cbs_00002.out").read_text() == "OPT OUTPUT"
    assert (folder / "cbs_00002.freq.out").read_text() == "FREQ OUTPUT"
    xyz = (folder / "cbs_00002.xyz").read_text().splitlines()
    assert xyz[0] == "5" and len(xyz) == 7


def test_json_names_the_chosen_conformer_hof_and_lowest_frequency(
    tmp_path: Path,
) -> None:
    folder = _export(tmp_path, _artifacts(tmp_path))
    data = json.loads((folder / "cbs_00002.json").read_text())
    assert data["conformer_index"] == 2
    assert data["conformers_evaluated"] == 3
    assert data["hof_kcal_mol"] == -17.9
    assert data["lowest_frequency_cm1"] == 1306.0
    assert data["protocol"] == "2" and data["hamiltonian"] == "PM7"
    assert data["smiles"] == "C" and data["calculation_id"] == "calc-1"
    assert data["signature"] == "abc123"


def test_a_second_export_replaces_the_first(tmp_path: Path) -> None:
    paths = _artifacts(tmp_path)
    _export(tmp_path, paths)
    (tmp_path / paths[0]).write_text("NEWER")
    folder = _export(tmp_path, paths)
    assert (folder / "cbs_00002.out").read_text() == "NEWER"
    assert not list((tmp_path / OUTPUT_DIRNAME).glob("*.tmp*"))


def test_a_missing_mopac_output_is_an_error_and_leaves_nothing(
    tmp_path: Path,
) -> None:
    with pytest.raises(FinalOutputError, match="MOPAC"):
        _export(tmp_path, ("artifacts/calc-1/mopac/opt.out",))
    assert not (tmp_path / OUTPUT_DIRNAME / "cbs_00002").exists()


def test_an_unsafe_mol_id_is_refused(tmp_path: Path) -> None:
    with pytest.raises(FinalOutputError, match="mol_id"):
        _export(tmp_path, _artifacts(tmp_path), mol_id="../escape")
