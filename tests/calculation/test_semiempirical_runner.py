from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from grimperium.calculation.contracts.enums import (
    OverallStatus,
    PropertyRole,
    StageExecutionStatus,
)
from grimperium.calculation.runners.semiempirical_runner import (
    SemiempiricalFormationEnthalpyRunner,
    _generate_initial_xyz,
)
from grimperium.crest_pm7.config import CRESTStatus, MOPACStatus, PM7Config
from grimperium.crest_pm7.conformer_generator import CRESTResult
from grimperium.crest_pm7.mopac_optimizer import MOPACResult


def _write_xyz(path: Path, title: str = "ethanol") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "3\n" f"{title}\n" "C 0.0 0.0 0.0\n" "O 1.0 0.0 0.0\n" "H 0.0 1.0 0.0\n",
        encoding="utf-8",
    )
    return path


def test_runner_returns_three_final_estimates_for_successful_method_a(
    tmp_path: Path,
) -> None:
    initial_xyz = _write_xyz(tmp_path / "initial" / "mol-001.xyz")
    preopt_xyz = _write_xyz(tmp_path / "xtb" / "mol-001_preopt.xyz")
    crest_conf = _write_xyz(tmp_path / "crest" / "mol-001_conf000.xyz")
    mopac_hof = {"AM1": -48.1, "PM3": -45.2, "PM7": -42.3}

    def fake_geometry(smiles: str, output_xyz: Path, name: str) -> Path:
        assert smiles == "CCO"
        assert name == "ethanol"
        output_xyz.parent.mkdir(parents=True, exist_ok=True)
        output_xyz.write_text(initial_xyz.read_text(encoding="utf-8"), encoding="utf-8")
        return output_xyz

    def fake_xtb(input_xyz: Path, work_dir: Path, mol_id: str) -> SimpleNamespace:
        assert input_xyz.name == "mol-001_initial.xyz"
        assert mol_id == "mol-001"
        return SimpleNamespace(
            success=True,
            output_xyz=preopt_xyz,
            error_message="",
            time_seconds=1.5,
        )

    def fake_crest(
        mol_id: str,
        input_xyz: Path,
        config: PM7Config,
        smiles: str | None = None,
    ) -> CRESTResult:
        assert mol_id == "mol-001"
        assert input_xyz == preopt_xyz
        assert smiles == "CCO"
        return CRESTResult(
            status=CRESTStatus.SUCCESS,
            conformers_found=1,
            conformer_files=[crest_conf],
            execution_time=2.5,
            work_dir=tmp_path / "crest",
        )

    def fake_mopac(**kwargs: object) -> MOPACResult:
        hamiltonian = str(kwargs["hamiltonian"])
        output_file = tmp_path / "mopac" / f"{hamiltonian.lower()}.out"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text("mopac output", encoding="utf-8")
        return MOPACResult(
            status=MOPACStatus.SUCCESS,
            hof=mopac_hof[hamiltonian],
            execution_time=3.0,
            output_file=output_file,
        )

    runner = SemiempiricalFormationEnthalpyRunner(
        config=PM7Config(temp_dir=tmp_path / "runtime"),
        work_root=tmp_path / "runs",
        geometry_generator=fake_geometry,
        xtb_preoptimizer=fake_xtb,
        crest_runner=fake_crest,
        mopac_runner=fake_mopac,
        descriptor_extractor=lambda path: {"mopac_homo_ev": -10.0},
    )

    result = runner.calculate_single_smiles(
        "CCO",
        molecule_id="mol-001",
        name="ethanol",
    )

    assert result.overall_status is OverallStatus.SUCCESS
    assert result.molecule.smiles == "CCO"
    assert result.molecule.name == "ethanol"
    assert result.run.method_ref.method_id == "semiempirical_am1_pm3_pm7"
    assert result.conformers[0].conformer_index == 0
    assert [estimate.hamiltonian for estimate in result.estimates] == [
        "AM1",
        "PM3",
        "PM7",
    ]
    assert {estimate.role for estimate in result.estimates} == {PropertyRole.FINAL}
    assert [estimate.value.value for estimate in result.estimates] == [
        -48.1,
        -45.2,
        -42.3,
    ]


def test_runner_uses_lowest_crest_conformer_and_passes_charge_state(
    tmp_path: Path,
) -> None:
    selected_conf = _write_xyz(tmp_path / "crest" / "mol-002_conf000.xyz")
    ignored_conf = _write_xyz(tmp_path / "crest" / "mol-002_conf001.xyz")
    mopac_calls: list[dict[str, object]] = []

    def fake_geometry(smiles: str, output_xyz: Path, name: str) -> Path:
        return _write_xyz(output_xyz, name)

    def fake_crest(
        mol_id: str,
        input_xyz: Path,
        config: PM7Config,
        smiles: str | None = None,
    ) -> CRESTResult:
        return CRESTResult(
            status=CRESTStatus.SUCCESS,
            conformers_found=2,
            conformer_files=[selected_conf, ignored_conf],
            execution_time=2.5,
        )

    def fake_mopac(**kwargs: object) -> MOPACResult:
        mopac_calls.append(kwargs)
        hamiltonian = str(kwargs["hamiltonian"])
        return MOPACResult(
            status=MOPACStatus.SUCCESS,
            hof={"AM1": -1.0, "PM3": -2.0, "PM7": -3.0}[hamiltonian],
            output_file=tmp_path / f"{hamiltonian.lower()}.out",
        )

    runner = SemiempiricalFormationEnthalpyRunner(
        config=PM7Config(temp_dir=tmp_path / "runtime"),
        work_root=tmp_path / "runs",
        xtb_enabled=False,
        geometry_generator=fake_geometry,
        crest_runner=fake_crest,
        mopac_runner=fake_mopac,
        descriptor_extractor=lambda path: {},
    )

    runner.calculate_single_smiles(
        "O",
        molecule_id="mol-002",
        charge=-1,
        multiplicity=3,
    )

    assert [call["xyz_file"] for call in mopac_calls] == [
        selected_conf,
        selected_conf,
        selected_conf,
    ]
    assert [call["hamiltonian"] for call in mopac_calls] == ["AM1", "PM3", "PM7"]
    assert {call["charge"] for call in mopac_calls} == {-1}
    assert {call["multiplicity"] for call in mopac_calls} == {3}
    assert ignored_conf not in {call["xyz_file"] for call in mopac_calls}


def test_runner_preserves_successful_estimates_when_one_hamiltonian_fails(
    tmp_path: Path,
) -> None:
    selected_conf = _write_xyz(tmp_path / "crest" / "mol-003_conf000.xyz")

    def fake_geometry(smiles: str, output_xyz: Path, name: str) -> Path:
        return _write_xyz(output_xyz, name)

    def fake_xtb(input_xyz: Path, work_dir: Path, mol_id: str) -> SimpleNamespace:
        return SimpleNamespace(
            success=True,
            output_xyz=input_xyz,
            error_message="",
            time_seconds=1.25,
        )

    def fake_crest(
        mol_id: str,
        input_xyz: Path,
        config: PM7Config,
        smiles: str | None = None,
    ) -> CRESTResult:
        return CRESTResult(
            status=CRESTStatus.SUCCESS,
            conformers_found=1,
            conformer_files=[selected_conf],
            execution_time=2.5,
        )

    def fake_mopac(**kwargs: object) -> MOPACResult:
        hamiltonian = str(kwargs["hamiltonian"])
        if hamiltonian == "PM3":
            return MOPACResult(
                status=MOPACStatus.ERROR,
                hof=None,
                error_message="PM3 failed",
            )
        return MOPACResult(
            status=MOPACStatus.SUCCESS,
            hof={"AM1": -10.0, "PM7": -30.0}[hamiltonian],
            execution_time=3.0,
            output_file=tmp_path / f"{hamiltonian.lower()}.out",
        )

    runner = SemiempiricalFormationEnthalpyRunner(
        config=PM7Config(temp_dir=tmp_path / "runtime"),
        work_root=tmp_path / "runs",
        geometry_generator=fake_geometry,
        xtb_preoptimizer=fake_xtb,
        crest_runner=fake_crest,
        mopac_runner=fake_mopac,
        descriptor_extractor=lambda path: {"mopac_homo_ev": -9.0},
    )

    result = runner.calculate_single_smiles("CCO", molecule_id="mol-003")

    assert result.overall_status is OverallStatus.PARTIAL
    assert [estimate.hamiltonian for estimate in result.estimates] == ["AM1", "PM7"]
    hamiltonian_results = result.conformers[0].hamiltonian_results
    assert [item.hamiltonian for item in hamiltonian_results] == ["AM1", "PM3", "PM7"]
    assert hamiltonian_results[1].status is OverallStatus.FAILED
    assert hamiltonian_results[1].energy_hof is None
    assert hamiltonian_results[1].error_message == "PM3 failed"

    stage_statuses = {stage.stage_id: stage.status for stage in result.stage_executions}
    assert stage_statuses["xtb_preoptimization"] is StageExecutionStatus.SUCCESS
    assert stage_statuses["crest_conformer_search"] is StageExecutionStatus.SUCCESS
    assert stage_statuses["mopac_am1"] is StageExecutionStatus.SUCCESS
    assert stage_statuses["mopac_pm3"] is StageExecutionStatus.FAILED
    assert stage_statuses["mopac_pm7"] is StageExecutionStatus.SUCCESS


def test_default_xtb_adapter_uses_project_preoptimizer_argument_order(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    selected_conf = _write_xyz(tmp_path / "crest" / "mol-004_conf000.xyz")
    preopt_xyz = _write_xyz(tmp_path / "xtb" / "mol-004_preopt.xyz")

    class FakeOptimizer:
        def __init__(self, *, enabled: bool) -> None:
            assert enabled is True

        def preoptimize_structure(
            self,
            mol_id: str,
            xyz_path: Path,
            work_dir: Path,
        ) -> SimpleNamespace:
            assert mol_id == "mol-004"
            assert xyz_path.name == "mol-004_initial.xyz"
            assert work_dir.name == "xtb"
            return SimpleNamespace(
                success=True,
                output_xyz=preopt_xyz,
                error_message="",
                time_seconds=1.0,
            )

    monkeypatch.setattr(
        "grimperium.calculation.runners.semiempirical_runner.xTBPreOptimizer",
        FakeOptimizer,
    )

    def fake_geometry(smiles: str, output_xyz: Path, name: str) -> Path:
        return _write_xyz(output_xyz, name)

    def fake_crest(
        mol_id: str,
        input_xyz: Path,
        config: PM7Config,
        smiles: str | None = None,
    ) -> CRESTResult:
        assert input_xyz == preopt_xyz
        return CRESTResult(
            status=CRESTStatus.SUCCESS,
            conformers_found=1,
            conformer_files=[selected_conf],
            execution_time=2.5,
        )

    def fake_mopac(**kwargs: object) -> MOPACResult:
        hamiltonian = str(kwargs["hamiltonian"])
        return MOPACResult(
            status=MOPACStatus.SUCCESS,
            hof={"AM1": -1.0, "PM3": -2.0, "PM7": -3.0}[hamiltonian],
        )

    runner = SemiempiricalFormationEnthalpyRunner(
        config=PM7Config(temp_dir=tmp_path / "runtime"),
        work_root=tmp_path / "runs",
        geometry_generator=fake_geometry,
        crest_runner=fake_crest,
        mopac_runner=fake_mopac,
        descriptor_extractor=lambda path: {},
    )

    result = runner.calculate_single_smiles("O", molecule_id="mol-004")

    assert result.overall_status is OverallStatus.SUCCESS


def _make_fake_rdkit(
    *,
    invalid_smiles: bool = False,
    embed_status: int = 0,
    num_atoms: int = 3,
    symbols: tuple[str, ...] = ("C", "O", "H"),
    positions: tuple[tuple[float, float, float], ...] | None = None,
) -> tuple[SimpleNamespace, SimpleNamespace]:
    if positions is None:
        positions = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))

    atoms: list[SimpleNamespace] = []
    for idx, symbol in enumerate(symbols[:num_atoms]):
        atoms.append(
            SimpleNamespace(
                GetSymbol=lambda sym=symbol: sym,
                GetIdx=lambda i=idx: i,
            )
        )

    def get_atom_position(atom_idx: int) -> SimpleNamespace:
        x, y, z = positions[atom_idx]
        return SimpleNamespace(x=x, y=y, z=z)

    fake_conf = SimpleNamespace(GetAtomPosition=get_atom_position)
    fake_mol = SimpleNamespace(
        GetNumAtoms=lambda: num_atoms,
        GetAtoms=lambda: atoms,
        GetConformer=lambda: fake_conf,
    )

    fake_chem = SimpleNamespace(
        MolFromSmiles=lambda smiles: None if invalid_smiles else fake_mol,
        AddHs=lambda mol: mol,
    )
    fake_allchem = SimpleNamespace(
        EmbedMolecule=lambda mol, **kwargs: embed_status,
        UFFOptimizeMolecule=lambda mol, maxIters=200: None,
    )
    return fake_chem, fake_allchem


def _patch_rdkit_import(
    monkeypatch: pytest.MonkeyPatch,
    fake_chem: SimpleNamespace,
    fake_allchem: SimpleNamespace,
) -> None:
    def fake_import_module(name: str) -> Any:
        if name == "rdkit.Chem":
            return fake_chem
        if name == "rdkit.Chem.AllChem":
            return fake_allchem
        raise ImportError(f"No module named {name!r}")

    monkeypatch.setattr(
        "grimperium.calculation.runners.semiempirical_runner.import_module",
        fake_import_module,
    )


def test_generate_initial_xyz_raises_for_invalid_smiles(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake_chem, fake_allchem = _make_fake_rdkit(invalid_smiles=True)
    _patch_rdkit_import(monkeypatch, fake_chem, fake_allchem)

    with pytest.raises(ValueError, match="RDKit cannot parse SMILES"):
        _generate_initial_xyz("invalid_smiles", tmp_path / "out.xyz", "bad")


def test_generate_initial_xyz_raises_when_embedding_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake_chem, fake_allchem = _make_fake_rdkit(embed_status=-1)
    _patch_rdkit_import(monkeypatch, fake_chem, fake_allchem)

    with pytest.raises(ValueError, match="RDKit failed to embed"):
        _generate_initial_xyz("CCO", tmp_path / "out.xyz", "ethanol")


def test_generate_initial_xyz_writes_valid_xyz_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake_chem, fake_allchem = _make_fake_rdkit()
    _patch_rdkit_import(monkeypatch, fake_chem, fake_allchem)
    output_xyz = tmp_path / "out.xyz"

    result = _generate_initial_xyz("CCO", output_xyz, "ethanol")

    assert result == output_xyz
    assert output_xyz.is_file()
    lines = output_xyz.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "3"
    assert lines[1] == "ethanol"
    assert len(lines) == 5
    for coord_line in lines[2:]:
        parts = coord_line.split()
        assert len(parts) == 4
        symbol, x, y, z = parts[0], parts[1], parts[2], parts[3]
        assert symbol in {"C", "O", "H"}
        float(x)
        float(y)
        float(z)
