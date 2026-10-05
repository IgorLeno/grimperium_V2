"""Contract for the subprocess CREST runner, driven by a fake process."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from semi_imperium.conformers import (
    Conformer,
    ConformerBackendError,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerRequest,
    ConformerSearchProvenance,
    CrestConformerSearch,
)
from semi_imperium.conformers.crest_runner import (
    CREST_CACHE_SETTINGS_MISMATCH,
    CREST_CACHE_UNRECORDED,
    CREST_MISSING_ENSEMBLE,
    CREST_PREOPT_FAILED,
    CREST_TIMEOUT,
    CREST_UNAVAILABLE,
    ENSEMBLE_FILENAME,
    RUN_RECORD_FILENAME,
    SubprocessCrestRunner,
    parse_crest_version,
)
from semi_imperium.domain import ConformerSearchSettings, ConformerSource

WATER = ConformerGeometry(
    elements=("O", "H", "H"),
    coordinates=((0.0, 0.0, 0.0), (0.96, 0.0, 0.0), (-0.24, 0.93, 0.0)),
)
ENSEMBLE = (
    "3\n -5.07\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"
    "3\n -5.06\nO 0.0 0.0 0.1\nH 0.96 0.0 0.1\nH -0.24 0.93 0.1\n"
)
BANNER = "       Version 3.0.2, Mon, 24 November 17:49:33, 11/24/2025\n"


class FixedInitialStructure:
    """Initial-structure double: always the same water geometry."""

    def __init__(self) -> None:
        self.calls = 0

    def build(
        self, request: ConformerRequest, settings: ConformerSearchSettings
    ) -> ConformerEnsemble:
        self.calls += 1
        return ConformerEnsemble(
            conformers=(Conformer(index=0, geometry=WATER, label="initial"),),
            provenance=ConformerSearchProvenance(
                source=ConformerSource.RDKIT_INITIAL_3D,
                program="test",
                program_version="0",
                settings=settings,
            ),
        )


class FakeProcess:
    """Records argv and imitates xTB/CREST by writing their output files."""

    def __init__(
        self,
        *,
        crest_exit: int = 0,
        write_ensemble: bool = True,
        xtb_exit: int = 0,
        raise_on: str | None = None,
        error: Exception | None = None,
    ) -> None:
        self.crest_exit = crest_exit
        self.write_ensemble = write_ensemble
        self.xtb_exit = xtb_exit
        self.raise_on = raise_on
        self.error = error
        self.calls: list[tuple[list[str], Path, float]] = []

    def run(
        self, argv: list[str], *, cwd: Path, timeout: float
    ) -> subprocess.CompletedProcess[str]:
        self.calls.append((argv, cwd, timeout))
        program = Path(argv[0]).name
        if self.error is not None and program == self.raise_on:
            raise self.error
        if program == "xtb":
            if self.xtb_exit == 0:
                (cwd / "xtbopt.xyz").write_text(
                    (cwd / argv[1]).read_text(encoding="utf-8"), encoding="utf-8"
                )
            return subprocess.CompletedProcess(argv, self.xtb_exit, "", "xtb err")
        if self.crest_exit == 0 and self.write_ensemble:
            (cwd / ENSEMBLE_FILENAME).write_text(ENSEMBLE, encoding="utf-8")
            (cwd / "crest_best.xyz").write_text(ENSEMBLE, encoding="utf-8")
        if not self.write_ensemble:
            (cwd / "crest_best.xyz").write_text(ENSEMBLE, encoding="utf-8")
        return subprocess.CompletedProcess(argv, self.crest_exit, BANNER, "boom")


def request(molecule_id: str = "cbs_00001") -> ConformerRequest:
    return ConformerRequest(molecule_id=molecule_id, smiles="O")


def runner(
    tmp_path: Path, process: FakeProcess, **kwargs: object
) -> SubprocessCrestRunner:
    return SubprocessCrestRunner(
        work_root=tmp_path,
        initial_structure=FixedInitialStructure(),
        process_runner=process,
        **kwargs,  # type: ignore[arg-type]
    )


def test_default_settings_build_the_approved_crest_command(tmp_path: Path) -> None:
    process = FakeProcess()
    crest = runner(tmp_path, process, threads=4, timeout_seconds=7200)

    run = crest.run(request(), ConformerSearchSettings())

    xtb_call, crest_call = process.calls
    assert xtb_call[0][:4] == ["xtb", "input.xyz", "--opt", "--gfn2"]
    argv, cwd, timeout = crest_call
    assert argv == [
        "crest",
        "preopt.xyz",
        "--gfn2",
        "--v3",
        "--ewin",
        "6.0",
        "--rthr",
        "0.125",
        "--opt",
        "2",
        "--chrg",
        "0",
        "--uhf",
        "0",
        "--T",
        "4",
    ]
    assert cwd == tmp_path / "cbs_00001"
    assert timeout == 7200
    assert run.exit_code == 0
    assert run.program_version == "3.0.2"
    assert run.command == tuple(argv)
    assert run.ensemble_xyz == ENSEMBLE


def test_quick_mode_nci_and_no_preoptimizer_change_the_command(
    tmp_path: Path,
) -> None:
    process = FakeProcess()
    settings = ConformerSearchSettings(
        quick_mode="squick", nci=True, use_v3=False, preoptimizer="none"
    )

    runner(tmp_path, process).run(request(), settings)

    ((argv, _, _),) = process.calls
    assert argv[:5] == ["crest", "input.xyz", "--gfn2", "--nci", "--squick"]
    assert "--v3" not in argv


def test_input_geometry_keeps_the_initial_structure_atom_order(
    tmp_path: Path,
) -> None:
    runner(tmp_path, FakeProcess()).run(
        request(), ConformerSearchSettings(preoptimizer="none")
    )

    lines = (tmp_path / "cbs_00001" / "input.xyz").read_text().splitlines()
    assert lines[0] == "3"
    assert [line.split()[0] for line in lines[2:]] == ["O", "H", "H"]


def test_finished_run_is_reused_without_running_crest_again(
    tmp_path: Path,
) -> None:
    process = FakeProcess()
    settings = ConformerSearchSettings()
    first = runner(tmp_path, process).run(request(), settings)
    calls = len(process.calls)

    again = runner(tmp_path, process, threads=8).run(request(), settings)

    assert len(process.calls) == calls
    assert again.ensemble_xyz == first.ensemble_xyz
    assert again.command == first.command
    assert again.program_version == "3.0.2"
    record = json.loads((tmp_path / "cbs_00001" / RUN_RECORD_FILENAME).read_text())
    assert record["settings"] == settings.to_dict()
    assert record["elapsed_seconds"] >= 0


def test_stored_run_under_other_settings_is_refused(tmp_path: Path) -> None:
    process = FakeProcess()
    runner(tmp_path, process).run(request(), ConformerSearchSettings())

    with pytest.raises(ConformerBackendError) as error:
        runner(tmp_path, process).run(
            request(), ConformerSearchSettings(energy_window_kcal_mol=3.0)
        )

    assert error.value.code == CREST_CACHE_SETTINGS_MISMATCH


def test_unrecorded_ensemble_is_neither_reused_nor_overwritten(
    tmp_path: Path,
) -> None:
    mol_dir = tmp_path / "cbs_00001"
    mol_dir.mkdir()
    (mol_dir / ENSEMBLE_FILENAME).write_text("stale", encoding="utf-8")
    process = FakeProcess()

    with pytest.raises(ConformerBackendError) as error:
        runner(tmp_path, process).run(request(), ConformerSearchSettings())

    assert error.value.code == CREST_CACHE_UNRECORDED
    assert process.calls == []
    assert (mol_dir / ENSEMBLE_FILENAME).read_text() == "stale"


def test_missing_ensemble_never_falls_back_to_crest_best(tmp_path: Path) -> None:
    process = FakeProcess(write_ensemble=False)

    with pytest.raises(ConformerBackendError) as error:
        runner(tmp_path, process).run(request(), ConformerSearchSettings())

    assert error.value.code == CREST_MISSING_ENSEMBLE
    assert not (tmp_path / "cbs_00001" / RUN_RECORD_FILENAME).exists()


def test_timeout_is_reported_with_a_stable_code(tmp_path: Path) -> None:
    process = FakeProcess(
        raise_on="crest", error=subprocess.TimeoutExpired("crest", 7200)
    )

    with pytest.raises(ConformerBackendError) as error:
        runner(tmp_path, process).run(request(), ConformerSearchSettings())

    assert error.value.code == CREST_TIMEOUT
    assert not (tmp_path / "cbs_00001" / RUN_RECORD_FILENAME).exists()


def test_missing_executable_is_reported_with_a_stable_code(tmp_path: Path) -> None:
    process = FakeProcess(raise_on="crest", error=FileNotFoundError("crest"))

    with pytest.raises(ConformerBackendError) as error:
        runner(tmp_path, process).run(request(), ConformerSearchSettings())

    assert error.value.code == CREST_UNAVAILABLE


def test_failed_preoptimization_stops_before_crest(tmp_path: Path) -> None:
    process = FakeProcess(xtb_exit=1)

    with pytest.raises(ConformerBackendError) as error:
        runner(tmp_path, process).run(request(), ConformerSearchSettings())

    assert error.value.code == CREST_PREOPT_FAILED
    assert [Path(argv[0]).name for argv, _, _ in process.calls] == ["xtb"]


def test_nonzero_exit_surfaces_through_the_search_adapter(tmp_path: Path) -> None:
    process = FakeProcess(crest_exit=3)
    search = CrestConformerSearch(runner=runner(tmp_path, process))

    with pytest.raises(ConformerBackendError) as error:
        search.search(request(), ConformerSearchSettings())

    assert error.value.code == "crest_failed"
    assert "boom" in str(error.value)
    assert not (tmp_path / "cbs_00001" / RUN_RECORD_FILENAME).exists()


def test_search_adapter_parses_the_runner_ensemble(tmp_path: Path) -> None:
    search = CrestConformerSearch(runner=runner(tmp_path, FakeProcess()))

    ensemble = search.search(request(), ConformerSearchSettings())

    assert ensemble.size == 2
    assert ensemble.provenance.program_version == "3.0.2"
    assert ensemble.provenance.command[0] == "crest"


@pytest.mark.parametrize(
    "settings",
    [
        ConformerSearchSettings(method="gfn7"),
        ConformerSearchSettings(quick_mode="fast"),
        ConformerSearchSettings(preoptimizer="mopac"),
    ],
)
def test_unknown_settings_are_refused(
    tmp_path: Path, settings: ConformerSearchSettings
) -> None:
    with pytest.raises(ValueError):
        runner(tmp_path, FakeProcess()).run(request(), settings)


def test_unsafe_molecule_id_cannot_escape_the_work_root(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        runner(tmp_path, FakeProcess()).run(
            request("../outside"), ConformerSearchSettings()
        )


def test_version_parser_reports_unknown_without_a_banner() -> None:
    assert parse_crest_version(BANNER) == "3.0.2"
    assert parse_crest_version("no banner") == "unknown"
