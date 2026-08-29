"""Opt-in scientific smoke test using installed CREST and MOPAC binaries.

This file is excluded from normal scientific claims unless explicitly enabled.
Unlike the deterministic integration suite, it starts both real executables on
one water molecule and leaves all artifacts under pytest's temporary directory.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, replace
from pathlib import Path

import pytest

from grimperium.crest_pm7.config import PM7Config
from semi_imperium.calculation import SemiImperiumCalculationWorkflow
from semi_imperium.conformers import (
    ConformerRequest,
    ConformerWorkflow,
)
from semi_imperium.conformers.crest import CrestConformerSearch, CrestRun
from semi_imperium.conformers.initial_structure import RDKitInitialStructure
from semi_imperium.domain import (
    ConformerSearchSettings,
    ConformerSelectionSettings,
    ConformerSource,
    EffectiveConfiguration,
    SemiempiricalSettings,
    VerificationPolicy,
    VerificationSettings,
)
from semi_imperium.mopac import CandidateState

REAL_SMOKE_ENV = "SEMI_IMPERIUM_REAL_SMOKE"


@dataclass(frozen=True)
class LocalCrestRunner:
    """Minimal real runner used only by the opt-in smoke test."""

    executable: str
    work_root: Path
    timeout_seconds: float = 180.0
    threads: int = 1

    def run(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> CrestRun:
        work_dir = self.work_root / request.molecule_id
        work_dir.mkdir(parents=True, exist_ok=True)
        input_path = work_dir / "input.xyz"
        initial = RDKitInitialStructure().build(
            request, replace(settings, enabled=False)
        )
        input_path.write_text(
            _xyz(
                initial.conformers[0].geometry.elements,
                initial.conformers[0].geometry.coordinates,
            ),
            encoding="utf-8",
        )

        method_flags = {
            "gfn2": "--gfn2",
            "gfnff": "--gfnff",
            "gfn2//gfnff": "--gfn2//gfnff",
        }
        argv = [self.executable, str(input_path), method_flags[settings.method]]
        quick_flags = {
            "quick": "--quick",
            "squick": "--squick",
            "mquick": "--mquick",
        }
        if settings.quick_mode in quick_flags:
            argv.append(quick_flags[settings.quick_mode])
        if settings.use_v3:
            argv.append("--v3")
        if settings.nci:
            argv.append("--nci")
        argv.extend(
            [
                "--ewin",
                str(settings.energy_window_kcal_mol),
                "--rthr",
                str(settings.rmsd_threshold),
                "--opt",
                str(settings.opt_level),
                "--T",
                str(self.threads),
                "--chrg",
                str(request.charge),
                "--uhf",
                str(request.multiplicity - 1),
            ]
        )
        completed = subprocess.run(  # noqa: S603 - explicit opt-in local tool
            argv,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=self.timeout_seconds,
            check=False,
        )
        ensemble = work_dir / "crest_conformers.xyz"
        if not ensemble.is_file():
            ensemble = work_dir / "crest_best.xyz"
        text = ensemble.read_text(encoding="utf-8") if ensemble.is_file() else ""
        return CrestRun(
            ensemble_xyz=text,
            program_version=_crest_version(self.executable),
            command=tuple(argv),
            exit_code=completed.returncode,
            stderr=completed.stderr,
        )


def _xyz(
    elements: tuple[str, ...],
    coordinates: tuple[tuple[float, float, float], ...],
) -> str:
    lines = [str(len(elements)), "RDKit initial geometry for real smoke"]
    lines.extend(
        f"{element:2s} {x:16.10f} {y:16.10f} {z:16.10f}"
        for element, (x, y, z) in zip(elements, coordinates)
    )
    return "\n".join(lines) + "\n"


def _crest_version(executable: str) -> str:
    completed = subprocess.run(  # noqa: S603 - explicit opt-in local tool
        [executable, "--version"],
        capture_output=True,
        text=True,
        timeout=20.0,
        check=False,
    )
    match = re.search(r"\bcrest\s+([0-9][^\s]*)", completed.stdout, re.IGNORECASE)
    return match.group(1) if match else "unknown"


@pytest.mark.integration
def test_real_water_crest_to_mopac_smoke(tmp_path: Path) -> None:
    if os.environ.get(REAL_SMOKE_ENV) != "1":
        pytest.skip(f"set {REAL_SMOKE_ENV}=1 to run installed scientific tools")
    crest = shutil.which("crest")
    mopac = shutil.which("mopac")
    missing = [name for name, path in (("CREST", crest), ("MOPAC", mopac)) if not path]
    if missing:
        pytest.skip("real scientific smoke unavailable: " + ", ".join(missing))
    assert crest is not None
    assert mopac is not None

    search = ConformerSearchSettings(
        enabled=True,
        quick_mode="quick",
        energy_window_kcal_mol=3.0,
        max_conformers=3,
    )
    workflow = SemiImperiumCalculationWorkflow.from_pm7_config(
        conformer_workflow=ConformerWorkflow(
            search_backend=CrestConformerSearch(
                runner=LocalCrestRunner(crest, tmp_path / "crest")
            ),
            initial_structure_backend=RDKitInitialStructure(),
        ),
        config=PM7Config(
            crest_executable=crest,
            mopac_executable=mopac,
            temp_dir=tmp_path / "mopac-temp",
            crest_threads=1,
            crest_timeout=180.0,
            mopac_timeout_base=180.0,
            mopac_timeout_margin=1.0,
        ),
        calculation_id="real-water-smoke",
        journal_path=tmp_path / "journal.json",
        work_dir=tmp_path / "mopac",
    )
    configuration = EffectiveConfiguration(
        method_id="semi_imperium_conformer_mopac",
        method_version="real-smoke",
        property_id="standard_enthalpy_of_formation",
        conformer_search=search,
        conformer_selection=ConformerSelectionSettings(top_n=1),
        semiempirical=SemiempiricalSettings(hamiltonian="PM7"),
        verification=VerificationSettings(
            policy=VerificationPolicy.REQUIRE_MINIMUM,
            max_displacement_reoptimizations=0,
        ),
    )

    result = workflow.run(
        ConformerRequest(molecule_id="water", smiles="O", run_id="real-smoke"),
        configuration,
        hamiltonians=("PM7",),
    )
    outcome = result.minima.for_hamiltonian("PM7")

    assert result.conformers.provenance.source is ConformerSource.CREST
    assert result.conformers.ensemble.size >= 1
    assert result.conformers.selection.selected_indices
    assert outcome.attempts
    assert outcome.state is not CandidateState.OPTIMIZATION_FAILED
    assert outcome.attempts[0].provisional_heat_of_formation_kcal_mol is not None
    assert outcome.attempts[0].optimization_output_path is not None
    assert Path(outcome.attempts[0].optimization_output_path).is_file()
    print(
        {
            "crest_version": result.conformers.provenance.program_version,
            "crest_ensemble_size": result.conformers.ensemble.size,
            "selected_indices": result.conformers.selection.selected_indices,
            "mopac_state": outcome.state.value,
            "provisional_hof_kcal_mol": (
                outcome.attempts[0].provisional_heat_of_formation_kcal_mol
            ),
            "verified_hof_kcal_mol": outcome.verified_heat_of_formation_kcal_mol,
        }
    )
