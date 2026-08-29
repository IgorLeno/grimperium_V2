"""Contract for independent MOPAC minima and bounded saddle recovery."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from grimperium.crest_pm7.config import PM7Config
from semi_imperium.calculation import SemiImperiumCalculationWorkflow
from semi_imperium.conformers import (
    Conformer,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerRequest,
    ConformerSearchProvenance,
    ConformerWorkflow,
    SelectionResult,
)
from semi_imperium.domain import (
    ConformerSearchSettings,
    ConformerSelectionSettings,
    ConformerSelectionStrategy,
    ConformerSource,
    EffectiveConfiguration,
    SemiempiricalSettings,
    VerificationPolicy,
    VerificationSettings,
)
from semi_imperium.mopac import (
    AttemptOrigin,
    CandidateState,
    DisplacementLineage,
    ForceRun,
    JsonWorkflowJournal,
    MopacExecutableBackend,
    MopacExecutableSettings,
    MopacMinimumWorkflow,
    OptimizationRun,
    classify_force_output,
    parse_last_cartesian_geometry,
    parse_normal_coordinate_vectors,
)


def geometry(offset: float = 0.0) -> ConformerGeometry:
    return ConformerGeometry(
        elements=("O", "H", "H"),
        coordinates=(
            (offset, 0.0, 0.0),
            (offset + 0.96, 0.0, 0.0),
            (offset - 0.24, 0.93, 0.0),
        ),
    )


def selection() -> SelectionResult:
    return SelectionResult(
        strategy=ConformerSelectionStrategy.CREST_ENERGY_TOP_N,
        selected=(
            Conformer(index=4, geometry=geometry(), energy_kcal_mol=0.0),
            Conformer(index=9, geometry=geometry(0.2), energy_kcal_mol=0.7),
        ),
        considered=6,
        ranking_basis="crest_search_energy",
    )


def force_output(*frequencies: float) -> str:
    vibrations = "\n".join(
        f"VIBRATION {mode} A1 FREQ. {frequency:.3f}"
        for mode, frequency in enumerate(frequencies, start=1)
    )
    return (
        "START OF FORCE CALCULATION OUTPUT\n"
        "GRADIENTS WERE INITIALLY ACCEPTABLY SMALL\n"
        "GRADIENT NORM = 0.041\n"
        "DESCRIPTION OF VIBRATIONS\n"
        f"{vibrations}\n"
        "== MOPAC DONE ==\n"
    )


def minimum_settings(**overrides: object) -> VerificationSettings:
    fields: dict[str, object] = {
        "policy": VerificationPolicy.REQUIRE_MINIMUM,
    }
    fields.update(overrides)
    return VerificationSettings(**fields)  # type: ignore[arg-type]


class FakeBackend:
    """Records calls while returning method-specific optimization surfaces."""

    energies = {
        "AM1": {0: -13.0, 4: -12.0, 9: -10.0},
        "PM3": {4: -20.0, 9: -24.0},
        "PM7": {4: -31.0, 9: -30.0},
    }

    def __init__(
        self,
        *,
        selected_verdict: str = "minimum",
        recovery_verdict: str = "saddle",
    ) -> None:
        self.selected_verdict = selected_verdict
        self.recovery_verdict = recovery_verdict
        self.optimize_calls: list[dict[str, object]] = []
        self.force_calls: list[dict[str, object]] = []

    def optimize(
        self,
        *,
        hamiltonian: str,
        geometry: ConformerGeometry,
        source_conformer_index: int,
        attempt_id: str,
        displacement: DisplacementLineage | None,
    ) -> OptimizationRun:
        self.optimize_calls.append(
            {
                "hamiltonian": hamiltonian,
                "source_conformer_index": source_conformer_index,
                "attempt_id": attempt_id,
                "displacement": displacement,
                "geometry": geometry,
            }
        )
        base = self.energies[hamiltonian][source_conformer_index]
        recovery_shift = 0.25 if displacement is not None else 0.0
        return OptimizationRun(
            converged=True,
            geometry=geometry,
            heat_of_formation_kcal_mol=base + recovery_shift,
            output_path=f"artifacts/{attempt_id}.out",
        )

    def verify_force(
        self,
        *,
        hamiltonian: str,
        optimization: OptimizationRun,
        attempt_id: str,
    ) -> ForceRun:
        self.force_calls.append(
            {
                "hamiltonian": hamiltonian,
                "attempt_id": attempt_id,
                "optimization": optimization,
            }
        )
        call = next(
            item for item in self.optimize_calls if item["attempt_id"] == attempt_id
        )
        verdict = (
            self.recovery_verdict
            if call["displacement"] is not None
            else self.selected_verdict
        )
        if verdict == "failure":
            return ForceRun(output="MOPAC FORCE CALCULATION FAILED")
        if verdict == "minimum":
            return ForceRun(output=force_output(-2.5, 56.0, 1240.0))
        mode = ((1.0, 0.0, 0.0), (-0.5, 0.0, 0.0), (-0.5, 0.0, 0.0))
        return ForceRun(
            output=force_output(-145.0, 56.0, 1240.0),
            normal_modes={1: mode},
        )


def test_minimum_workflow_rejects_a_policy_that_would_skip_verification() -> None:
    with pytest.raises(ValueError, match="require_minimum"):
        MopacMinimumWorkflow(
            FakeBackend(),
            verification=VerificationSettings(policy=VerificationPolicy.NONE),
        )


def test_each_hamiltonian_optimizes_selected_and_picks_its_own_lowest() -> None:
    backend = FakeBackend()
    result = MopacMinimumWorkflow(backend).run(selection())

    assert [
        (call["hamiltonian"], call["source_conformer_index"])
        for call in backend.optimize_calls
    ] == [
        ("AM1", 4),
        ("AM1", 9),
        ("PM3", 4),
        ("PM3", 9),
        ("PM7", 4),
        ("PM7", 9),
    ]
    assert result.for_hamiltonian("AM1").verified_heat_of_formation_kcal_mol == -12.0
    assert result.for_hamiltonian("PM3").verified_heat_of_formation_kcal_mol == -24.0
    assert result.for_hamiltonian("PM7").verified_heat_of_formation_kcal_mol == -31.0
    assert result.for_hamiltonian("AM1").verified_attempt_id == "am1-attempt-000"
    assert result.for_hamiltonian("PM3").verified_attempt_id == "pm3-attempt-001"
    assert result.selection_lineage.to_dict() == {
        "strategy": "crest_energy_top_n",
        "considered": 6,
        "ranking_basis": "crest_search_energy",
        "selected_indices": [4, 9],
        "evidence": [],
        "experimental": False,
    }
    assert all(
        attempt.state is CandidateState.MINIMUM_VERIFIED
        for outcome in result.hamiltonian_results
        for attempt in outcome.attempts
    )


def test_converged_candidate_is_journaled_unverified_before_force_classification(
    tmp_path: Path,
) -> None:
    journal_path = tmp_path / "mopac-workflow.json"
    workflow = MopacMinimumWorkflow(
        FakeBackend(),
        journal=JsonWorkflowJournal(journal_path),
    )

    workflow.run(selection(), hamiltonians=("AM1",))
    payload = json.loads(journal_path.read_text(encoding="utf-8"))
    events = payload["attempt_events"]

    assert [event["state"] for event in events[:2]] == [
        CandidateState.OPTIMIZED_UNVERIFIED.value,
        CandidateState.OPTIMIZED_UNVERIFIED.value,
    ]
    assert [event["state"] for event in events[2:]] == [
        CandidateState.MINIMUM_VERIFIED.value,
        CandidateState.MINIMUM_VERIFIED.value,
    ]
    assert events[2]["state_history"] == [
        CandidateState.OPTIMIZED_UNVERIFIED.value,
        CandidateState.MINIMUM_VERIFIED.value,
    ]
    assert (
        payload["hamiltonian_results"]["AM1"]["verified_heat_of_formation_kcal_mol"]
        == -12.0
    )


def test_force_parser_ignores_documented_trivial_and_numerical_low_modes() -> None:
    output = (
        "START OF FORCE CALCULATION OUTPUT\n"
        "STATIONARY POINT CONFIRMED\n"
        "NON-LINEAR MOLECULE\n"
        "TRIVIAL MODES INCLUDED\n"
        "ALL 3N FREQUENCIES: -7.0 -4.0 -1.0 0.0 3.0 8.0 44.0 320.0 1600.0\n"
        "FORCE CALCULATION COMPLETED\n"
    )

    classification = classify_force_output(output, VerificationSettings())

    assert classification.state is CandidateState.MINIMUM_VERIFIED
    assert len(classification.diagnostics.trivial_frequencies_cm1) == 6
    assert classification.diagnostics.nontrivial_frequencies_cm1 == (
        44.0,
        320.0,
        1600.0,
    )
    assert classification.diagnostics.imaginary_frequencies_cm1 == ()


def test_projected_low_negative_vibration_is_recorded_as_numerical_not_saddle() -> None:
    classification = classify_force_output(
        force_output(-42.5, 45.0, 1200.0),
        VerificationSettings(),
    )

    assert classification.state is CandidateState.MINIMUM_VERIFIED
    assert classification.diagnostics.numerical_low_frequencies_cm1 == (-42.5,)
    assert classification.diagnostics.frequency_source == (
        "mopac_vibration_descriptions"
    )


def test_ambiguous_all_3n_low_modes_fail_closed_instead_of_hiding_a_saddle() -> None:
    output = (
        "START OF FORCE CALCULATION OUTPUT\n"
        "STATIONARY POINT CONFIRMED\n"
        "NON-LINEAR MOLECULE\n"
        "TRIVIAL MODES INCLUDED\n"
        "ALL 3N FREQUENCIES: -18.0 -7.0 -4.0 -1.0 0.0 3.0 8.0 320.0 1600.0\n"
        "FORCE CALCULATION COMPLETED\n"
    )

    classification = classify_force_output(output, VerificationSettings())

    assert classification.state is CandidateState.VERIFICATION_FAILED
    assert classification.diagnostics.failure_reason is not None
    assert "ambiguous" in classification.diagnostics.failure_reason


def test_real_imaginary_nontrivial_mode_is_a_saddle_and_never_a_final_value() -> None:
    backend = FakeBackend(selected_verdict="saddle")
    result = MopacMinimumWorkflow(
        backend,
        verification=minimum_settings(max_displacement_reoptimizations=0),
    ).run(selection(), hamiltonians=("PM7",))
    outcome = result.for_hamiltonian("PM7")

    assert outcome.state is CandidateState.SADDLE_DETECTED
    assert outcome.verified_heat_of_formation_kcal_mol is None
    assert outcome.verified_attempt_id is None
    assert outcome.provisional_lowest_attempt_id == "pm7-attempt-000"
    assert all(
        attempt.diagnostics is not None
        and attempt.diagnostics.imaginary_frequencies_cm1 == (-145.0,)
        for attempt in outcome.attempts
    )


def test_displacement_reoptimization_is_explicit_and_bounded() -> None:
    backend = FakeBackend(selected_verdict="saddle", recovery_verdict="saddle")
    result = MopacMinimumWorkflow(
        backend,
        verification=minimum_settings(
            max_displacement_reoptimizations=2,
            displacement_step_angstrom=0.08,
        ),
    ).run(selection(), hamiltonians=("PM7",))
    outcome = result.for_hamiltonian("PM7")

    assert outcome.state is CandidateState.SADDLE_DETECTED
    assert outcome.recovery_attempts_used == 2
    assert len(outcome.attempts) == len(selection().selected) + 2
    assert len(backend.optimize_calls) == len(selection().selected) + 2
    assert outcome.verified_heat_of_formation_kcal_mol is None
    recovery = outcome.attempts[-2:]
    assert {attempt.origin for attempt in recovery} == {
        AttemptOrigin.DISPLACEMENT_REOPTIMIZATION
    }
    assert [
        attempt.displacement.direction for attempt in recovery if attempt.displacement
    ] == [
        1,
        -1,
    ]
    assert all(
        attempt.displacement is not None
        and attempt.displacement.amplitude_angstrom == pytest.approx(0.08)
        for attempt in recovery
    )


def test_recovery_promotes_only_the_force_verified_geometry() -> None:
    backend = FakeBackend(selected_verdict="saddle", recovery_verdict="minimum")
    result = MopacMinimumWorkflow(
        backend,
        verification=minimum_settings(max_displacement_reoptimizations=3),
    ).run(selection(), hamiltonians=("AM1",))
    outcome = result.for_hamiltonian("AM1")

    assert outcome.state is CandidateState.MINIMUM_VERIFIED
    assert outcome.recovery_attempts_used == 1
    assert outcome.verified_attempt_id == "am1-attempt-002"
    assert outcome.attempts[-1].state_history == (
        CandidateState.OPTIMIZED_UNVERIFIED,
        CandidateState.MINIMUM_VERIFIED,
    )
    assert outcome.verified_heat_of_formation_kcal_mol == pytest.approx(-11.75)


def test_force_failure_never_invents_a_verified_hof() -> None:
    backend = FakeBackend(selected_verdict="failure")
    result = MopacMinimumWorkflow(
        backend,
        verification=minimum_settings(max_displacement_reoptimizations=5),
    ).run(selection(), hamiltonians=("PM3",))
    outcome = result.for_hamiltonian("PM3")

    assert outcome.state is CandidateState.VERIFICATION_FAILED
    assert outcome.recovery_attempts_used == 0
    assert len(outcome.attempts) == len(selection().selected)
    assert all(
        attempt.state is CandidateState.VERIFICATION_FAILED
        and attempt.diagnostics is not None
        and attempt.diagnostics.failure_reason is not None
        for attempt in outcome.attempts
    )
    assert outcome.verified_heat_of_formation_kcal_mol is None
    assert outcome.to_dict()["verified_heat_of_formation_kcal_mol"] is None


# ---------------------------------------------------------------------------
# Concrete executable adapter (the production MopacMinimumBackend)
# ---------------------------------------------------------------------------


MOPAC_OPTIMIZATION_OUTPUT = """\
 GRADIENTS WERE INITIALLY ACCEPTABLY SMALL
 FINAL HEAT OF FORMATION = -21.500 KCAL/MOL
 CARTESIAN COORDINATES

 NO. ATOM X Y Z
 1 O 0.0100 0.0000 0.0000
 2 H 0.9700 0.0000 0.0000
 3 H -0.2300 0.9300 0.0000

 == MOPAC DONE ==
"""

MOPAC_FORCE_OUTPUT = """\
 START OF FORCE CALCULATION OUTPUT
 THE GRADIENT NORM IS ACCEPTABLY LOW
 NORMAL COORDINATE ANALYSIS
 Root No.
 1 2 3
 1 A1 2 A1 3 A1
 -7.5 56.0 1240.0
 1  0.10  0.00  0.00
 2  0.00  0.10  0.00
 3  0.00  0.00  0.10
 4 -0.05  0.00  0.00
 5  0.00 -0.05  0.00
 6  0.00  0.00 -0.05
 7 -0.05  0.00  0.00
 8  0.00 -0.05  0.00
 9  0.00  0.00 -0.05
 DESCRIPTION OF VIBRATIONS
 VIBRATION 1 A1 ATOM PAIR ENERGY CONTRIBUTION RADIAL
 FREQ. -7.5
 VIBRATION 2 A1 ATOM PAIR ENERGY CONTRIBUTION RADIAL FREQ. 56.0
 VIBRATION 3 A1 ATOM PAIR ENERGY CONTRIBUTION RADIAL FREQ. 1240.0
 == MOPAC DONE ==
"""


class RecordingMopacRunner:
    """Emulates file-producing MOPAC while exercising the real adapter."""

    def __init__(self) -> None:
        self.inputs: list[str] = []

    def run(
        self, argv: list[str], *, cwd: Path, timeout: float
    ) -> subprocess.CompletedProcess[str]:
        del timeout
        input_path = Path(argv[1])
        content = input_path.read_text(encoding="utf-8")
        self.inputs.append(content)
        output = (
            MOPAC_FORCE_OUTPUT
            if " FORCE " in f" {content.splitlines()[0]} "
            else MOPAC_OPTIMIZATION_OUTPUT
        )
        input_path.with_suffix(".out").write_text(output, encoding="utf-8")
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")


@pytest.mark.parametrize("hamiltonian", ["AM1", "PM3", "PM7"])
def test_concrete_backend_runs_optimization_and_force_with_real_artifacts(
    tmp_path: Path, hamiltonian: str
) -> None:
    runner = RecordingMopacRunner()
    backend = MopacExecutableBackend(
        MopacExecutableSettings(
            executable="mopac",
            work_dir=tmp_path / "jobs",
            timeout_seconds=30.0,
            calculation_id="water-run-1",
        ),
        runner=runner,
    )

    optimization = backend.optimize(
        hamiltonian=hamiltonian,
        geometry=geometry(),
        source_conformer_index=4,
        attempt_id=f"{hamiltonian.lower()}-attempt-000",
        displacement=None,
    )
    force = backend.verify_force(
        hamiltonian=hamiltonian,
        optimization=optimization,
        attempt_id=f"{hamiltonian.lower()}-attempt-000",
    )

    assert optimization.converged is True
    assert optimization.heat_of_formation_kcal_mol == pytest.approx(-21.5)
    assert optimization.geometry is not None
    assert optimization.geometry.elements == geometry(0.01).elements
    for actual, expected in zip(
        optimization.geometry.coordinates, geometry(0.01).coordinates
    ):
        assert actual == pytest.approx(expected)
    assert optimization.output_path is not None
    assert Path(optimization.output_path).is_file()
    assert "water-run-1" in Path(optimization.output_path).parts
    assert force.output_path is not None
    assert Path(force.output_path).is_file()
    assert force.execution_error is None
    assert force.normal_modes[1] == (
        (0.1, 0.0, 0.0),
        (-0.05, 0.0, 0.0),
        (-0.05, 0.0, 0.0),
    )
    assert runner.inputs[0].splitlines()[0].startswith(f"{hamiltonian} EF PRECISE")
    assert runner.inputs[1].splitlines()[0] == f"{hamiltonian} FORCE NOREOR"
    assert all(
        line.rstrip().endswith("0") for line in runner.inputs[1].splitlines()[3:]
    )


def test_concrete_backend_rejects_an_unsafe_artifact_namespace(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="calculation_id must be one safe"):
        MopacExecutableSettings(
            executable="mopac",
            work_dir=tmp_path,
            timeout_seconds=30.0,
            calculation_id="../other-calculation",
        )


def test_mopac_output_parsers_read_documented_multiline_sections() -> None:
    parsed = parse_last_cartesian_geometry(
        MOPAC_OPTIMIZATION_OUTPUT,
        expected_elements=("O", "H", "H"),
    )
    assert parsed is not None
    assert parsed.elements == ("O", "H", "H")
    for actual, expected in zip(parsed.coordinates, geometry(0.01).coordinates):
        assert actual == pytest.approx(expected)
    modes = parse_normal_coordinate_vectors(MOPAC_FORCE_OUTPUT, atom_count=3)
    assert set(modes) == {1, 2, 3}

    classification = classify_force_output(
        MOPAC_FORCE_OUTPUT,
        VerificationSettings(imaginary_frequency_threshold_cm1=-10.0),
        atom_count=3,
    )

    assert classification.state is CandidateState.MINIMUM_VERIFIED
    assert classification.diagnostics.numerical_low_frequencies_cm1 == (-7.5,)


class NoSearch:
    def search(self, request: object, settings: object) -> ConformerEnsemble:
        raise AssertionError("disabled CREST route unexpectedly ran")


class OneInitialConformer:
    def build(
        self, request: ConformerRequest, settings: ConformerSearchSettings
    ) -> ConformerEnsemble:
        return ConformerEnsemble(
            conformers=(Conformer(index=0, geometry=geometry()),),
            provenance=ConformerSearchProvenance(
                source=ConformerSource.RDKIT_INITIAL_3D,
                program="rdkit",
                program_version="test",
                settings=settings,
                run_id=request.run_id,
            ),
        )


def test_integrated_calculation_connects_selection_to_persisted_mopac(
    tmp_path: Path,
) -> None:
    journal_path = tmp_path / "integrated-mopac.json"
    runner = RecordingMopacRunner()
    workflow = SemiImperiumCalculationWorkflow.from_pm7_config(
        conformer_workflow=ConformerWorkflow(
            search_backend=NoSearch(),
            initial_structure_backend=OneInitialConformer(),
        ),
        config=PM7Config(temp_dir=tmp_path / "mopac"),
        calculation_id="calc-water-run-1",
        journal_path=journal_path,
        runner=runner,
    )
    configuration = EffectiveConfiguration(
        method_id="semi_imperium_minima",
        method_version="1",
        property_id="standard_enthalpy_of_formation",
        conformer_search=ConformerSearchSettings(enabled=False),
        conformer_selection=ConformerSelectionSettings(top_n=1),
        semiempirical=SemiempiricalSettings(),
        verification=minimum_settings(max_displacement_reoptimizations=0),
    )

    result = workflow.run(
        ConformerRequest(molecule_id="water", smiles="O", run_id="run-1"),
        configuration,
        hamiltonians=("AM1",),
    )

    assert result.conformers.selection.selected_indices == (0,)
    assert result.minima.for_hamiltonian("AM1").state is CandidateState.MINIMUM_VERIFIED
    assert runner.inputs[0].splitlines()[0].startswith("AM1 EF PRECISE")
    assert runner.inputs[1].splitlines()[0] == "AM1 FORCE NOREOR"
    persisted = json.loads(journal_path.read_text(encoding="utf-8"))
    assert persisted["workflow_result"]["selection_lineage"]["selected_indices"] == [0]
