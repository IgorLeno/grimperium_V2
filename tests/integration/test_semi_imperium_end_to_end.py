"""End-to-end contracts for the independently usable Semi-Imperium workflow.

These tests cross the package's public boundaries with deterministic adapters.
They prove orchestration and persistence behavior; they are deliberately not
presented as real CREST/MOPAC scientific evidence.  The opt-in real smoke test
lives in ``test_semi_imperium_real_smoke.py``.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import tomllib

from grimperium.cli.app import GrimperiumCLI
from semi_imperium.app import SemiImperiumCLI
from semi_imperium.calculation import SemiImperiumCalculationWorkflow
from semi_imperium.conformers import (
    Conformer,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerRequest,
    ConformerSearchProvenance,
    ConformerWorkflow,
    ConfPassCandidate,
    ConfPassRanking,
    MoleculeTopology,
)
from semi_imperium.domain import (
    CalculationResultData,
    CalculationState,
    ConformerSearchSettings,
    ConformerSelectionSettings,
    ConformerSelectionStrategy,
    ConformerSource,
    EffectiveConfiguration,
    MolecularIdentity,
    MoleculeInputType,
    SemiempiricalSettings,
    VerificationOutcome,
    VerificationPolicy,
    VerificationSettings,
)
from semi_imperium.molecules.service import ResolutionOutcome, ResolutionStatus
from semi_imperium.mopac import (
    CandidateState,
    DisplacementLineage,
    ForceRun,
    JsonWorkflowJournal,
    OptimizationRun,
)
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.settings import RuntimeSettings, SemiImperiumSettings
from semi_imperium.workflows.calculation import (
    CalculateSession,
    ExecutionOutcome,
    ExecutionRequest,
    ReuseAction,
    build_provenance,
    execute_prepared,
    prepare_runs,
    review,
    save_prepared_runs,
)
from semi_imperium.workflows.execution import outcome_from_hamiltonian_result
from semi_imperium.workspace import SemiImperiumWorkspace

ETHANOL = MolecularIdentity(
    canonical_smiles="CCO",
    display_name="ethanol",
    original_input="ethanol",
    input_type=MoleculeInputType.CHEMICAL_NAME,
    resolved_name="Ethanol",
    resolver="integration-resolver",
)
WATER = MolecularIdentity(
    canonical_smiles="O",
    display_name="water",
    original_input="O",
    input_type=MoleculeInputType.SMILES,
    resolver="direct_smiles",
)


@dataclass
class IntegrationResolver:
    """Offline resolver that still preserves the name/SMILES distinction."""

    identities: dict[str, MolecularIdentity] = field(
        default_factory=lambda: {"ethanol": ETHANOL, "O": WATER}
    )

    def resolve(
        self,
        value: str,
        input_type: MoleculeInputType | str,
        *,
        path: Any = None,
        charge: int | None = None,
        multiplicity: int = 1,
    ) -> ResolutionOutcome:
        del path, charge, multiplicity
        parsed_type = MoleculeInputType(input_type)
        identity = self.identities.get(value)
        if identity is None:
            return ResolutionOutcome(
                original_input=value,
                input_type=parsed_type,
                status=(
                    ResolutionStatus.INVALID
                    if parsed_type is MoleculeInputType.SMILES
                    else ResolutionStatus.UNRESOLVED
                ),
                message=(
                    "SMILES failed local validation"
                    if parsed_type is MoleculeInputType.SMILES
                    else "Chemical name did not resolve to a structure"
                ),
            )
        return ResolutionOutcome(
            original_input=value,
            input_type=parsed_type,
            status=ResolutionStatus.RESOLVED,
            identity=identity,
        )


@dataclass
class LifecycleExecutor:
    """Return a different terminal verification outcome per Hamiltonian."""

    requests: list[ExecutionRequest] = field(default_factory=list)

    def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
        self.requests.append(request)
        result = CalculationResultData(
            energy_hof_kcal_mol={"AM1": -52.1, "PM3": -55.0, "PM7": -56.3}[
                request.hamiltonian
            ],
            conformer_index=0,
            conformers_evaluated=2,
        )
        if request.hamiltonian == "AM1":
            return ExecutionOutcome(
                state=CalculationState.VERIFIED,
                verification=VerificationOutcome.MINIMUM_CONFIRMED,
                result=result,
            )
        if request.hamiltonian == "PM3":
            return ExecutionOutcome(
                state=CalculationState.UNVERIFIED,
                verification=VerificationOutcome.FAILED,
                result=result,
                error_message="FORCE output was incomplete",
            )
        return ExecutionOutcome(
            state=CalculationState.SADDLE,
            verification=VerificationOutcome.SADDLE_POINT,
            result=result,
        )


@dataclass
class VerifiedExecutor:
    requests: list[ExecutionRequest] = field(default_factory=list)

    def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
        self.requests.append(request)
        return ExecutionOutcome(
            state=CalculationState.VERIFIED,
            verification=VerificationOutcome.MINIMUM_CONFIRMED,
            result=CalculationResultData(
                energy_hof_kcal_mol=-50.0 - len(self.requests),
                conformer_index=0,
                conformers_evaluated=1,
            ),
        )


def integration_settings(tmp_path: Path) -> SemiImperiumSettings:
    return replace(
        SemiImperiumSettings(),
        verification=VerificationSettings(
            policy=VerificationPolicy.REQUIRE_MINIMUM,
            max_displacement_reoptimizations=0,
        ),
        runtime=RuntimeSettings(
            work_dir=tmp_path / "work",
            store_root=tmp_path / "store",
        ),
    )


def test_multi_molecule_name_and_smiles_session_persists_explicit_outcomes(
    tmp_path: Path,
) -> None:
    settings = integration_settings(tmp_path)
    store = SemiImperiumStore(settings.runtime.store_root)
    executor = LifecycleExecutor()
    session = CalculateSession(
        resolution_service=IntegrationResolver()  # type: ignore[arg-type]
    )
    name_entry = session.add(
        "ethanol",
        MoleculeInputType.CHEMICAL_NAME,
        hamiltonians=("AM1", "PM7"),
        crest_enabled=True,
    )
    smiles_entry = session.add(
        "O",
        MoleculeInputType.SMILES,
        hamiltonians=("PM3",),
        crest_enabled=False,
    )
    invalid_entry = session.add("not-a-smiles", MoleculeInputType.SMILES)
    unresolved_entry = session.add("unknown compound", MoleculeInputType.CHEMICAL_NAME)

    plan = review(session, store=store, settings=settings)

    assert [(item.entry_id, item.identity.input_type) for item in plan.molecules] == [
        (name_entry.entry_id, MoleculeInputType.CHEMICAL_NAME),
        (smiles_entry.entry_id, MoleculeInputType.SMILES),
    ]
    assert [(item.entry_id, item.status) for item in plan.blocked] == [
        (invalid_entry.entry_id, ResolutionStatus.INVALID.value),
        (unresolved_entry.entry_id, ResolutionStatus.UNRESOLVED.value),
    ]
    assert plan.requested_count == 3

    provenance = build_provenance(
        settings,
        hostname="integration-host",
        notes="deterministic adapter contract; not scientific evidence",
        program_versions={"crest": "adapter", "mopac": "adapter"},
    )
    work = prepare_runs(
        plan,
        ReuseAction.CALCULATE_MISSING,
        provenance=provenance,
        run_id_factory=lambda: "run-integration-session",
    )
    paths = save_prepared_runs(work, store)
    summary = execute_prepared(work, store=store, executor=executor)

    assert len(paths) == 6  # three run manifests plus three calculations
    assert summary.executed_count == 3
    assert [
        (request.hamiltonian, request.configuration.conformer_search.enabled)
        for request in executor.requests
    ] == [
        ("AM1", True),
        ("PM7", True),
        ("PM3", False),
    ]
    records = store.list_calculations()
    assert {record.state for record in records} == {
        CalculationState.VERIFIED,
        CalculationState.UNVERIFIED,
        CalculationState.SADDLE,
    }
    assert {record.verification for record in records} == {
        VerificationOutcome.MINIMUM_CONFIRMED,
        VerificationOutcome.FAILED,
        VerificationOutcome.SADDLE_POINT,
    }
    assert {record.provenance.source for record in records} == {"semi_imperium"}
    assert {record.provenance.hostname for record in records} == {"integration-host"}
    for record in records:
        payload = json.loads(store.calculation_path(record).read_text(encoding="utf-8"))
        assert payload["signature"]["digest"] == record.signature.digest
        assert payload["provenance"]["program_versions"] == {
            "crest": "adapter",
            "mopac": "adapter",
        }

    # Under require_minimum, only AM1 is reusable.  A reproducible saddle or a
    # FORCE failure remains valuable provenance, but cannot answer a request
    # that explicitly requires a confirmed local minimum.
    repeated = CalculateSession(
        resolution_service=IntegrationResolver()  # type: ignore[arg-type]
    )
    repeated.add("ethanol", hamiltonians=("AM1", "PM7"), crest_enabled=True)
    repeated.add(
        "O", MoleculeInputType.SMILES, hamiltonians=("PM3",), crest_enabled=False
    )
    repeated_plan = review(repeated, store=store, settings=settings)
    assert [
        item.hamiltonian
        for molecule in repeated_plan.molecules
        for item in molecule.reusable
    ] == ["AM1"]
    assert {
        item.hamiltonian
        for molecule in repeated_plan.molecules
        for item in molecule.missing
    } == {"PM3", "PM7"}


def test_signature_reuse_calculates_only_the_missing_hamiltonian(
    tmp_path: Path,
) -> None:
    settings = integration_settings(tmp_path)
    store = SemiImperiumStore(settings.runtime.store_root)
    resolver = IntegrationResolver()

    first_session = CalculateSession(resolution_service=resolver)  # type: ignore[arg-type]
    first_session.add("ethanol", hamiltonians=("PM7",), crest_enabled=True)
    first_plan = review(first_session, store=store, settings=settings)
    first_work = prepare_runs(
        first_plan,
        ReuseAction.CALCULATE_MISSING,
        provenance=build_provenance(settings),
        run_id_factory=lambda: "run-first",
    )
    first_executor = VerifiedExecutor()
    save_prepared_runs(first_work, store)
    execute_prepared(first_work, store=store, executor=first_executor)
    original = store.list_calculations()[0]
    original_path = store.calculation_path(original)
    original_text = original_path.read_text(encoding="utf-8")

    second_session = CalculateSession(resolution_service=resolver)  # type: ignore[arg-type]
    second_session.add("ethanol", hamiltonians=("AM1", "PM7"), crest_enabled=True)
    second_plan = review(second_session, store=store, settings=settings)

    assert [item.hamiltonian for item in second_plan.molecules[0].reusable] == ["PM7"]
    assert [item.hamiltonian for item in second_plan.molecules[0].missing] == ["AM1"]
    assert (
        second_plan.molecules[0].hamiltonian_plans[0].signature
        != second_plan.molecules[0].hamiltonian_plans[1].signature
    )

    second_work = prepare_runs(
        second_plan,
        ReuseAction.CALCULATE_MISSING,
        provenance=build_provenance(settings),
        run_id_factory=lambda: "run-second",
    )
    second_executor = VerifiedExecutor()
    save_prepared_runs(second_work, store)
    second_summary = execute_prepared(
        second_work, store=store, executor=second_executor
    )

    assert [request.hamiltonian for request in second_executor.requests] == ["AM1"]
    assert second_summary.reused == (original,)
    assert len(store.list_calculations()) == 2
    assert original_path.read_text(encoding="utf-8") == original_text


def water_geometry(offset: float = 0.0) -> ConformerGeometry:
    return ConformerGeometry(
        elements=("O", "H", "H"),
        coordinates=(
            (offset, 0.0, 0.0),
            (offset + 0.96, 0.0, 0.0),
            (offset - 0.24, 0.93, 0.0),
        ),
    )


def crest_ensemble(settings: ConformerSearchSettings) -> ConformerEnsemble:
    return ConformerEnsemble(
        conformers=(
            Conformer(index=0, geometry=water_geometry(), energy_kcal_mol=3.0),
            Conformer(index=1, geometry=water_geometry(0.1), energy_kcal_mol=-1.0),
            Conformer(index=2, geometry=water_geometry(0.2), energy_kcal_mol=0.5),
        ),
        provenance=ConformerSearchProvenance(
            source=ConformerSource.CREST,
            program="crest",
            program_version="integration-adapter",
            settings=settings,
            run_id="run-conformers",
        ),
    )


@dataclass
class RecordingSearch:
    calls: int = 0

    def search(
        self, request: ConformerRequest, settings: ConformerSearchSettings
    ) -> ConformerEnsemble:
        del request
        self.calls += 1
        return crest_ensemble(settings)


@dataclass
class RecordingInitialStructure:
    calls: int = 0

    def build(
        self, request: ConformerRequest, settings: ConformerSearchSettings
    ) -> ConformerEnsemble:
        self.calls += 1
        return ConformerEnsemble(
            conformers=(Conformer(index=0, geometry=water_geometry()),),
            provenance=ConformerSearchProvenance(
                source=ConformerSource.RDKIT_INITIAL_3D,
                program="rdkit",
                program_version="integration-adapter",
                settings=settings,
                run_id=request.run_id,
            ),
        )


@dataclass
class RecordingConfPass:
    candidates: tuple[ConfPassCandidate, ...] = ()

    def prioritize(
        self, candidates: tuple[ConfPassCandidate, ...]
    ) -> tuple[ConfPassRanking, ...]:
        self.candidates = candidates
        priorities = {0: 1, 1: 2, 2: 0}
        return tuple(
            ConfPassRanking(
                index=candidate.index,
                priority=priorities[candidate.index],
                pas_completeness_class=("complete" if candidate.index == 2 else None),
            )
            for candidate in candidates
        )


def force_output(frequency: float) -> str:
    return (
        "START OF FORCE CALCULATION OUTPUT\n"
        "GRADIENTS WERE INITIALLY ACCEPTABLY SMALL\n"
        f"VIBRATION 1 A1 FREQ. {frequency:.3f}\n"
        "VIBRATION 2 A1 FREQ. 56.000\n"
        "== MOPAC DONE ==\n"
    )


@dataclass
class MixedMopacBackend:
    optimize_calls: list[tuple[str, int]] = field(default_factory=list)

    def optimize(
        self,
        *,
        hamiltonian: str,
        geometry: ConformerGeometry,
        source_conformer_index: int,
        attempt_id: str,
        displacement: DisplacementLineage | None,
    ) -> OptimizationRun:
        del attempt_id, displacement
        self.optimize_calls.append((hamiltonian, source_conformer_index))
        energy = {
            "AM1": {1: -10.0, 2: -12.0},
            "PM3": {1: -22.0, 2: -20.0},
            "PM7": {1: -30.0, 2: -31.0},
        }[hamiltonian][source_conformer_index]
        return OptimizationRun(
            converged=True,
            geometry=geometry,
            heat_of_formation_kcal_mol=energy,
        )

    def verify_force(
        self,
        *,
        hamiltonian: str,
        optimization: OptimizationRun,
        attempt_id: str,
    ) -> ForceRun:
        del optimization, attempt_id
        if hamiltonian == "AM1":
            return ForceRun(output=force_output(-2.0))
        if hamiltonian == "PM3":
            return ForceRun(output="FORCE terminated before frequencies")
        return ForceRun(output=force_output(-140.0))


def scientific_configuration(
    *,
    search: ConformerSearchSettings,
    selection: ConformerSelectionSettings,
) -> EffectiveConfiguration:
    return EffectiveConfiguration(
        method_id="semi_imperium_conformer_mopac",
        method_version="integration",
        property_id="standard_enthalpy_of_formation",
        conformer_search=search,
        conformer_selection=selection,
        semiempirical=SemiempiricalSettings(),
        verification=VerificationSettings(
            policy=VerificationPolicy.REQUIRE_MINIMUM,
            max_displacement_reoptimizations=0,
        ),
    )


def test_conformer_boundaries_and_hamiltonian_surfaces_remain_independent(
    tmp_path: Path,
) -> None:
    search = RecordingSearch()
    initial = RecordingInitialStructure()
    confpass = RecordingConfPass()
    conformers = ConformerWorkflow(
        search_backend=search,
        initial_structure_backend=initial,
        confpass_backend=confpass,
    )
    backend = MixedMopacBackend()
    calculation = SemiImperiumCalculationWorkflow(
        conformer_workflow=conformers,
        mopac_backend=backend,
        journal=JsonWorkflowJournal(tmp_path / "workflow.json"),
    )
    request = ConformerRequest(molecule_id="water", smiles="O", run_id="run-conformers")
    energy_configuration = scientific_configuration(
        search=ConformerSearchSettings(enabled=True),
        selection=ConformerSelectionSettings(top_n=2),
    )

    result = calculation.run(request, energy_configuration)

    assert search.calls == 1
    assert initial.calls == 0
    assert result.conformers.selection.selected_indices == (1, 2)
    assert result.conformers.selection.ranking_basis == "crest_search_energy"
    assert backend.optimize_calls == [
        (hamiltonian, conformer_index)
        for hamiltonian in ("AM1", "PM3", "PM7")
        for conformer_index in (1, 2)
    ]
    assert result.minima.for_hamiltonian("AM1").state is CandidateState.MINIMUM_VERIFIED
    assert result.minima.for_hamiltonian("AM1").verified_attempt_id == (
        "am1-attempt-001"
    )
    assert (
        result.minima.for_hamiltonian("PM3").state is CandidateState.VERIFICATION_FAILED
    )
    assert result.minima.for_hamiltonian("PM7").state is CandidateState.SADDLE_DETECTED
    projected = {
        name: outcome_from_hamiltonian_result(
            result.minima.for_hamiltonian(name),
            policy=VerificationPolicy.REQUIRE_MINIMUM,
        )
        for name in ("AM1", "PM3", "PM7")
    }
    assert {name: item.state for name, item in projected.items()} == {
        "AM1": CalculationState.VERIFIED,
        "PM3": CalculationState.UNVERIFIED,
        "PM7": CalculationState.SADDLE,
    }
    journal = json.loads((tmp_path / "workflow.json").read_text(encoding="utf-8"))
    assert set(journal["hamiltonian_results"]) == {"AM1", "PM3", "PM7"}

    no_search = replace(energy_configuration.conformer_search, enabled=False)
    initial_preparation = conformers.prepare(
        request,
        search_settings=no_search,
        selection_settings=ConformerSelectionSettings(top_n=10),
    )
    assert search.calls == 1
    assert initial.calls == 1
    assert initial_preparation.provenance.source is ConformerSource.RDKIT_INITIAL_3D
    assert initial_preparation.selection.selected_indices == (0,)

    experimental = conformers.prepare(
        request,
        search_settings=energy_configuration.conformer_search,
        selection_settings=ConformerSelectionSettings(
            strategy=ConformerSelectionStrategy.CONFPASS_PRIORITIZATION.value,
            top_n=2,
        ),
        topology=MoleculeTopology(
            atom_count=3,
            bonds=((0, 1, 1), (0, 2, 1)),
        ),
    )
    assert [candidate.index for candidate in confpass.candidates] == [0, 1, 2]
    assert experimental.selection.selected_indices == (2, 0)
    assert experimental.selection.is_experimental is True
    assert experimental.selection.advisory_labels == {
        "conf002": "pas_completeness_class=complete"
    }
    assert all(
        "completeness" not in evidence.lower()
        for evidence in experimental.selection.evidence
    )


def test_independent_launch_keeps_the_legacy_grimperium_entrypoint(
    tmp_path: Path,
) -> None:
    settings = integration_settings(tmp_path)
    semi = SemiImperiumCLI(
        SemiImperiumWorkspace(
            settings=settings,
            store=SemiImperiumStore(settings.runtime.store_root),
            session=CalculateSession(
                resolution_service=IntegrationResolver()  # type: ignore[arg-type]
            ),
        )
    )
    legacy = GrimperiumCLI()

    assert set(semi.controller._views) == {
        "calc",
        "methods",
        "databases",
        "settings",
    }
    assert legacy.controller.get_view("models") is not None
    assert legacy.controller.get_view("results") is not None

    project_root = Path(__file__).resolve().parents[2]
    metadata = tomllib.loads((project_root / "pyproject.toml").read_text("utf-8"))
    scripts = metadata["tool"]["poetry"]["scripts"]
    assert scripts["semi-imperium"] == "semi_imperium.app:main"
    assert scripts["grimperium"] == "grimperium.cli:main"

    with (
        patch.object(sys, "argv", ["semi-imperium", "--skip-preflight"]),
        patch("grimperium.utils.logging.setup_logging"),
        patch("semi_imperium.app.SemiImperiumCLI") as cli_class,
    ):
        instance = MagicMock()
        instance.run.return_value = 0
        cli_class.return_value = instance
        from semi_imperium.app import main

        assert main() == 0
    instance.run.assert_called_once_with(skip_preflight=True)
