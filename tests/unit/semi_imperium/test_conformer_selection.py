"""Optional CREST search, bounded selection strategies and their limits.

The contract exercised here is:

* CREST is configured independently of AM1/PM3/PM7 and records its
  effective settings, executable version and run provenance;
* disabling the search still leaves a valid initial-3D route for MOPAC;
* Energy Top-N is the default strategy, ranks at the CREST search level
  and defaults to N=10;
* CONFPASS prioritization is experimental, sees the whole ensemble
  before the cut, survives the XYZ-to-SDF adaptation intact, and never
  turns the PAS completeness label into scientific evidence;
* every external program sits behind an adapter, so these tests run
  with in-memory doubles.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import pytest

from semi_imperium.conformers import (
    HARTREE_TO_KCAL_MOL,
    PAS_COMPLETENESS_LABEL_KEY,
    Conformer,
    ConformerBackendError,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerRequest,
    ConformerSearchProvenance,
    ConformerWorkflow,
    ConfPassCandidate,
    ConfPassRanking,
    ConfPassSelector,
    CrestConformerSearch,
    CrestRun,
    EnergyTopNSelector,
    MoleculeTopology,
    SelectionResult,
    UnavailableConfPass,
    parse_crest_ensemble,
    read_sd_record,
    to_sd_record,
)
from semi_imperium.conformers.initial_structure import RDKitInitialStructure
from semi_imperium.domain import (
    DEFAULT_CONFORMER_TOP_N,
    ConformerSearchSettings,
    ConformerSelectionSettings,
    ConformerSelectionStrategy,
    ConformerSource,
    EffectiveConfiguration,
    SemiempiricalSettings,
)
from semi_imperium.domain.configuration import SIGNATURE_VERSION

CREST_XYZ = """\
3
       -76.400000
O    0.000000    0.000000    0.000000
H    0.960000    0.000000    0.000000
H   -0.240000    0.930000    0.000000
3
       -76.399000
O    0.010000    0.000000    0.000000
H    0.970000    0.000000    0.000000
H   -0.230000    0.930000    0.000000
"""

TRUNCATED_XYZ = """\
3
       -76.400000
O    0.000000    0.000000    0.000000
H    0.960000    0.000000    0.000000
"""

WATER_BONDS = ((0, 1, 1), (0, 2, 1))


# ---------------------------------------------------------------------------
# Helpers and test doubles
# ---------------------------------------------------------------------------


def geometry(offset: float = 0.0) -> ConformerGeometry:
    return ConformerGeometry(
        elements=("O", "H", "H"),
        coordinates=(
            (offset, 0.0, 0.0),
            (0.96, 0.0, 0.0),
            (-0.24, 0.93, 0.0),
        ),
    )


def provenance(
    *,
    source: ConformerSource = ConformerSource.CREST,
    program: str = "crest",
    settings: ConformerSearchSettings | None = None,
) -> ConformerSearchProvenance:
    return ConformerSearchProvenance(
        source=source,
        program=program,
        program_version="3.0.1",
        settings=settings if settings is not None else ConformerSearchSettings(),
        run_id="run-0001",
        command=("crest", "input.xyz", "--gfn2"),
    )


def ensemble(energies: Sequence[float | None]) -> ConformerEnsemble:
    conformers = tuple(
        Conformer(
            index=position,
            geometry=geometry(0.01 * position),
            energy_kcal_mol=energy,
        )
        for position, energy in enumerate(energies)
    )
    return ConformerEnsemble(conformers=conformers, provenance=provenance())


def configuration(
    *,
    hamiltonian: str = "PM7",
    search: ConformerSearchSettings | None = None,
    selection: ConformerSelectionSettings | None = None,
) -> EffectiveConfiguration:
    return EffectiveConfiguration(
        method_id="crest_pm7",
        method_version="1.0",
        property_id="standard_enthalpy_of_formation",
        conformer_search=search if search is not None else ConformerSearchSettings(),
        conformer_selection=(
            selection if selection is not None else ConformerSelectionSettings()
        ),
        semiempirical=SemiempiricalSettings(hamiltonian=hamiltonian),
    )


def request(molecule_id: str = "water", smiles: str = "O") -> ConformerRequest:
    return ConformerRequest(molecule_id=molecule_id, smiles=smiles, run_id="run-0001")


@dataclass
class FakeCrestRunner:
    """Stands in for the CREST executable with recorded output."""

    ensemble_xyz: str = CREST_XYZ
    program_version: str = "3.0.1"
    exit_code: int = 0
    stderr: str = ""
    calls: list[str] = field(default_factory=list)

    def run(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> CrestRun:
        self.calls.append(request.molecule_id)
        return CrestRun(
            ensemble_xyz=self.ensemble_xyz,
            program_version=self.program_version,
            command=("crest", "input.xyz", f"--{settings.method}"),
            exit_code=self.exit_code,
            stderr=self.stderr,
        )


class ExplodingSearch:
    """Fails the test if the CREST route runs when it should not."""

    def search(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        raise AssertionError("the CREST search ran while it was disabled")


@dataclass
class StubInitialStructure:
    """Initial-3D route double: one conformer, no energy, no RDKit."""

    calls: int = 0

    def build(
        self,
        request: ConformerRequest,
        settings: ConformerSearchSettings,
    ) -> ConformerEnsemble:
        self.calls += 1
        return ConformerEnsemble(
            conformers=(Conformer(index=0, geometry=geometry(), label="initial"),),
            provenance=provenance(
                source=ConformerSource.RDKIT_INITIAL_3D,
                program="rdkit",
                settings=settings,
            ),
        )


@dataclass
class FakeConfPass:
    """CONFPASS double that records exactly what it was allowed to see."""

    priorities: dict[int, int] = field(default_factory=dict)
    pas_classes: dict[int, str] = field(default_factory=dict)
    drop_index: int | None = None
    seen: list[ConfPassCandidate] = field(default_factory=list)

    def prioritize(
        self,
        candidates: Sequence[ConfPassCandidate],
    ) -> Sequence[ConfPassRanking]:
        self.seen = list(candidates)
        return [
            ConfPassRanking(
                index=candidate.index,
                priority=self.priorities.get(candidate.index, candidate.index),
                pas_completeness_class=self.pas_classes.get(candidate.index),
            )
            for candidate in candidates
            if candidate.index != self.drop_index
        ]


# ---------------------------------------------------------------------------
# 1. CREST is its own concept, with its own provenance
# ---------------------------------------------------------------------------


def test_crest_configuration_is_separate_from_the_hamiltonian_choice() -> None:
    pm7 = configuration(hamiltonian="PM7")
    am1 = configuration(hamiltonian="AM1")

    assert pm7.conformer_search == am1.conformer_search
    assert pm7.conformer_selection == am1.conformer_selection
    assert pm7.signature() != am1.signature()
    assert "hamiltonian" not in pm7.to_dict()["conformer_search"]


def test_disabling_the_crest_search_is_a_signature_relevant_choice() -> None:
    enabled = configuration()
    disabled = configuration(search=ConformerSearchSettings(enabled=False))

    assert enabled.signature() != disabled.signature()

    payload = disabled.to_dict()["conformer_search"]
    assert payload["enabled"] is False
    assert ConformerSearchSettings.from_dict(payload).enabled is False


def test_conformer_contract_change_bumps_the_signature_version() -> None:
    assert SIGNATURE_VERSION == 3
    assert configuration().signature().version == 3


def test_crest_ensemble_records_settings_version_and_run_provenance() -> None:
    runner = FakeCrestRunner()
    settings = ConformerSearchSettings(method="gfnff")
    search = CrestConformerSearch(runner=runner, energy_unit="kcal/mol")

    produced = search.search(request(), settings)
    recorded = produced.provenance.to_dict()

    assert runner.calls == ["water"]
    assert recorded["source"] == ConformerSource.CREST.value
    assert recorded["program_version"] == "3.0.1"
    assert recorded["settings"] == settings.to_dict()
    assert recorded["run_id"] == "run-0001"
    assert recorded["command"] == ["crest", "input.xyz", "--gfnff"]


def test_provenance_refuses_to_leave_the_program_version_blank() -> None:
    with pytest.raises(ValueError, match="program_version must not be empty"):
        ConformerSearchProvenance(
            source=ConformerSource.CREST,
            program="crest",
            program_version="   ",
            settings=ConformerSearchSettings(),
        )


def test_crest_adapter_refuses_to_run_a_disabled_search() -> None:
    search = CrestConformerSearch(runner=FakeCrestRunner())

    with pytest.raises(ConformerBackendError) as failure:
        search.search(request(), ConformerSearchSettings(enabled=False))

    assert failure.value.code == "crest_disabled"


def test_crest_adapter_surfaces_a_failed_execution() -> None:
    runner = FakeCrestRunner(exit_code=2, stderr="segmentation fault")
    search = CrestConformerSearch(runner=runner)

    with pytest.raises(ConformerBackendError) as failure:
        search.search(request(), ConformerSearchSettings())

    assert failure.value.code == "crest_failed"
    assert "segmentation fault" in str(failure.value)


def test_crest_parsing_preserves_order_and_converts_energies() -> None:
    conformers = parse_crest_ensemble(CREST_XYZ)

    assert len(conformers) == 2
    assert conformers[0].geometry.elements == ("O", "H", "H")
    assert conformers[0].geometry.coordinates[1] == (0.96, 0.0, 0.0)

    expected = -76.4 * HARTREE_TO_KCAL_MOL
    assert conformers[0].require_energy() == pytest.approx(expected)
    difference = conformers[1].require_energy() - conformers[0].require_energy()
    assert difference == pytest.approx(0.001 * HARTREE_TO_KCAL_MOL)


def test_malformed_crest_ensemble_is_rejected_instead_of_shortened() -> None:
    with pytest.raises(ConformerBackendError) as failure:
        parse_crest_ensemble(TRUNCATED_XYZ)

    assert failure.value.code == "crest_parse_failed"


# ---------------------------------------------------------------------------
# 2. A disabled search still leaves a valid initial-3D route
# ---------------------------------------------------------------------------


def test_workflow_uses_the_initial_structure_route_when_crest_is_disabled() -> None:
    initial = StubInitialStructure()
    workflow = ConformerWorkflow(
        search_backend=ExplodingSearch(),
        initial_structure_backend=initial,
    )

    prepared = workflow.prepare(
        request(),
        search_settings=ConformerSearchSettings(enabled=False),
        selection_settings=ConformerSelectionSettings(),
    )

    assert initial.calls == 1
    assert prepared.provenance.source is ConformerSource.RDKIT_INITIAL_3D
    assert prepared.provenance.settings.enabled is False
    assert len(prepared.selected) == 1

    basis = prepared.selection.ranking_basis
    assert basis == "single_structure_without_search_energies"


def test_rdkit_initial_structure_builds_a_geometry_without_any_search() -> None:
    settings = ConformerSearchSettings(enabled=False)

    produced = RDKitInitialStructure().build(
        ConformerRequest(molecule_id="ethanol", smiles="CCO"),
        settings,
    )
    conformer = produced.conformers[0]

    assert produced.size == 1
    assert produced.source is ConformerSource.RDKIT_INITIAL_3D
    assert produced.provenance.program_version.strip() != ""
    assert "rdkit.EmbedMolecule" in produced.provenance.command
    assert conformer.atom_count == 9
    assert conformer.energy_kcal_mol is None
    assert any(
        abs(value) > 1e-6
        for position in conformer.geometry.coordinates
        for value in position
    )


def test_rdkit_initial_structure_reports_an_unusable_smiles() -> None:
    with pytest.raises(ConformerBackendError) as failure:
        RDKitInitialStructure().build(
            ConformerRequest(molecule_id="broken", smiles="not-a-smiles"),
            ConformerSearchSettings(enabled=False),
        )

    assert failure.value.code == "initial_structure_parse_failed"


# ---------------------------------------------------------------------------
# 3. Energy Top-N is the default, bounded strategy
# ---------------------------------------------------------------------------


def test_energy_top_n_is_the_default_strategy_with_n_of_ten() -> None:
    settings = ConformerSelectionSettings()

    assert DEFAULT_CONFORMER_TOP_N == 10
    assert settings.top_n == 10
    assert settings.resolved_strategy is ConformerSelectionStrategy.CREST_ENERGY_TOP_N
    assert settings.is_experimental is False


def test_energy_top_n_ranks_the_search_ensemble_and_keeps_ten() -> None:
    energies = [5.0, 1.0, 9.0, 2.0, 8.0, 3.0, 7.0, 4.0, 6.0, 0.0, 10.0, 11.0]
    result = EnergyTopNSelector().select(
        ensemble(energies),
        ConformerSelectionSettings(),
    )

    assert result.considered == 12
    assert len(result.selected) == 10
    assert result.ranking_basis == "crest_search_energy"
    assert result.selected[0].index == 9
    picked = [conformer.require_energy() for conformer in result.selected]
    assert picked == sorted(picked)
    assert max(picked) == 9.0


def test_energy_top_n_uses_the_configured_n() -> None:
    result = EnergyTopNSelector().select(
        ensemble([4.0, 1.0, 3.0, 2.0]),
        ConformerSelectionSettings(top_n=2),
    )

    assert result.selected_indices == (1, 3)
    assert result.considered == 4


def test_energy_top_n_breaks_ties_in_the_search_order() -> None:
    result = EnergyTopNSelector().select(
        ensemble([1.0, 1.0, 1.0]),
        ConformerSelectionSettings(top_n=2),
    )

    assert result.selected_indices == (0, 1)


def test_energy_top_n_applies_an_optional_energy_window() -> None:
    result = EnergyTopNSelector().select(
        ensemble([0.0, 1.5, 9.0, 0.5]),
        ConformerSelectionSettings(top_n=10, energy_window_kcal_mol=2.0),
    )

    assert result.selected_indices == (0, 3, 1)
    assert result.considered == 4


def test_energy_top_n_refuses_an_ensemble_it_cannot_rank() -> None:
    with pytest.raises(ValueError, match="needs an energy for every conformer"):
        EnergyTopNSelector().select(
            ensemble([1.0, None]),
            ConformerSelectionSettings(),
        )


def test_selection_settings_reject_an_unknown_strategy() -> None:
    with pytest.raises(ValueError, match="Unknown conformer selection strategy"):
        ConformerSelectionSettings(strategy="lowest_pm7_hof")


def test_selection_settings_reject_a_non_positive_n() -> None:
    with pytest.raises(ValueError, match="top_n must be >= 1"):
        ConformerSelectionSettings(top_n=0)


def test_selection_settings_reject_an_empty_energy_window() -> None:
    with pytest.raises(ValueError, match="energy_window_kcal_mol must be > 0"):
        ConformerSelectionSettings(energy_window_kcal_mol=0.0)


# ---------------------------------------------------------------------------
# 4. CONFPASS prioritization is experimental and strictly bounded
# ---------------------------------------------------------------------------


def confpass_selection(
    backend: FakeConfPass,
    *,
    conformers: int = 6,
    top_n: int = 2,
) -> SelectionResult:
    energies = [float(index) for index in range(conformers)]
    selector = ConfPassSelector(
        backend=backend,
        topology=MoleculeTopology(atom_count=3, bonds=WATER_BONDS),
        molecule_id="water",
    )
    return selector.select(
        ensemble(energies),
        ConformerSelectionSettings(
            strategy=ConformerSelectionStrategy.CONFPASS_PRIORITIZATION.value,
            top_n=top_n,
        ),
    )


def test_confpass_strategy_is_marked_experimental() -> None:
    settings = ConformerSelectionSettings(
        strategy=ConformerSelectionStrategy.CONFPASS_PRIORITIZATION.value
    )

    assert settings.is_experimental is True

    result = confpass_selection(FakeConfPass())
    assert result.is_experimental is True
    assert "experimental_strategy" in result.evidence


def test_confpass_consumes_the_whole_ensemble_before_the_cut() -> None:
    backend = FakeConfPass(priorities={0: 5, 1: 4, 2: 3, 3: 2, 4: 1, 5: 0})

    result = confpass_selection(backend, conformers=6, top_n=2)

    assert [candidate.index for candidate in backend.seen] == [0, 1, 2, 3, 4, 5]
    assert result.considered == 6
    assert result.selected_indices == (5, 4)


def test_confpass_pas_label_is_advisory_and_never_evidence() -> None:
    backend = FakeConfPass(
        priorities={0: 0, 1: 1, 2: 2},
        pas_classes={0: "complete", 1: "incomplete"},
    )

    result = confpass_selection(backend, conformers=3, top_n=3)

    assert result.advisory_labels == {
        "conf000": f"{PAS_COMPLETENESS_LABEL_KEY}=complete",
        "conf001": f"{PAS_COMPLETENESS_LABEL_KEY}=incomplete",
    }
    assert all("completeness" not in entry.lower() for entry in result.evidence)


def test_pas_completeness_cannot_be_recorded_as_evidence() -> None:
    with pytest.raises(ValueError, match="advisory metadata"):
        SelectionResult(
            strategy=ConformerSelectionStrategy.CONFPASS_PRIORITIZATION,
            selected=ensemble([1.0]).conformers,
            considered=1,
            ranking_basis="confpass_priority",
            evidence=("pas_completeness=complete",),
        )


def test_confpass_rejects_a_ranking_that_skips_conformers() -> None:
    backend = FakeConfPass(drop_index=2)

    with pytest.raises(ConformerBackendError) as failure:
        confpass_selection(backend, conformers=4, top_n=2)

    assert failure.value.code == "confpass_invalid_ranking"


def test_missing_confpass_backend_fails_loudly() -> None:
    selector = ConfPassSelector(
        backend=UnavailableConfPass(),
        topology=MoleculeTopology(atom_count=3, bonds=WATER_BONDS),
    )

    with pytest.raises(ConformerBackendError) as failure:
        selector.select(
            ensemble([1.0, 2.0]),
            ConformerSelectionSettings(
                strategy=ConformerSelectionStrategy.CONFPASS_PRIORITIZATION.value
            ),
        )

    assert failure.value.code == "confpass_unavailable"


def test_xyz_to_sdf_adaptation_preserves_structure_and_provenance() -> None:
    source = ensemble([-3.5, -2.0])
    topology = MoleculeTopology(atom_count=3, bonds=WATER_BONDS)
    conformer = source.conformers[1]

    record = to_sd_record(
        conformer,
        topology,
        source.provenance,
        molecule_id="water",
    )
    adapted = read_sd_record(record)

    assert adapted.geometry.elements == conformer.geometry.elements
    for restored, original in zip(
        adapted.geometry.coordinates, conformer.geometry.coordinates
    ):
        assert restored == pytest.approx(original, abs=1e-4)
    assert adapted.topology == topology
    assert adapted.title == "water conf001"
    assert adapted.data["CONFORMER_INDEX"] == "1"
    assert adapted.data["CONFORMER_SOURCE"] == ConformerSource.CREST.value
    assert adapted.data["SEARCH_PROGRAM_VERSION"] == "3.0.1"
    assert adapted.data["RUN_ID"] == "run-0001"
    assert adapted.data["ENERGY_KCAL_MOL"] == "-2.0"


def test_adaptation_refuses_connectivity_from_another_molecule() -> None:
    source = ensemble([-1.0])

    with pytest.raises(ValueError, match="disagree on the molecule"):
        to_sd_record(
            source.conformers[0],
            MoleculeTopology(atom_count=5),
            source.provenance,
            molecule_id="water",
        )


# ---------------------------------------------------------------------------
# 5. The workflow keeps every external program behind an adapter
# ---------------------------------------------------------------------------


def test_workflow_runs_the_search_then_the_default_strategy() -> None:
    runner = FakeCrestRunner()
    workflow = ConformerWorkflow(
        search_backend=CrestConformerSearch(runner=runner, energy_unit="kcal/mol"),
        initial_structure_backend=StubInitialStructure(),
    )

    prepared = workflow.prepare(
        request(),
        search_settings=ConformerSearchSettings(),
        selection_settings=ConformerSelectionSettings(top_n=1),
    )
    payload = prepared.to_dict()

    assert runner.calls == ["water"]
    assert prepared.ensemble.size == 2
    assert prepared.selected[0].index == 0
    assert prepared.is_experimental is False
    assert payload["ensemble_size"] == 2
    assert payload["selection"]["strategy"] == "crest_energy_top_n"
    assert payload["provenance"]["program"] == "crest"


def test_workflow_requires_a_confpass_backend_and_a_topology() -> None:
    confpass_settings = ConformerSelectionSettings(
        strategy=ConformerSelectionStrategy.CONFPASS_PRIORITIZATION.value
    )
    without_backend = ConformerWorkflow(
        search_backend=CrestConformerSearch(runner=FakeCrestRunner()),
        initial_structure_backend=StubInitialStructure(),
    )

    with pytest.raises(ConformerBackendError) as failure:
        without_backend.prepare(
            request(),
            search_settings=ConformerSearchSettings(),
            selection_settings=confpass_settings,
        )
    assert failure.value.code == "confpass_unavailable"

    with_backend = ConformerWorkflow(
        search_backend=CrestConformerSearch(
            runner=FakeCrestRunner(),
            energy_unit="kcal/mol",
        ),
        initial_structure_backend=StubInitialStructure(),
        confpass_backend=FakeConfPass(),
    )

    with pytest.raises(ValueError, match="needs the molecule topology"):
        with_backend.prepare(
            request(),
            search_settings=ConformerSearchSettings(),
            selection_settings=confpass_settings,
        )
