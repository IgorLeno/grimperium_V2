"""Topology derivation from the SMILES and its use in production.

The topology CONFPASS and the folding filter read must address atoms in
the ensemble's own order. These tests pin the SMILES-derived topology to
the order the initial-3D route embeds, show that a reordered ensemble is
refused with a code instead of silently mislabelled, and check that the
production executor wires the provider and the CONFPASS port in.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from semi_imperium.calculation import SemiImperiumCalculationWorkflow
from semi_imperium.conformers import (
    ConformerBackendError,
    ConformerEnsemble,
    ConformerGeometry,
    ConformerPreparation,
    ConformerRequest,
    ConformerWorkflow,
    ConfPassRanking,
    MoleculeTopology,
)
from semi_imperium.conformers.confpass_port import PortedConfPass
from semi_imperium.conformers.initial_structure import RDKitInitialStructure
from semi_imperium.conformers.topology import (
    TOPOLOGY_ATOM_ORDER_MISMATCH,
    TOPOLOGY_PARSE_FAILED,
    SmilesTopology,
    topology_from_smiles,
)
from semi_imperium.domain import (
    ConformerSearchSettings,
    ConformerSelectionSettings,
    ConformerSelectionStrategy,
    FoldingFilterSettings,
    MolecularIdentity,
)
from semi_imperium.settings import SemiImperiumSettings
from semi_imperium.workflows.calculation import ExecutionRequest
from semi_imperium.workflows.execution import ScientificCalculationExecutor

CONFPASS = ConformerSelectionStrategy.CONFPASS_PRIORITIZATION.value


def request(smiles: str = "CCO") -> ConformerRequest:
    return ConformerRequest(molecule_id="mol", smiles=smiles, run_id="run-1")


def embedded(smiles: str = "CCO") -> ConformerEnsemble:
    return RDKitInitialStructure().build(
        request(smiles), ConformerSearchSettings(enabled=False)
    )


def with_geometry(
    ensemble: ConformerEnsemble, geometry: ConformerGeometry
) -> ConformerEnsemble:
    conformer = replace(ensemble.conformers[0], geometry=geometry)
    return ConformerEnsemble(conformers=(conformer,), provenance=ensemble.provenance)


def swap(geometry: ConformerGeometry, first: int, second: int) -> ConformerGeometry:
    elements = list(geometry.elements)
    coordinates = list(geometry.coordinates)
    elements[first], elements[second] = elements[second], elements[first]
    coordinates[first], coordinates[second] = coordinates[second], coordinates[first]
    return ConformerGeometry(elements=tuple(elements), coordinates=tuple(coordinates))


# ---------------------------------------------------------------------------
# Derivation and the atom-order check
# ---------------------------------------------------------------------------


def test_smiles_topology_follows_the_initial_structure_atom_order() -> None:
    structure = embedded("CCO")

    topology, elements = topology_from_smiles("CCO")

    assert elements == structure.conformers[0].geometry.elements
    assert topology.atom_count == 9
    assert SmilesTopology()(request("CCO"), structure) == topology


def test_aromatic_bonds_are_kekulized_to_integer_orders() -> None:
    topology, _ = topology_from_smiles("c1ccccc1")

    ring = [
        order for first, second, order in topology.bonds if first < 6 and second < 6
    ]
    assert sorted(ring) == [1, 1, 1, 2, 2, 2]


def test_an_unreadable_smiles_is_refused_with_a_code() -> None:
    with pytest.raises(ConformerBackendError) as failure:
        topology_from_smiles("C1CC")
    assert failure.value.code == TOPOLOGY_PARSE_FAILED


def test_an_ensemble_with_other_elements_per_position_is_refused() -> None:
    structure = embedded("CCO")
    # Atoms 1 (C) and 2 (O) trade places: same atom count, wrong order.
    reordered = with_geometry(structure, swap(structure.conformers[0].geometry, 1, 2))

    with pytest.raises(ConformerBackendError) as failure:
        SmilesTopology()(request("CCO"), reordered)
    assert failure.value.code == TOPOLOGY_ATOM_ORDER_MISMATCH


def test_hydrogens_moved_between_carbons_are_caught_by_bond_lengths() -> None:
    structure = embedded("CCO")
    topology, elements = topology_from_smiles("CCO")
    hydrogens_of = {
        heavy: [
            first if second == heavy else second
            for first, second, _ in topology.bonds
            if heavy in (first, second) and "H" in (elements[first], elements[second])
        ]
        for heavy in (0, 1)
    }
    # Same elements everywhere, but one H now sits on the other carbon.
    reordered = with_geometry(
        structure,
        swap(structure.conformers[0].geometry, hydrogens_of[0][0], hydrogens_of[1][0]),
    )

    with pytest.raises(ConformerBackendError) as failure:
        SmilesTopology()(request("CCO"), reordered)
    assert failure.value.code == TOPOLOGY_ATOM_ORDER_MISMATCH
    assert "bond" in str(failure.value)


# ---------------------------------------------------------------------------
# The workflow asks the provider only when something reads connectivity
# ---------------------------------------------------------------------------


class RecordingProvider:
    def __init__(self) -> None:
        self.calls = 0

    def __call__(
        self, request: ConformerRequest, ensemble: ConformerEnsemble
    ) -> MoleculeTopology:
        self.calls += 1
        return SmilesTopology()(request, ensemble)


class OrderPreservingConfPass:
    def prioritize(
        self, candidates: Any
    ) -> tuple[ConfPassRanking, ...]:  # pragma: no cover - single conformer
        return tuple(
            ConfPassRanking(index=item.index, priority=item.index)
            for item in candidates
        )


def rdkit_workflow(provider: RecordingProvider) -> ConformerWorkflow:
    return ConformerWorkflow(
        search_backend=_NoSearch(),
        initial_structure_backend=RDKitInitialStructure(),
        confpass_backend=OrderPreservingConfPass(),
        topology_provider=provider,
    )


class _NoSearch:
    def search(self, request: ConformerRequest, settings: Any) -> ConformerEnsemble:
        raise AssertionError("the CREST route must not run with CREST disabled")


def test_the_default_strategy_never_asks_for_a_topology() -> None:
    provider = RecordingProvider()

    rdkit_workflow(provider).prepare(
        request(),
        search_settings=ConformerSearchSettings(enabled=False),
        selection_settings=ConformerSelectionSettings(),
    )

    assert provider.calls == 0


def test_confpass_and_the_folding_filter_get_the_provided_topology() -> None:
    provider = RecordingProvider()
    workflow = rdkit_workflow(provider)
    off = ConformerSearchSettings(enabled=False)

    confpass = workflow.prepare(
        request(),
        search_settings=off,
        selection_settings=ConformerSelectionSettings(strategy=CONFPASS),
    )
    folded = workflow.prepare(
        request(),
        search_settings=off,
        selection_settings=ConformerSelectionSettings(
            folding_filter=FoldingFilterSettings(enabled=True)
        ),
    )

    assert provider.calls == 2
    assert confpass.selection.ranking_basis == "single_conformer_ensemble"
    assert folded.folding is not None


def test_an_explicit_topology_takes_precedence_over_the_provider() -> None:
    provider = RecordingProvider()
    topology, _ = topology_from_smiles("CCO")

    rdkit_workflow(provider).prepare(
        request(),
        search_settings=ConformerSearchSettings(enabled=False),
        selection_settings=ConformerSelectionSettings(strategy=CONFPASS),
        topology=topology,
    )

    assert provider.calls == 0


# ---------------------------------------------------------------------------
# Production executor wiring
# ---------------------------------------------------------------------------


class _Prepared(Exception):
    """Carries the conformer stage's output out of the stubbed workflow."""

    def __init__(self, prepared: ConformerPreparation) -> None:
        super().__init__("prepared")
        self.prepared = prepared


def capture_preparation(monkeypatch: pytest.MonkeyPatch) -> list[ConformerWorkflow]:
    """Stop the executor after conformer preparation; MOPAC never runs."""
    built: list[ConformerWorkflow] = []

    class StubCalculation:
        def __init__(self, conformer_workflow: ConformerWorkflow) -> None:
            self.conformer_workflow = conformer_workflow

        def run(self, request: ConformerRequest, configuration: Any, **_: Any) -> None:
            raise _Prepared(
                self.conformer_workflow.prepare(
                    request,
                    search_settings=configuration.conformer_search,
                    selection_settings=configuration.conformer_selection,
                )
            )

    def from_pm7_config(**kwargs: Any) -> StubCalculation:
        built.append(kwargs["conformer_workflow"])
        return StubCalculation(kwargs["conformer_workflow"])

    monkeypatch.setattr(
        SemiImperiumCalculationWorkflow, "from_pm7_config", from_pm7_config
    )
    return built


def execution_request(settings: SemiImperiumSettings) -> ExecutionRequest:
    return ExecutionRequest(
        run_id="run-1",
        calculation_id="calc-1",
        identity=MolecularIdentity.from_smiles("CCO"),
        configuration=settings.configuration_for("PM7", crest_enabled=False),
        hamiltonian="PM7",
    )


def settings_with(tmp_path: Path, strategy: str) -> SemiImperiumSettings:
    base = SemiImperiumSettings()
    return replace(
        base,
        conformer_selection=ConformerSelectionSettings(strategy=strategy),
        runtime=replace(base.runtime, store_root=tmp_path),
    )


def test_executor_runs_confpass_with_the_port_and_a_smiles_topology(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    built = capture_preparation(monkeypatch)
    settings = settings_with(tmp_path, CONFPASS)

    with pytest.raises(_Prepared) as captured:
        ScientificCalculationExecutor(settings).execute(execution_request(settings))

    assert isinstance(built[0]._confpass_backend, PortedConfPass)
    assert isinstance(built[0]._topology_provider, SmilesTopology)
    selection = captured.value.prepared.selection
    assert selection.ranking_basis == "single_conformer_ensemble"
    assert selection.is_experimental


def test_executor_keeps_an_injected_backend_and_skips_the_port_by_default(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    built = capture_preparation(monkeypatch)
    injected = OrderPreservingConfPass()
    confpass_settings = settings_with(tmp_path, CONFPASS)
    default_settings = settings_with(
        tmp_path, ConformerSelectionStrategy.CREST_ENERGY_TOP_N.value
    )

    with pytest.raises(_Prepared):
        ScientificCalculationExecutor(
            confpass_settings, confpass_backend=injected
        ).execute(execution_request(confpass_settings))
    with pytest.raises(_Prepared):
        ScientificCalculationExecutor(default_settings).execute(
            execution_request(default_settings)
        )

    assert built[0]._confpass_backend is injected
    assert built[1]._confpass_backend is None
