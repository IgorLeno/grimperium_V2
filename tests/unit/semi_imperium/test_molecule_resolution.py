"""Resolver, disambiguation and early chemistry-validation contract."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import pytest

from semi_imperium.domain import MolecularIdentity, MoleculeInputType
from semi_imperium.molecules import (
    ComputationalPath,
    MoleculeResolutionError,
    MoleculeResolutionService,
    MoleculeValidationError,
    MoleculeValidator,
    PubChemResolver,
    ResolutionCandidate,
    ResolutionStatus,
    ResolverUnavailableError,
)


@dataclass
class FakeResolver:
    candidates: Sequence[ResolutionCandidate] = ()
    error: Exception | None = None
    calls: int = 0

    @property
    def resolver_id(self) -> str:
        return "fake"

    def resolve(self, chemical_name: str) -> Sequence[ResolutionCandidate]:
        self.calls += 1
        if self.error is not None:
            raise self.error
        return self.candidates


def candidate(
    cid: int,
    smiles: str,
    name: str,
    *,
    inchi: str | None = None,
    inchikey: str | None = None,
) -> ResolutionCandidate:
    return ResolutionCandidate(
        canonical_smiles=smiles,
        resolved_name=name,
        resolver="pubchem",
        resolver_identifier=str(cid),
        cid=cid,
        inchi=inchi,
        inchikey=inchikey,
    )


def test_pubchem_adapter_maps_pug_rest_without_live_network() -> None:
    requested: list[tuple[str, float]] = []

    def fake_request(url: str, timeout: float) -> Mapping[str, Any]:
        requested.append((url, timeout))
        return {
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 702,
                        "Title": "Ethanol",
                        "SMILES": "CCO",
                        "ConnectivitySMILES": "CCO",
                        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                    }
                ]
            }
        }

    resolver = PubChemResolver(request_json=fake_request, timeout_seconds=3.5)

    assert resolver.resolve("ethyl alcohol") == [
        candidate(
            702,
            "CCO",
            "Ethanol",
            inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        )
    ]
    assert requested == [
        (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            "ethyl%20alcohol/property/Title,IUPACName,SMILES,ConnectivitySMILES,"
            "InChI,InChIKey/JSON",
            3.5,
        )
    ]


def test_resolved_name_preserves_identity_and_resolution_provenance() -> None:
    resolver = FakeResolver(
        [
            candidate(
                702,
                "OCC",
                "Ethanol",
                inchi="pubchem-inchi",
                inchikey="PUBCHEM-INCHIKEY",
            )
        ]
    )

    outcome = MoleculeResolutionService(resolver).resolve_name("ethyl alcohol")
    identity = outcome.require_identity()

    assert outcome.status is ResolutionStatus.RESOLVED
    assert identity.original_input == "ethyl alcohol"
    assert identity.input_type is MoleculeInputType.CHEMICAL_NAME
    assert identity.resolved_name == "Ethanol"
    assert identity.canonical_smiles == "CCO"
    assert identity.inchi == "pubchem-inchi"
    assert identity.inchikey == "PUBCHEM-INCHIKEY"
    assert identity.resolver == "pubchem"
    assert identity.resolver_identifier == "702"
    assert identity.cid == 702
    assert MolecularIdentity.from_dict(identity.to_dict()) == identity


def test_different_structures_require_explicit_disambiguation() -> None:
    service = MoleculeResolutionService(
        FakeResolver(
            [
                candidate(702, "CCO", "ethanol"),
                candidate(8071, "COC", "dimethyl ether"),
            ]
        )
    )

    outcome = service.resolve_name("C2H6O")

    assert outcome.status is ResolutionStatus.AMBIGUOUS
    assert outcome.identity is None
    assert outcome.requires_disambiguation is True
    assert outcome.available_actions == ("select_candidate", "remove")
    assert {item.resolver_identifier for item in outcome.candidates} == {"702", "8071"}

    selected = service.select_candidate(outcome, "8071")
    assert selected.require_identity().canonical_smiles == "COC"
    assert selected.require_identity().cid == 8071


def test_path_failure_does_not_silently_eliminate_an_ambiguous_candidate() -> None:
    service = MoleculeResolutionService(
        FakeResolver(
            [
                candidate(702, "CCO", "ethanol"),
                candidate(23990, "[U]", "uranium"),
            ]
        )
    )

    outcome = service.resolve_name("ambiguous name")

    assert outcome.status is ResolutionStatus.AMBIGUOUS
    assert {item.resolver_identifier for item in outcome.candidates} == {
        "702",
        "23990",
    }
    unsupported = service.select_candidate(outcome, "23990")
    assert unsupported.status is ResolutionStatus.INVALID
    assert "UFF lacks parameters" in (unsupported.message or "")


def test_equivalent_candidate_spellings_are_not_chemical_ambiguity() -> None:
    service = MoleculeResolutionService(
        FakeResolver(
            [
                candidate(1, "CCO", "ethanol"),
                candidate(2, "OCC", "ethyl alcohol"),
            ]
        )
    )

    outcome = service.resolve_name("ethanol")

    assert outcome.status is ResolutionStatus.RESOLVED
    assert outcome.require_identity().canonical_smiles == "CCO"


def test_unresolved_name_offers_manual_smiles_or_removal() -> None:
    service = MoleculeResolutionService(FakeResolver())
    unresolved = service.resolve_name("not in resolver")

    assert unresolved.status is ResolutionStatus.UNRESOLVED
    assert unresolved.can_enter_manual_smiles is True
    assert unresolved.available_actions == ("enter_manual_smiles", "remove")

    manual = service.enter_manual_smiles(unresolved, "OCC")
    identity = manual.require_identity()
    assert identity.original_input == "not in resolver"
    assert identity.input_type is MoleculeInputType.CHEMICAL_NAME
    assert identity.canonical_smiles == "CCO"
    assert identity.resolver == "manual"

    removed = service.remove(unresolved)
    assert removed.status is ResolutionStatus.REMOVED
    assert removed.identity is None


def test_resolver_outage_remains_recoverable_without_claiming_no_match() -> None:
    service = MoleculeResolutionService(
        FakeResolver(error=ResolverUnavailableError("PubChem timed out; retry later"))
    )

    outcome = service.resolve_name("ethanol")

    assert outcome.status is ResolutionStatus.UNRESOLVED
    assert "timed out" in (outcome.message or "")
    assert outcome.can_enter_manual_smiles is True


def test_direct_smiles_never_contacts_name_resolver_and_is_canonicalized() -> None:
    resolver = FakeResolver()
    outcome = MoleculeResolutionService(resolver).resolve(
        " OCC ", MoleculeInputType.SMILES
    )

    identity = outcome.require_identity()
    assert resolver.calls == 0
    assert identity.original_input == " OCC "
    assert identity.input_type is MoleculeInputType.SMILES
    assert identity.canonical_smiles == "CCO"
    assert identity.inchi is not None
    assert identity.inchikey is not None


@pytest.mark.parametrize(
    ("smiles", "message"),
    [
        ("not-a-smiles", "RDKit could not parse"),
        ("CC.O", "Disconnected fragments"),
        ("*CC", "Wildcard/query atoms"),
    ],
)
def test_malformed_or_unsupported_smiles_fails_before_calculation(
    smiles: str, message: str
) -> None:
    outcome = MoleculeResolutionService(FakeResolver()).resolve_smiles(smiles)

    assert outcome.status is ResolutionStatus.INVALID
    assert message in (outcome.message or "")
    with pytest.raises(MoleculeResolutionError):
        outcome.require_identity()


def test_requested_charge_must_match_the_explicit_structure() -> None:
    outcome = MoleculeResolutionService(FakeResolver()).resolve_smiles(
        "[NH4+]", charge=0
    )

    assert outcome.status is ResolutionStatus.INVALID
    assert "does not match" in (outcome.message or "")


def test_multiplicity_must_match_electron_count_parity() -> None:
    outcome = MoleculeResolutionService(FakeResolver()).resolve_smiles(
        "[CH3]", multiplicity=1
    )

    assert outcome.status is ResolutionStatus.INVALID
    assert "electron count" in (outcome.message or "")

    doublet = MoleculeResolutionService(FakeResolver()).resolve_smiles(
        "[CH3]", multiplicity=2
    )
    assert doublet.status is ResolutionStatus.RESOLVED


def test_computational_path_rejects_unknown_hamiltonian_immediately() -> None:
    with pytest.raises(MoleculeValidationError, match="Choose AM1, PM3 or PM7"):
        ComputationalPath(hamiltonians=("PM6",))


def test_local_validator_rejects_atoms_without_runner_3d_parameters() -> None:
    with pytest.raises(MoleculeValidationError, match="UFF lacks parameters"):
        MoleculeValidator().validate_smiles("[U]")
