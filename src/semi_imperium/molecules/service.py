"""Application service for traceable molecule input and resolution."""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum

from semi_imperium.domain.identity import MolecularIdentity, MoleculeInputType
from semi_imperium.molecules.resolvers import (
    MoleculeResolver,
    ResolutionCandidate,
    ResolverError,
)
from semi_imperium.molecules.validation import (
    ComputationalPath,
    MoleculeValidationError,
    MoleculeValidator,
    ValidatedStructure,
)


class ResolutionStatus(str, Enum):
    """Explicit state of one pending molecule input."""

    RESOLVED = "resolved"
    AMBIGUOUS = "ambiguous"
    UNRESOLVED = "unresolved"
    INVALID = "invalid"
    REMOVED = "removed"


@dataclass(frozen=True)
class ResolutionOutcome:
    """One input's state plus the actions needed to make it calculable."""

    original_input: str
    input_type: MoleculeInputType
    status: ResolutionStatus
    identity: MolecularIdentity | None = None
    candidates: tuple[ResolutionCandidate, ...] = ()
    message: str | None = None

    def __post_init__(self) -> None:
        if self.status is ResolutionStatus.RESOLVED and self.identity is None:
            raise ValueError("A resolved outcome must carry a MolecularIdentity")
        if self.status is not ResolutionStatus.RESOLVED and self.identity is not None:
            raise ValueError("Only a resolved outcome may carry a MolecularIdentity")
        if self.status is ResolutionStatus.AMBIGUOUS and len(self.candidates) < 2:
            raise ValueError("An ambiguous outcome must carry at least two candidates")

    @property
    def requires_disambiguation(self) -> bool:
        return self.status is ResolutionStatus.AMBIGUOUS

    @property
    def can_enter_manual_smiles(self) -> bool:
        return (
            self.input_type is MoleculeInputType.CHEMICAL_NAME
            and self.status is ResolutionStatus.UNRESOLVED
        )

    @property
    def can_remove(self) -> bool:
        return self.status not in {ResolutionStatus.RESOLVED, ResolutionStatus.REMOVED}

    @property
    def available_actions(self) -> tuple[str, ...]:
        if self.status is ResolutionStatus.AMBIGUOUS:
            return ("select_candidate", "remove")
        if self.can_enter_manual_smiles:
            return ("enter_manual_smiles", "remove")
        if self.can_remove:
            return ("remove",)
        return ()

    def require_identity(self) -> MolecularIdentity:
        if self.identity is None:
            raise MoleculeResolutionError(
                self.message or f"molecule is {self.status.value}"
            )
        return self.identity


class MoleculeResolutionError(ValueError):
    """Raised when a caller requires an identity from an unresolved outcome."""


class MoleculeResolutionService:
    """Coordinate external resolution and mandatory local validation."""

    def __init__(
        self,
        resolver: MoleculeResolver | None = None,
        validator: MoleculeValidator | None = None,
    ) -> None:
        if resolver is None:
            from semi_imperium.molecules.pubchem import PubChemResolver

            resolver = PubChemResolver()
        self._resolver = resolver
        self._validator = validator or MoleculeValidator()

    @classmethod
    def with_pubchem(cls) -> MoleculeResolutionService:
        """Build the production service without leaking PubChem into callers."""
        from semi_imperium.molecules.pubchem import PubChemResolver

        return cls(PubChemResolver())

    def resolve(
        self,
        value: str,
        input_type: MoleculeInputType | str,
        *,
        path: ComputationalPath | None = None,
        charge: int | None = None,
        multiplicity: int = 1,
    ) -> ResolutionOutcome:
        """Resolve either supported input type through one uniform entry point."""
        parsed_type = MoleculeInputType(input_type)
        if parsed_type is MoleculeInputType.CHEMICAL_NAME:
            return self.resolve_name(
                value, path=path, charge=charge, multiplicity=multiplicity
            )
        return self.resolve_smiles(
            value, path=path, charge=charge, multiplicity=multiplicity
        )

    def resolve_name(
        self,
        chemical_name: str,
        *,
        path: ComputationalPath | None = None,
        charge: int | None = None,
        multiplicity: int = 1,
        selected_identifier: str | None = None,
    ) -> ResolutionOutcome:
        """Resolve a name, refusing to choose among different structures."""
        original_input = chemical_name if isinstance(chemical_name, str) else ""
        if not original_input.strip():
            return ResolutionOutcome(
                original_input=original_input,
                input_type=MoleculeInputType.CHEMICAL_NAME,
                status=ResolutionStatus.INVALID,
                message="Chemical name is empty; enter a name or remove this item",
            )
        try:
            raw_candidates = tuple(self._resolver.resolve(original_input.strip()))
        except ResolverError as exc:
            return ResolutionOutcome(
                original_input=original_input,
                input_type=MoleculeInputType.CHEMICAL_NAME,
                status=ResolutionStatus.UNRESOLVED,
                message=str(exc),
            )

        if not raw_candidates:
            return self._unresolved_name(
                original_input,
                f"No structure was found for {original_input.strip()!r}; "
                "enter SMILES manually or remove this molecule",
            )

        if selected_identifier is not None:
            selected = next(
                (
                    candidate
                    for candidate in raw_candidates
                    if candidate.resolver_identifier == selected_identifier
                ),
                None,
            )
            if selected is not None:
                return self._validate_candidate(
                    original_input,
                    selected,
                    path=path,
                    charge=charge,
                    multiplicity=multiplicity,
                    invalid_status=ResolutionStatus.INVALID,
                )

        candidates = self._distinct_structures(raw_candidates)
        if selected_identifier is not None:
            return ResolutionOutcome(
                original_input=original_input,
                input_type=MoleculeInputType.CHEMICAL_NAME,
                status=(
                    ResolutionStatus.AMBIGUOUS
                    if len(candidates) > 1
                    else ResolutionStatus.INVALID
                ),
                candidates=candidates if len(candidates) > 1 else (),
                message=(
                    f"Candidate {selected_identifier!r} is not available; "
                    "select one of the listed resolver identifiers"
                ),
            )
        if len(candidates) > 1:
            return ResolutionOutcome(
                original_input=original_input,
                input_type=MoleculeInputType.CHEMICAL_NAME,
                status=ResolutionStatus.AMBIGUOUS,
                candidates=candidates,
                message=(
                    f"{len(candidates)} different structures match "
                    f"{original_input.strip()!r}; explicitly select one candidate"
                ),
            )
        return self._validate_candidate(
            original_input,
            candidates[0],
            path=path,
            charge=charge,
            multiplicity=multiplicity,
            invalid_status=ResolutionStatus.UNRESOLVED,
        )

    def resolve_smiles(
        self,
        smiles: str,
        *,
        path: ComputationalPath | None = None,
        charge: int | None = None,
        multiplicity: int = 1,
    ) -> ResolutionOutcome:
        """Validate a direct SMILES input without contacting a resolver."""
        original_input = smiles if isinstance(smiles, str) else ""
        try:
            validated = self._validator.validate_smiles(
                smiles,
                path=path,
                expected_charge=charge,
                multiplicity=multiplicity,
            )
        except MoleculeValidationError as exc:
            return ResolutionOutcome(
                original_input=original_input,
                input_type=MoleculeInputType.SMILES,
                status=ResolutionStatus.INVALID,
                message=str(exc),
            )
        identity = MolecularIdentity(
            canonical_smiles=validated.canonical_smiles,
            charge=validated.formal_charge if charge is None else charge,
            multiplicity=multiplicity,
            inchikey=validated.inchikey,
            original_input=original_input,
            input_type=MoleculeInputType.SMILES,
            inchi=validated.inchi,
            resolver="manual",
        )
        return ResolutionOutcome(
            original_input=original_input,
            input_type=MoleculeInputType.SMILES,
            status=ResolutionStatus.RESOLVED,
            identity=identity,
        )

    def select_candidate(
        self,
        outcome: ResolutionOutcome,
        resolver_identifier: str,
        *,
        path: ComputationalPath | None = None,
        charge: int | None = None,
        multiplicity: int = 1,
    ) -> ResolutionOutcome:
        """Validate an explicit choice from an ambiguous outcome."""
        if not outcome.requires_disambiguation:
            raise MoleculeResolutionError(
                "Candidate selection requires an ambiguous resolution outcome"
            )
        selected = next(
            (
                candidate
                for candidate in outcome.candidates
                if candidate.resolver_identifier == resolver_identifier
            ),
            None,
        )
        if selected is None:
            raise MoleculeResolutionError(
                f"Candidate {resolver_identifier!r} is not present in this outcome"
            )
        return self._validate_candidate(
            outcome.original_input,
            selected,
            path=path,
            charge=charge,
            multiplicity=multiplicity,
            invalid_status=ResolutionStatus.INVALID,
        )

    def enter_manual_smiles(
        self,
        outcome: ResolutionOutcome,
        smiles: str,
        *,
        path: ComputationalPath | None = None,
        charge: int | None = None,
        multiplicity: int = 1,
    ) -> ResolutionOutcome:
        """Recover an unresolved name with a user-supplied explicit structure."""
        if not outcome.can_enter_manual_smiles:
            raise MoleculeResolutionError(
                "Manual SMILES entry is only available for an unresolved chemical name"
            )
        try:
            validated = self._validator.validate_smiles(
                smiles,
                path=path,
                expected_charge=charge,
                multiplicity=multiplicity,
            )
        except MoleculeValidationError as exc:
            return replace(outcome, message=f"Manual SMILES is invalid: {exc}")
        identity = MolecularIdentity(
            canonical_smiles=validated.canonical_smiles,
            charge=validated.formal_charge if charge is None else charge,
            multiplicity=multiplicity,
            inchikey=validated.inchikey,
            original_input=outcome.original_input,
            input_type=MoleculeInputType.CHEMICAL_NAME,
            display_name=outcome.original_input.strip(),
            inchi=validated.inchi,
            resolver="manual",
        )
        return ResolutionOutcome(
            original_input=outcome.original_input,
            input_type=MoleculeInputType.CHEMICAL_NAME,
            status=ResolutionStatus.RESOLVED,
            identity=identity,
        )

    @staticmethod
    def remove(outcome: ResolutionOutcome) -> ResolutionOutcome:
        """Explicitly remove an unresolved/invalid item from the pending set."""
        if not outcome.can_remove:
            raise MoleculeResolutionError(
                f"A molecule in state {outcome.status.value!r} cannot be removed here"
            )
        return ResolutionOutcome(
            original_input=outcome.original_input,
            input_type=outcome.input_type,
            status=ResolutionStatus.REMOVED,
            message="Molecule removed from the pending calculation",
        )

    def _distinct_structures(
        self, candidates: tuple[ResolutionCandidate, ...]
    ) -> tuple[ResolutionCandidate, ...]:
        """Group equivalent SMILES before any path-specific rejection."""
        distinct: dict[str, ResolutionCandidate] = {}
        for candidate in candidates:
            try:
                canonical = self._validator.canonicalize_smiles(
                    candidate.canonical_smiles
                )
                normalized = replace(candidate, canonical_smiles=canonical)
                key = f"structure:{canonical}"
            except MoleculeValidationError:
                normalized = candidate
                key = f"unparseable:{candidate.canonical_smiles.strip()}"
            distinct.setdefault(key, normalized)
        return tuple(distinct.values())

    def _validate_candidate(
        self,
        original_input: str,
        candidate: ResolutionCandidate,
        *,
        path: ComputationalPath | None,
        charge: int | None,
        multiplicity: int,
        invalid_status: ResolutionStatus,
    ) -> ResolutionOutcome:
        try:
            validated = self._validator.validate_smiles(
                candidate.canonical_smiles,
                path=path,
                expected_charge=charge,
                multiplicity=multiplicity,
            )
        except MoleculeValidationError as exc:
            message = (
                "The selected structure is not supported by the requested "
                f"computational path: {exc}."
            )
            if invalid_status is ResolutionStatus.UNRESOLVED:
                message += " Enter SMILES manually or remove it."
                return self._unresolved_name(original_input, message)
            return ResolutionOutcome(
                original_input=original_input,
                input_type=MoleculeInputType.CHEMICAL_NAME,
                status=invalid_status,
                message=message,
            )
        return self._resolved_candidate(
            original_input,
            candidate,
            validated,
            charge=charge,
            multiplicity=multiplicity,
        )

    @staticmethod
    def _unresolved_name(original_input: str, message: str) -> ResolutionOutcome:
        return ResolutionOutcome(
            original_input=original_input,
            input_type=MoleculeInputType.CHEMICAL_NAME,
            status=ResolutionStatus.UNRESOLVED,
            message=message,
        )

    @staticmethod
    def _resolved_candidate(
        original_input: str,
        candidate: ResolutionCandidate,
        validated: ValidatedStructure,
        *,
        charge: int | None,
        multiplicity: int,
    ) -> ResolutionOutcome:
        identity = MolecularIdentity(
            canonical_smiles=validated.canonical_smiles,
            charge=validated.formal_charge if charge is None else charge,
            multiplicity=multiplicity,
            inchikey=candidate.inchikey or validated.inchikey,
            display_name=candidate.resolved_name,
            original_input=original_input,
            input_type=MoleculeInputType.CHEMICAL_NAME,
            resolved_name=candidate.resolved_name,
            inchi=candidate.inchi or validated.inchi,
            resolver=candidate.resolver,
            resolver_identifier=candidate.resolver_identifier,
            cid=candidate.cid,
        )
        return ResolutionOutcome(
            original_input=original_input,
            input_type=MoleculeInputType.CHEMICAL_NAME,
            status=ResolutionStatus.RESOLVED,
            identity=identity,
        )


__all__ = [
    "MoleculeResolutionError",
    "MoleculeResolutionService",
    "ResolutionOutcome",
    "ResolutionStatus",
]
