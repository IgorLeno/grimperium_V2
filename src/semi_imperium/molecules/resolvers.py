"""Resolver-independent contracts for turning names into structures."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol, runtime_checkable


class ResolverError(RuntimeError):
    """Base error raised by an external molecule resolver adapter."""


class ResolverUnavailableError(ResolverError):
    """Raised when a resolver cannot be contacted or returns unusable data."""


@dataclass(frozen=True)
class ResolutionCandidate:
    """Resolver-neutral candidate structure and its external provenance."""

    canonical_smiles: str
    resolved_name: str | None
    resolver: str
    resolver_identifier: str
    inchi: str | None = None
    inchikey: str | None = None
    cid: int | None = None

    def __post_init__(self) -> None:
        for field_name in ("canonical_smiles", "resolver", "resolver_identifier"):
            value = getattr(self, field_name)
            if not value or not value.strip():
                raise ValueError(f"ResolutionCandidate.{field_name} must not be empty")
        if self.cid is not None and self.cid <= 0:
            raise ValueError(
                f"ResolutionCandidate.cid must be positive, got {self.cid}"
            )

    @property
    def identifier(self) -> str:
        """Short alias used by disambiguation callers."""
        return self.resolver_identifier


@runtime_checkable
class MoleculeResolver(Protocol):
    """Boundary implemented by chemical-name resolution providers."""

    @property
    def resolver_id(self) -> str:
        """Stable identifier recorded as scientific provenance."""

    def resolve(self, chemical_name: str) -> Sequence[ResolutionCandidate]:
        """Return zero or more candidates without silently selecting one."""


__all__ = [
    "MoleculeResolver",
    "ResolutionCandidate",
    "ResolverError",
    "ResolverUnavailableError",
]
