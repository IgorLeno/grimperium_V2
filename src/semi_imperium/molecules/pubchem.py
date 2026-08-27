"""PubChem PUG REST adapter for chemical-name resolution."""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

from semi_imperium.molecules.resolvers import (
    ResolutionCandidate,
    ResolverUnavailableError,
)

PUBCHEM_RESOLVER_ID = "pubchem"
PUBCHEM_PUG_REST_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_PROPERTIES = (
    "Title",
    "IUPACName",
    "SMILES",
    "ConnectivitySMILES",
    "InChI",
    "InChIKey",
)

JsonRequest = Callable[[str, float], Mapping[str, Any]]


class PubChemResolver:
    """Resolve names with PubChem while exposing only neutral candidates."""

    def __init__(
        self,
        *,
        request_json: JsonRequest | None = None,
        timeout_seconds: float = 10.0,
        base_url: str = PUBCHEM_PUG_REST_BASE_URL,
    ) -> None:
        if timeout_seconds <= 0:
            raise ValueError("PubChem timeout_seconds must be positive")
        self._request_json = request_json or _request_json
        self._timeout_seconds = timeout_seconds
        self._base_url = base_url.rstrip("/")

    @property
    def resolver_id(self) -> str:
        """Stable provenance identifier for this adapter."""
        return PUBCHEM_RESOLVER_ID

    def resolve(self, chemical_name: str) -> list[ResolutionCandidate]:
        """Return every PubChem structure candidate for ``chemical_name``.

        A missing PubChem record is represented by an empty list. Transport,
        decoding and schema failures remain distinguishable resolver failures.
        """
        name = chemical_name.strip()
        if not name:
            return []
        properties = ",".join(PUBCHEM_PROPERTIES)
        url = (
            f"{self._base_url}/compound/name/{quote(name, safe='')}"
            f"/property/{properties}/JSON"
        )
        try:
            payload = self._request_json(url, self._timeout_seconds)
        except HTTPError as exc:
            if exc.code in {400, 404}:
                return []
            raise ResolverUnavailableError(
                f"PubChem request failed with HTTP {exc.code}; "
                "enter SMILES manually or retry later"
            ) from exc
        except (OSError, TimeoutError, URLError, ValueError) as exc:
            raise ResolverUnavailableError(
                "PubChem is unavailable or returned invalid data; "
                "enter SMILES manually or retry later"
            ) from exc

        try:
            raw_properties = payload["PropertyTable"]["Properties"]
        except (KeyError, TypeError) as exc:
            fault = payload.get("Fault") if isinstance(payload, Mapping) else None
            if fault:
                return []
            raise ResolverUnavailableError(
                "PubChem response did not contain a property table; "
                "enter SMILES manually or retry later"
            ) from exc
        if not isinstance(raw_properties, list):
            raise ResolverUnavailableError(
                "PubChem property table has an unexpected format; "
                "enter SMILES manually or retry later"
            )

        candidates: list[ResolutionCandidate] = []
        for raw in raw_properties:
            if not isinstance(raw, Mapping):
                continue
            candidate = _candidate_from_property(raw)
            if candidate is not None:
                candidates.append(candidate)
        return candidates

    def resolve_name(self, chemical_name: str) -> list[ResolutionCandidate]:
        """Explicit alias for callers that prefer the longer method name."""
        return self.resolve(chemical_name)


def _candidate_from_property(raw: Mapping[str, Any]) -> ResolutionCandidate | None:
    cid_value = raw.get("CID")
    try:
        cid = int(str(cid_value))
    except (TypeError, ValueError):
        return None
    smiles = _first_string(
        raw,
        "SMILES",
        "IsomericSMILES",
        "CanonicalSMILES",
        "ConnectivitySMILES",
    )
    if smiles is None:
        return None
    return ResolutionCandidate(
        canonical_smiles=smiles,
        resolved_name=_first_string(raw, "Title", "IUPACName"),
        resolver=PUBCHEM_RESOLVER_ID,
        resolver_identifier=str(cid),
        inchi=_first_string(raw, "InChI"),
        inchikey=_first_string(raw, "InChIKey"),
        cid=cid,
    )


def _first_string(raw: Mapping[str, Any], *keys: str) -> str | None:
    for key in keys:
        value = raw.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return None


def _request_json(url: str, timeout_seconds: float) -> Mapping[str, Any]:
    """Fetch one PUG REST JSON object using the Python standard library."""
    request = Request(  # noqa: S310 - fixed HTTPS base URL by default
        url,
        headers={"Accept": "application/json", "User-Agent": "Semi-Imperium/0.2"},
    )
    with urlopen(request, timeout=timeout_seconds) as response:  # noqa: S310
        payload = json.load(response)
    if not isinstance(payload, Mapping):
        raise ValueError("PubChem response root is not a JSON object")
    return payload


__all__ = [
    "PUBCHEM_PUG_REST_BASE_URL",
    "PUBCHEM_RESOLVER_ID",
    "PubChemResolver",
]
