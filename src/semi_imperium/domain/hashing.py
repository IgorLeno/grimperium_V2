"""Deterministic hashing shared by Semi-Imperium identities and signatures.

Every identifier in this package is a SHA-256 digest over a JSON payload
encoded with sorted keys and no insignificant whitespace. Keeping that in
one place is what makes ids reproducible across processes, machines and
Python versions.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any

DIGEST_ALGORITHM = "sha256"


def canonical_json(payload: Any) -> str:
    """Return the canonical JSON encoding used for every digest."""
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def stable_digest(payload: Any) -> str:
    """Return the SHA-256 hex digest of the canonical encoding of ``payload``."""
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


__all__ = ["DIGEST_ALGORITHM", "canonical_json", "stable_digest"]
