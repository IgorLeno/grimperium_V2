"""Filesystem persistence for Semi-Imperium runs and calculations.

The store owns a directory of its own and writes nothing outside it. In
particular it never touches Grimperium's datasets (``data/*.csv``), never
rewrites Grimperium's batch CSVs and never reinterprets values that were
produced by another tool: it records what Semi-Imperium itself computed,
together with enough identity and provenance to explain each number.

Layout::

    <root>/runs/<run_id>/run.json
    <root>/calculations/<molecule_id>/<signature_digest>/<calculation_id>.json

Calculations are stored under their reuse key, so answering "have we
already computed this molecule under this configuration?" is a directory
listing rather than a scan. Listing a single run's calculations is the
rarer question and does walk the tree.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from collections.abc import Collection, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from semi_imperium.domain.configuration import CalculationSignature
from semi_imperium.domain.enums import DEFAULT_REUSABLE_STATES, CalculationState
from semi_imperium.domain.identity import MolecularIdentity
from semi_imperium.domain.records import (
    CalculationRecord,
    RunRecord,
    build_calculation_id,
)

RUNS_DIRNAME = "runs"
CALCULATIONS_DIRNAME = "calculations"
RUN_FILENAME = "run.json"

_SAFE_SEGMENT = re.compile(r"^[A-Za-z0-9._-]+$")


class StoreIntegrityError(RuntimeError):
    """Raised when persisted state would be corrupted or contradicted."""


@dataclass(frozen=True)
class ReuseDecision:
    """Why a requested calculation was or was not served from the store."""

    reuse_key: str
    record: CalculationRecord | None
    reason: str

    @property
    def can_reuse(self) -> bool:
        """Whether a stored record satisfied the request."""
        return self.record is not None


class SemiImperiumStore:
    """Durable, append-oriented store for runs and individual calculations."""

    def __init__(self, root: Path | str) -> None:
        self._root = Path(root)

    @property
    def root(self) -> Path:
        """Root directory owned by this store."""
        return self._root

    @property
    def runs_dir(self) -> Path:
        """Directory holding one subdirectory per run."""
        return self._root / RUNS_DIRNAME

    @property
    def calculations_dir(self) -> Path:
        """Directory holding calculations partitioned by reuse key."""
        return self._root / CALCULATIONS_DIRNAME

    # ------------------------------------------------------------------
    # Runs
    # ------------------------------------------------------------------

    def run_path(self, run_id: str) -> Path:
        """Return the manifest path for ``run_id``."""
        return self.runs_dir / _safe_segment(run_id, "run_id") / RUN_FILENAME

    def save_run(self, run: RunRecord) -> Path:
        """Persist a run manifest atomically and return its path.

        Raises:
            StoreIntegrityError: If a terminal run already on disk would be
                moved to a different state by this write.
        """
        path = self.run_path(run.run_id)
        existing = self._read_run(path)
        if (
            existing is not None
            and existing.state.is_terminal
            and existing.state is not run.state
        ):
            raise StoreIntegrityError(
                f"Run {run.run_id!r} is already terminal in state "
                f"{existing.state.value!r}; refusing to overwrite with "
                f"{run.state.value!r}"
            )
        _atomic_write_json(path, run.to_dict())
        return path

    def load_run(self, run_id: str) -> RunRecord:
        """Load one run manifest.

        Raises:
            FileNotFoundError: If no manifest exists for ``run_id``.
        """
        path = self.run_path(run_id)
        run = self._read_run(path)
        if run is None:
            raise FileNotFoundError(f"No Semi-Imperium run manifest at {path}")
        return run

    def list_runs(self) -> list[RunRecord]:
        """Return every stored run, most recently created first."""
        if not self.runs_dir.is_dir():
            return []
        runs = [
            run
            for path in sorted(self.runs_dir.glob(f"*/{RUN_FILENAME}"))
            if (run := self._read_run(path)) is not None
        ]
        return sorted(runs, key=lambda item: item.timestamps.created_at, reverse=True)

    # ------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------

    def reuse_dir(
        self, molecule: MolecularIdentity, signature: CalculationSignature
    ) -> Path:
        """Return the directory that holds every attempt for one reuse key."""
        return (
            self.calculations_dir
            / _safe_segment(molecule.molecule_id, "molecule_id")
            / _safe_segment(signature.digest, "signature digest")
        )

    def calculation_path(self, record: CalculationRecord) -> Path:
        """Return the file path for one calculation record."""
        directory = self.reuse_dir(record.molecule, record.signature)
        return directory / f"{_safe_segment(record.calculation_id, 'id')}.json"

    def save_calculation(self, record: CalculationRecord) -> Path:
        """Persist one calculation atomically and return its path.

        Re-saving the same terminal record is idempotent; moving an
        already-terminal record to a different state is refused, because
        a finished scientific result must not be silently rewritten.

        Raises:
            StoreIntegrityError: If the write would contradict stored state.
        """
        path = self.calculation_path(record)
        existing = self._read_calculation(path)
        if existing is not None and existing.state.is_terminal:
            if existing.state is not record.state:
                raise StoreIntegrityError(
                    f"Calculation {record.calculation_id!r} is already terminal "
                    f"in state {existing.state.value!r}; refusing to overwrite "
                    f"with {record.state.value!r}"
                )
            if existing.verification is not record.verification:
                raise StoreIntegrityError(
                    f"Calculation {record.calculation_id!r} is already terminal "
                    f"with verification {existing.verification.value!r}; refusing "
                    f"to overwrite with {record.verification.value!r}"
                )
        _atomic_write_json(path, record.to_dict())
        return path

    def get_calculation(
        self,
        *,
        run_id: str,
        molecule: MolecularIdentity,
        signature: CalculationSignature,
    ) -> CalculationRecord | None:
        """Return the record for one molecule in one run, if it exists."""
        calculation_id = build_calculation_id(
            run_id=run_id,
            molecule_id=molecule.molecule_id,
            signature_digest=signature.digest,
        )
        directory = self.reuse_dir(molecule, signature)
        return self._read_calculation(directory / f"{calculation_id}.json")

    def list_calculations(
        self,
        *,
        molecule: MolecularIdentity | None = None,
        signature: CalculationSignature | None = None,
        run_id: str | None = None,
    ) -> list[CalculationRecord]:
        """Return stored calculations, most recently created first.

        All filters are optional and combine with AND.
        """
        records = [
            record
            for record in self._iter_calculations(molecule, signature)
            if run_id is None or record.run_id == run_id
        ]
        return sorted(
            records, key=lambda item: item.timestamps.created_at, reverse=True
        )

    def find_reusable(
        self,
        molecule: MolecularIdentity,
        signature: CalculationSignature,
        *,
        accepted_states: Collection[CalculationState] = DEFAULT_REUSABLE_STATES,
    ) -> ReuseDecision:
        """Decide whether a stored result can stand in for a new calculation.

        A record qualifies only when it was produced for the *same*
        molecular identity under the *same* configuration signature, and
        its state is one the caller accepts. Callers that require a proven
        minimum pass ``accepted_states={CalculationState.VERIFIED}``.
        """
        reuse_key = f"{molecule.molecule_id}/{signature.digest}"
        candidates = self._iter_calculations(molecule, signature)
        accepted = [record for record in candidates if record.state in accepted_states]
        if not accepted:
            return ReuseDecision(
                reuse_key=reuse_key,
                record=None,
                reason="no stored calculation matches this identity and signature",
            )

        # Most recently completed wins: it reflects the newest software
        # versions recorded in provenance for an otherwise identical setup.
        best = max(
            accepted,
            key=lambda item: (
                item.timestamps.completed_at or item.timestamps.created_at
            ),
        )
        return ReuseDecision(
            reuse_key=reuse_key,
            record=best,
            reason=f"reusing {best.calculation_id} in state {best.state.value}",
        )

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _iter_calculations(
        self,
        molecule: MolecularIdentity | None,
        signature: CalculationSignature | None,
    ) -> Iterator[CalculationRecord]:
        """Yield stored calculations under the narrowest matching subtree."""
        if molecule is not None and signature is not None:
            roots = [self.reuse_dir(molecule, signature)]
        elif molecule is not None:
            molecule_dir = self.calculations_dir / _safe_segment(
                molecule.molecule_id, "molecule_id"
            )
            roots = sorted(p for p in molecule_dir.glob("*") if p.is_dir())
        else:
            roots = [self.calculations_dir]

        for directory in roots:
            if not directory.is_dir():
                continue
            for path in sorted(directory.rglob("*.json")):
                record = self._read_calculation(path)
                if record is None:
                    continue
                if (
                    signature is not None
                    and record.signature.digest != signature.digest
                ):
                    continue
                yield record

    def _read_run(self, path: Path) -> RunRecord | None:
        payload = _read_json(path)
        return None if payload is None else RunRecord.from_dict(payload)

    def _read_calculation(self, path: Path) -> CalculationRecord | None:
        payload = _read_json(path)
        return None if payload is None else CalculationRecord.from_dict(payload)


def _read_json(path: Path) -> dict[str, Any] | None:
    """Read a JSON object, returning ``None`` when the file is absent."""
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise StoreIntegrityError(f"Stored record is not a JSON object: {path}")
    return payload


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write ``payload`` durably: temp file, fsync, then atomic rename."""
    path.parent.mkdir(parents=True, exist_ok=True)
    content = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    temp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.stem}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temp_path = Path(handle.name)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, path)
        temp_path = None
    finally:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink()


def _safe_segment(value: str, label: str) -> str:
    """Reject identifiers that would escape or confuse the store layout."""
    if value in {".", ".."} or not _SAFE_SEGMENT.match(value):
        raise ValueError(
            f"Unsafe {label} for a store path: {value!r}; only letters, digits, "
            "'.', '_' and '-' are allowed"
        )
    return value


__all__ = [
    "CALCULATIONS_DIRNAME",
    "RUNS_DIRNAME",
    "ReuseDecision",
    "SemiImperiumStore",
    "StoreIntegrityError",
]
