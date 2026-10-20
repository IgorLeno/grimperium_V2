"""Results waiting for the central, on disk, resent until it answers.

Every finished attempt is written here *before* it is sent, one file per
``attempt_id`` (atomic write: a power cut leaves the old file or the new one,
never half of it). A file leaves the outbox only when the central settled it
(any 2xx, ``duplicate`` included: the central's acceptance is idempotent by
``attempt_id``) or refused it for good (4xx other than 401), in which case it
moves to the dead letter with the refusal so nothing is lost silently.
"""

from __future__ import annotations

import json
import os
import random
import re
import tempfile
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

OUTBOX_DIRNAME = "outbox"
DEAD_LETTER_DIRNAME = "dead-letter"

# attempt_id comes from the network and names a file: a plain name only.
_SAFE_ID = re.compile(r"[A-Za-z0-9_-]{1,64}")


class Verdict(Enum):
    SENT = "sent"  # settled by the central: remove
    RETRY = "retry"  # no answer or a passing error: keep, wait, resend
    DEAD = "dead"  # refused for good: dead letter
    STOP = "stop"  # token refused: keep, and stop talking to the central


@dataclass(frozen=True)
class Reply:
    """What one request got: ``status`` is ``None`` when nothing answered."""

    status: int | None
    body: dict[str, Any] = field(default_factory=dict)
    error: str = ""

    @property
    def detail(self) -> str:
        text = self.body.get("detail")
        return str(text) if text is not None else self.error


def classify(status: int | None) -> Verdict:
    """Whether an answer settles a result, retries it, kills it or stops."""
    if status is None or status >= 500 or status in (408, 429):
        return Verdict.RETRY
    if 200 <= status < 300:
        return Verdict.SENT
    if status == 401:
        return Verdict.STOP
    return Verdict.DEAD


@dataclass(frozen=True)
class Entry:
    attempt_id: str
    mol_id: str
    request: dict[str, Any]
    queued_ns: int


@dataclass(frozen=True)
class DeadLetter:
    attempt_id: str
    mol_id: str
    request: dict[str, Any]
    status: int | None
    detail: str


@dataclass(frozen=True)
class FlushReport:
    sent: int = 0
    dead: int = 0
    retry: bool = False
    stop: bool = False
    last_error: str = ""


def _atomic_write(path: Path, payload: dict[str, Any]) -> None:
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise
    # The rename itself must survive a power cut too.
    directory = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


class Outbox:
    """``<root>/outbox/<attempt_id>.json`` and ``<root>/dead-letter/``."""

    def __init__(self, root: Path) -> None:
        self.pending_dir = Path(root) / OUTBOX_DIRNAME
        self.dead_dir = Path(root) / DEAD_LETTER_DIRNAME
        self.pending_dir.mkdir(parents=True, exist_ok=True)
        self.dead_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _name(attempt_id: str) -> str:
        if not _SAFE_ID.fullmatch(attempt_id):
            raise ValueError(f"attempt_id inválido: {attempt_id!r}")
        return f"{attempt_id}.json"

    def put(self, attempt_id: str, mol_id: str, request: dict[str, Any]) -> None:
        """Keep ``request`` (the ``/api/cell/result`` body) until it is settled."""
        _atomic_write(
            self.pending_dir / self._name(attempt_id),
            {
                "attempt_id": attempt_id,
                "mol_id": mol_id,
                "queued_ns": time.time_ns(),
                "request": request,
            },
        )

    def pending(self) -> list[Entry]:
        entries = []
        for path in self.pending_dir.glob("*.json"):
            data = json.loads(path.read_text(encoding="utf-8"))
            entries.append(
                Entry(
                    attempt_id=str(data["attempt_id"]),
                    mol_id=str(data["mol_id"]),
                    request=dict(data["request"]),
                    queued_ns=int(data["queued_ns"]),
                )
            )
        return sorted(entries, key=lambda e: (e.queued_ns, e.attempt_id))

    def remove(self, attempt_id: str) -> None:
        (self.pending_dir / self._name(attempt_id)).unlink(missing_ok=True)

    def kill(self, entry: Entry, reply: Reply) -> None:
        """Move ``entry`` to the dead letter with the central's refusal."""
        _atomic_write(
            self.dead_dir / self._name(entry.attempt_id),
            {
                "attempt_id": entry.attempt_id,
                "mol_id": entry.mol_id,
                "request": entry.request,
                "status": reply.status,
                "detail": reply.detail,
                "dead_at_ns": time.time_ns(),
            },
        )
        self.remove(entry.attempt_id)

    def dead_letters(self) -> list[DeadLetter]:
        letters = []
        for path in sorted(self.dead_dir.glob("*.json")):
            data = json.loads(path.read_text(encoding="utf-8"))
            letters.append(
                DeadLetter(
                    attempt_id=str(data["attempt_id"]),
                    mol_id=str(data["mol_id"]),
                    request=dict(data["request"]),
                    status=data["status"],
                    detail=str(data["detail"]),
                )
            )
        return letters

    def counts(self) -> tuple[int, int]:
        """(waiting to be sent, in the dead letter)."""
        return (
            sum(1 for _ in self.pending_dir.glob("*.json")),
            sum(1 for _ in self.dead_dir.glob("*.json")),
        )


Send = Callable[[dict[str, Any]], Reply]


def flush(outbox: Outbox, send: Send) -> FlushReport:
    """Send every pending result, oldest first.

    Stops at the first result that must be retried (the central is down or
    failing: the rest would fail too) or at a refused token.
    """
    sent = dead = 0
    for entry in outbox.pending():
        reply = send(entry.request)
        verdict = classify(reply.status)
        if verdict is Verdict.SENT:
            outbox.remove(entry.attempt_id)
            sent += 1
        elif verdict is Verdict.DEAD:
            outbox.kill(entry, reply)
            dead += 1
        else:
            return FlushReport(
                sent=sent,
                dead=dead,
                retry=verdict is Verdict.RETRY,
                stop=verdict is Verdict.STOP,
                last_error=f"{reply.status or 'sem resposta'}: {reply.detail}",
            )
    return FlushReport(sent=sent, dead=dead)


@dataclass
class Backoff:
    """Growing wait between retries: ``initial_s`` doubling up to ``cap_s``.

    Each delay is drawn from the upper half of its ceiling, so cells cut off
    together by the same outage do not all come back in the same second.
    """

    initial_s: float = 5.0
    cap_s: float = 300.0
    rng: random.Random = field(default_factory=random.Random)

    def delay(self, failures: int) -> float:
        ceiling = min(self.cap_s, self.initial_s * 2.0 ** max(failures - 1, 0))
        return ceiling * (0.5 + 0.5 * self.rng.random())
