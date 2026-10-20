"""The cell's offline queue of results: on disk, resent until the central answers."""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from semi_imperium.cell.outbox import (
    Backoff,
    Outbox,
    Reply,
    Verdict,
    classify,
    flush,
)

ATTEMPT = "0123456789abcdef0123456789abcdef"


def _request(attempt_id: str = ATTEMPT) -> dict[str, object]:
    return {"attempt_id": attempt_id, "failure": "crest: timeout"}


@pytest.mark.parametrize(
    ("status", "verdict"),
    [
        (200, Verdict.SENT),
        (None, Verdict.RETRY),  # no answer at all: network down
        (500, Verdict.RETRY),
        (502, Verdict.RETRY),
        (503, Verdict.RETRY),
        (408, Verdict.RETRY),
        (429, Verdict.RETRY),
        (401, Verdict.STOP),  # token revoked: nothing will ever be accepted
        (403, Verdict.DEAD),
        (404, Verdict.DEAD),
        (409, Verdict.DEAD),
        (422, Verdict.DEAD),
    ],
)
def test_each_answer_has_one_verdict(status: int | None, verdict: Verdict) -> None:
    assert classify(status) is verdict


def test_a_queued_result_survives_a_new_process(tmp_path: Path) -> None:
    Outbox(tmp_path).put(ATTEMPT, "a", _request())

    (entry,) = Outbox(tmp_path).pending()

    assert entry.attempt_id == ATTEMPT and entry.mol_id == "a"
    assert entry.request == _request()


def test_queueing_again_replaces_and_leaves_no_temporary_file(tmp_path: Path) -> None:
    box = Outbox(tmp_path)
    box.put(ATTEMPT, "a", _request())
    box.put(ATTEMPT, "a", {"attempt_id": ATTEMPT, "failure": "other"})

    assert [e.request["failure"] for e in box.pending()] == ["other"]
    assert [p.name for p in (tmp_path / "outbox").iterdir()] == [f"{ATTEMPT}.json"]


def test_pending_comes_out_oldest_first(tmp_path: Path) -> None:
    box = Outbox(tmp_path)
    ids = [f"{n:032x}" for n in (3, 1, 2)]
    for attempt_id in ids:
        box.put(attempt_id, "m", _request(attempt_id))

    assert [e.attempt_id for e in box.pending()] == ids


@pytest.mark.parametrize("bad", ["../escape", "a/b", "", "x" * 65, "a.json"])
def test_an_attempt_id_that_is_not_a_plain_name_is_refused(
    tmp_path: Path, bad: str
) -> None:
    with pytest.raises(ValueError):
        Outbox(tmp_path).put(bad, "a", _request())


def test_flush_removes_what_the_central_settled(tmp_path: Path) -> None:
    box = Outbox(tmp_path)
    box.put(ATTEMPT, "a", _request())
    sent: list[dict[str, object]] = []

    def send(request: dict[str, object]) -> Reply:
        sent.append(request)
        return Reply(200, {"outcome": "duplicate"})

    report = flush(box, send)

    assert sent == [_request()] and report.sent == 1 and box.pending() == []
    assert not report.retry and not report.stop


def test_flush_keeps_everything_and_stops_at_the_first_network_error(
    tmp_path: Path,
) -> None:
    box = Outbox(tmp_path)
    for n in range(3):
        box.put(f"{n:032x}", "m", _request(f"{n:032x}"))
    calls = 0

    def send(_request: dict[str, object]) -> Reply:
        nonlocal calls
        calls += 1
        return Reply(None, {}, "ConnectError: refused")

    report = flush(box, send)

    assert calls == 1 and report.retry and len(box.pending()) == 3


def test_a_refused_result_moves_to_the_dead_letter_with_the_reason(
    tmp_path: Path,
) -> None:
    box = Outbox(tmp_path)
    box.put(ATTEMPT, "a", _request())
    box.put("f" * 32, "b", _request("f" * 32))

    def send(request: dict[str, object]) -> Reply:
        if request["attempt_id"] == ATTEMPT:
            return Reply(422, {"detail": "Assinatura diferente"})
        return Reply(200, {"outcome": "accepted"})

    report = flush(box, send)

    assert report.dead == 1 and report.sent == 1 and box.pending() == []
    (dead,) = box.dead_letters()
    assert dead.attempt_id == ATTEMPT and dead.status == 422
    assert "Assinatura diferente" in dead.detail and dead.request == _request()


def test_a_revoked_token_keeps_the_result_and_stops(tmp_path: Path) -> None:
    box = Outbox(tmp_path)
    box.put(ATTEMPT, "a", _request())

    report = flush(box, lambda _r: Reply(401, {"detail": "revogado"}))

    assert report.stop and len(box.pending()) == 1 and box.dead_letters() == []


def test_counts_are_read_from_disk(tmp_path: Path) -> None:
    box = Outbox(tmp_path)
    box.put(ATTEMPT, "a", _request())
    box.put("f" * 32, "b", _request("f" * 32))
    flush(box, lambda r: Reply(404 if r["attempt_id"] == ATTEMPT else 503, {}))

    assert Outbox(tmp_path).counts() == (1, 1)
    files = json.loads((tmp_path / "dead-letter" / f"{ATTEMPT}.json").read_text())
    assert files["status"] == 404


def test_backoff_grows_to_a_cap_with_jitter() -> None:
    backoff = Backoff(initial_s=5.0, cap_s=300.0, rng=random.Random(0))  # noqa: S311

    delays = [backoff.delay(n) for n in range(1, 12)]

    for failures, delay in enumerate(delays, start=1):
        ceiling = min(300.0, 5.0 * 2 ** (failures - 1))
        assert ceiling / 2 <= delay <= ceiling
    assert max(delays) <= 300.0 and delays[-1] >= 150.0
