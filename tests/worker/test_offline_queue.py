"""Tests for worker offline result queue with stable result_id."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import pytest

from grimperium.worker.local_store import LocalStore
from grimperium.worker.offline_queue import (
    DeadLetterQueue,
    DeadLetterRecord,
    OfflineResultQueue,
    PendingAbortQueue,
    compute_dead_letter_id,
    dead_letter_path_for,
    pending_abort_path_for,
)


def test_local_store_assigns_stable_result_id() -> None:
    store = LocalStore()
    store.add("m1", "CCO")
    store.mark_success("m1", {"H298_pm7": -1.0})
    record = store.get("m1")
    assert record is not None
    assert record.result_id is not None
    first_id = record.result_id
    # Re-marking must not rotate the ID.
    store.mark_success("m1", {"H298_pm7": -1.0})
    assert store.get("m1") is not None
    assert store.get("m1").result_id == first_id  # type: ignore[union-attr]


def test_offline_queue_persists_and_reloads(tmp_path: Path) -> None:
    path = tmp_path / "offline.jsonl"
    queue = OfflineResultQueue(path)
    entry = queue.enqueue(
        mol_id="m1",
        success=False,
        error="timeout",
        result_id="fixed-id",
        completed_at="2026-04-21T10:00:00Z",
    )
    assert entry.result_id == "fixed-id"
    reloaded = OfflineResultQueue(path)
    pending = reloaded.pending()
    assert len(pending) == 1
    assert pending[0].result_id == "fixed-id"
    assert pending[0].to_sync_dict()["result_id"] == "fixed-id"
    reloaded.confirm("fixed-id")
    assert OfflineResultQueue(path).pending() == []


def _sample_dead_letter_record(result_id: str = "rid-1") -> DeadLetterRecord:
    payload = {
        "result_id": result_id,
        "mol_id": "m1",
        "success": False,
        "error": "x",
        "completed_at": "2026-07-13T00:00:00Z",
    }
    return DeadLetterRecord(
        result_id=result_id,
        mol_id="m1",
        attempt_id="att-1",
        original_payload=payload,
        returned_status="conflict",
        detail="x",
        worker_id="w1",
        rejected_at="2026-07-13T00:00:00Z",
        rejection_origin="sync_results",
    )


def test_dead_letter_append_persist_failure_leaves_memory_and_disk_unchanged(
    tmp_path: Path,
) -> None:
    path = tmp_path / "dead_letter.jsonl"
    queue = DeadLetterQueue(path)
    record = _sample_dead_letter_record()
    with (
        patch(
            "grimperium.worker.offline_queue._atomic_write_text",
            side_effect=OSError("disk full"),
        ),
        pytest.raises(OSError),
    ):
        queue.append(record)
    assert queue.entries() == []
    assert not path.exists() or path.read_text(encoding="utf-8") == ""


def test_dead_letter_append_retry_after_persist_failure(tmp_path: Path) -> None:
    path = tmp_path / "dead_letter.jsonl"
    queue = DeadLetterQueue(path)
    record = _sample_dead_letter_record()
    attempts = {"count": 0}
    real_write = __import__(
        "grimperium.worker.offline_queue", fromlist=["_atomic_write_text"]
    )._atomic_write_text

    def flaky_write(target_path: Path, content: str) -> None:
        attempts["count"] += 1
        if attempts["count"] == 1:
            raise OSError("disk full")
        real_write(target_path, content)

    with patch(
        "grimperium.worker.offline_queue._atomic_write_text",
        side_effect=flaky_write,
    ):
        with pytest.raises(OSError):
            queue.append(record)
        assert queue.entries() == []
        queue.append(record)
    assert len(queue.entries()) == 1
    reloaded = DeadLetterQueue(path)
    assert len(reloaded.entries()) == 1


def test_dead_letter_append_preserves_existing_on_failure(tmp_path: Path) -> None:
    path = tmp_path / "dead_letter.jsonl"
    queue = DeadLetterQueue(path)
    queue.append(_sample_dead_letter_record("existing"))
    with (
        patch(
            "grimperium.worker.offline_queue._atomic_write_text",
            side_effect=OSError("disk full"),
        ),
        pytest.raises(OSError),
    ):
        queue.append(_sample_dead_letter_record("new"))
    assert len(queue.entries()) == 1
    assert queue.entries()[0].result_id == "existing"
    reloaded = DeadLetterQueue(path)
    assert len(reloaded.entries()) == 1


def test_dead_letter_append_idempotent_after_success(tmp_path: Path) -> None:
    path = tmp_path / "dead_letter.jsonl"
    queue = DeadLetterQueue(path)
    record = _sample_dead_letter_record()
    queue.append(record)
    queue.append(record)
    assert len(queue.entries()) == 1


def test_pending_abort_queue_persists_and_reloads(tmp_path: Path) -> None:
    queue_path = tmp_path / "offline.jsonl"
    path = pending_abort_path_for(queue_path)
    queue = PendingAbortQueue(path)
    record = replace(
        _sample_dead_letter_record("pending-1"), rejection_origin="lease_lost"
    )
    queue.append(record)
    reloaded = PendingAbortQueue(path)
    assert len(reloaded.entries()) == 1
    dead_letter_id = compute_dead_letter_id(
        result_id=record.result_id,
        returned_status=record.returned_status,
        original_payload=record.original_payload,
    )
    reloaded.remove(dead_letter_id)
    assert PendingAbortQueue(path).entries() == []


def test_pending_abort_path_for_sibling(tmp_path: Path) -> None:
    queue_path = tmp_path / "w1_offline_results.jsonl"
    assert (
        pending_abort_path_for(queue_path).name
        == "w1_offline_results_pending_aborts.jsonl"
    )
    assert (
        dead_letter_path_for(queue_path).name == "w1_offline_results_dead_letter.jsonl"
    )
