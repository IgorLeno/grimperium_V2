"""Tests for worker offline result queue with stable result_id."""

from __future__ import annotations

from pathlib import Path

from grimperium.worker.local_store import LocalStore
from grimperium.worker.offline_queue import OfflineResultQueue


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
