"""Tests for worker/local_store.py — LocalStore in-memory tracking."""

import pytest

from grimperium.worker.local_store import LocalRecord, LocalStore


class TestLocalRecord:
    def test_defaults(self) -> None:
        r = LocalRecord(mol_id="m1", smiles="CCO")
        assert r.mol_id == "m1"
        assert r.smiles == "CCO"
        assert r.completed is False
        assert r.success is None
        assert r.error is None
        assert r.result_update == {}
        assert r.completed_at is None
        assert r.claimed_at is not None


class TestLocalStore:
    def test_add_and_get(self) -> None:
        store = LocalStore()
        record = store.add("m1", "CCO")
        assert record.mol_id == "m1"
        assert store.get("m1") is record

    def test_get_missing_returns_none(self) -> None:
        store = LocalStore()
        assert store.get("ghost") is None

    def test_add_duplicate_raises(self) -> None:
        store = LocalStore()
        store.add("m1", "CCO")
        with pytest.raises(ValueError, match="m1"):
            store.add("m1", "CCO")

    def test_len(self) -> None:
        store = LocalStore()
        assert len(store) == 0
        store.add("m1", "CCO")
        store.add("m2", "CCC")
        assert len(store) == 2

    def test_pending_excludes_completed(self) -> None:
        store = LocalStore()
        store.add("m1", "CCO")
        store.add("m2", "CCC")
        store.mark_success("m1", {"H298_pm7": -1.0})
        pending = store.pending()
        assert len(pending) == 1
        assert pending[0].mol_id == "m2"

    def test_pending_empty_when_all_done(self) -> None:
        store = LocalStore()
        store.add("m1", "CCO")
        store.mark_failure("m1", "timeout")
        assert store.pending() == []

    def test_mark_success(self) -> None:
        store = LocalStore()
        store.add("m1", "CCO")
        store.mark_success("m1", {"H298_pm7": -42.0})
        r = store.get("m1")
        assert r is not None
        assert r.completed is True
        assert r.success is True
        assert r.result_update == {"H298_pm7": -42.0}
        assert r.error is None
        assert r.completed_at is not None

    def test_mark_failure(self) -> None:
        store = LocalStore()
        store.add("m1", "CCO")
        store.mark_failure("m1", "CREST timeout")
        r = store.get("m1")
        assert r is not None
        assert r.completed is True
        assert r.success is False
        assert r.error == "CREST timeout"
        assert r.completed_at is not None

    def test_mark_success_missing_key_raises(self) -> None:
        store = LocalStore()
        with pytest.raises(KeyError):
            store.mark_success("ghost", {})

    def test_mark_failure_missing_key_raises(self) -> None:
        store = LocalStore()
        with pytest.raises(KeyError):
            store.mark_failure("ghost", "error")

    def test_completed_returns_only_done(self) -> None:
        store = LocalStore()
        store.add("m1", "CCO")
        store.add("m2", "CCC")
        store.mark_success("m2", {})
        done = store.completed()
        assert len(done) == 1
        assert done[0].mol_id == "m2"

    def test_clear_completed_removes_and_returns_count(self) -> None:
        store = LocalStore()
        store.add("m1", "CCO")
        store.add("m2", "CCC")
        store.mark_success("m1", {})
        n = store.clear_completed()
        assert n == 1
        assert len(store) == 1
        assert store.get("m1") is None
        assert store.get("m2") is not None

    def test_clear_completed_when_none_done(self) -> None:
        store = LocalStore()
        store.add("m1", "CCO")
        assert store.clear_completed() == 0
        assert len(store) == 1

    def test_clear_all_completed(self) -> None:
        store = LocalStore()
        store.add("m1", "CCO")
        store.add("m2", "CCC")
        store.mark_success("m1", {})
        store.mark_failure("m2", "err")
        n = store.clear_completed()
        assert n == 2
        assert len(store) == 0
