"""Tests for WorkerRegistry."""

import time

from grimperium.server.worker_registry import WorkerRegistry, make_worker_registry


class TestWorkerRegistryRegister:
    def test_register_new_worker(self) -> None:
        reg = WorkerRegistry()
        entry = reg.register("IQ", "iq-host")
        assert entry.worker_id == "IQ"
        assert entry.hostname == "iq-host"
        assert entry.processed == 0
        assert entry.shutdown_requested is False

    def test_register_twice_updates_hostname(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "old-host")
        entry = reg.register("IQ", "new-host")
        assert entry.hostname == "new-host"
        assert len(reg) == 1

    def test_evict_removes_worker(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        reg.evict("IQ")
        assert reg.get_worker("IQ") is None
        assert len(reg) == 0

    def test_evict_unknown_is_noop(self) -> None:
        reg = WorkerRegistry()
        reg.evict("nonexistent")


class TestWorkerRegistryHeartbeat:
    def test_heartbeat_updates_last_seen(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        first = reg.get_worker("IQ")
        assert first is not None
        t0 = first.last_seen
        time.sleep(0.01)
        reg.heartbeat("IQ")
        assert first.last_seen > t0

    def test_heartbeat_unknown_returns_false(self) -> None:
        reg = WorkerRegistry()
        assert reg.heartbeat("ghost") is False

    def test_heartbeat_known_returns_true(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        assert reg.heartbeat("IQ") is True


class TestWorkerRegistryMolTracking:
    def test_set_and_clear_current_mol(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        reg.set_current_mol("IQ", "mol_001")
        entry = reg.get_worker("IQ")
        assert entry is not None
        assert entry.current_mol_id == "mol_001"
        reg.clear_current_mol("IQ")
        assert entry.current_mol_id is None

    def test_get_worker_for_mol(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        reg.set_current_mol("IQ", "mol_007")
        assert reg.get_worker_for_mol("mol_007") == "IQ"
        assert reg.get_worker_for_mol("mol_999") is None

    def test_set_current_mol_unknown_worker_noop(self) -> None:
        reg = WorkerRegistry()
        reg.set_current_mol("ghost", "mol_001")


class TestWorkerRegistryMetrics:
    def test_record_success_increments(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        reg.set_current_mol("IQ", "mol_001")
        reg.record_success("IQ")
        e = reg.get_worker("IQ")
        assert e is not None
        assert e.processed == 1
        assert e.successful == 1
        assert e.failed == 0
        assert e.current_mol_id is None

    def test_record_failure_increments(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        reg.record_failure("IQ")
        e = reg.get_worker("IQ")
        assert e is not None
        assert e.processed == 1
        assert e.failed == 1

    def test_record_skip_increments(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        reg.record_skip("IQ")
        e = reg.get_worker("IQ")
        assert e is not None
        assert e.processed == 1
        assert e.skipped == 1


class TestWorkerRegistryConfig:
    def test_set_and_get_config(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        cfg = {"batch_size": 5, "crest_timeout_minutes": 30}
        assert reg.set_config("IQ", cfg) is True
        assert reg.get_config("IQ") == cfg

    def test_set_config_unknown_returns_false(self) -> None:
        reg = WorkerRegistry()
        assert reg.set_config("ghost", {}) is False

    def test_get_config_unknown_returns_none(self) -> None:
        reg = WorkerRegistry()
        assert reg.get_config("ghost") is None

    def test_get_config_before_set_returns_none(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        assert reg.get_config("IQ") is None


class TestWorkerRegistryShutdown:
    def test_request_shutdown_known(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        assert reg.request_shutdown("IQ") is True
        assert reg.is_shutdown_requested("IQ") is True

    def test_request_shutdown_unknown(self) -> None:
        reg = WorkerRegistry()
        assert reg.request_shutdown("ghost") is False

    def test_is_shutdown_not_set_by_default(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        assert reg.is_shutdown_requested("IQ") is False

    def test_is_shutdown_unknown_returns_false(self) -> None:
        reg = WorkerRegistry()
        assert reg.is_shutdown_requested("ghost") is False

    def test_request_shutdown_all(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h1")
        reg.register("remote1", "h2")
        evicted = reg.request_shutdown_all()
        assert set(evicted) == {"IQ", "remote1"}
        assert reg.is_shutdown_requested("IQ") is True
        assert reg.is_shutdown_requested("remote1") is True


class TestWorkerRegistryQueries:
    def test_all_workers_empty(self) -> None:
        reg = WorkerRegistry()
        assert reg.all_workers() == []

    def test_all_workers_returns_snapshot(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h1")
        reg.register("remote1", "h2")
        workers = reg.all_workers()
        assert len(workers) == 2
        ids = {w.worker_id for w in workers}
        assert ids == {"IQ", "remote1"}

    def test_len(self) -> None:
        reg = WorkerRegistry()
        assert len(reg) == 0
        reg.register("IQ", "h")
        assert len(reg) == 1


class TestWorkerRegistryCompatViews:
    def test_heartbeat_registry_view(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "iq-host")
        hr = reg.heartbeat_registry
        assert "IQ" in hr
        hostname, last_seen = hr["IQ"]
        assert hostname == "iq-host"

    def test_running_molecules_view(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        reg.set_current_mol("IQ", "mol_001")
        rm = reg.running_molecules
        assert rm == {"mol_001": "IQ"}

    def test_running_molecules_empty_when_none_active(self) -> None:
        reg = WorkerRegistry()
        reg.register("IQ", "h")
        assert reg.running_molecules == {}


def test_make_worker_registry() -> None:
    reg = make_worker_registry()
    assert isinstance(reg, WorkerRegistry)
    assert len(reg) == 0
