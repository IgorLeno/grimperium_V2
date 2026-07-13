"""Tests for worker/runner.py — WorkerRunner main processing loop."""

import json
import tempfile
import time
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pandas as pd

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.worker.client import WorkerClient
from grimperium.worker.offline_queue import dead_letter_path_for
from grimperium.worker.runner import WorkerConfig, WorkerRunner

# ── Helpers ───────────────────────────────────────────────────────────────────


def _make_config(**kwargs: Any) -> WorkerConfig:
    defaults: dict[str, Any] = {
        "server_url": "http://test-server",
        "worker_id": "w1",
        "poll_interval_s": 0.0,
        "heartbeat_interval_s": 999.0,  # prevent actual heartbeat in tests
        "max_idle_polls": 3,
        "crest_timeout_minutes": 30,
        "mopac_timeout_minutes": 10,
        "batch_id": "test-batch",
        # Isolar fila offline por teste para evitar contaminação entre casos.
        "offline_queue_path": str(
            Path(tempfile.mkdtemp(prefix="grimperium-worker-test-"))
            / "offline_results.jsonl"
        ),
    }
    defaults.update(kwargs)
    return WorkerConfig(**defaults)


def _echo_sync_results(results: list[dict[str, Any]]) -> dict[str, Any]:
    """Mock de sync que confirma exatamente os result_ids enviados."""
    return {
        "accepted": len(results),
        "rejected": 0,
        "duplicate": False,
        "items": [
            {
                "result_id": str(item.get("result_id", "")),
                "mol_id": str(item.get("mol_id", "")),
                "status": "applied",
            }
            for item in results
        ],
    }


def _mock_client(
    claim_returns: tuple[str, str, str | None] | None = None,
) -> MagicMock:
    client = MagicMock(spec=WorkerClient)
    client.claim.return_value = claim_returns
    client.sync_results.side_effect = _echo_sync_results
    return client


def _mock_pipeline(success: bool = True, error_msg: str | None = None) -> MagicMock:
    from grimperium.crest_pm7 import CRESTPM7Pipeline, PM7Result

    result = MagicMock(spec=PM7Result)
    result.mol_id = "m1"
    result.success = success
    result.error_message = error_msg
    result.most_stable_hof = -42.0 if success else None

    pipeline = MagicMock(spec=CRESTPM7Pipeline)
    pipeline.process_molecule.return_value = result
    return pipeline


def _write_minimal_batch_csv(csv_path: Path, mol_id: str = "m1") -> None:
    schema = BatchCSVManager(csv_path=None).get_schema()
    row: dict[str, Any] = dict.fromkeys(schema)
    row.update(
        {
            "mol_id": mol_id,
            "status": "Pending",
            "smiles": "CCO",
            "nheavy": 3,
            "reruns": 0,
        }
    )
    pd.DataFrame([row], columns=schema).to_csv(csv_path, index=False)


# ── run_one ───────────────────────────────────────────────────────────────────


class TestRunOne:
    def test_returns_false_when_queue_empty(self) -> None:
        client = _mock_client(claim_returns=None)
        pipeline = _mock_pipeline()
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        assert runner.run_one() is False
        pipeline.process_molecule.assert_not_called()

    @patch(
        "grimperium.worker.runner._pm7result_to_update",
        return_value={"H298_pm7": -42.0},
    )
    def test_returns_true_on_success(self, _mock_update: MagicMock) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = _mock_pipeline(success=True)
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        assert runner.run_one() is True

    @patch(
        "grimperium.worker.runner._pm7result_to_update",
        return_value={"H298_pm7": -42.0},
    )
    def test_calls_process_molecule_with_claimed_args(
        self, _mock_update: MagicMock
    ) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = _mock_pipeline(success=True)
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        runner.run_one()
        pipeline.process_molecule.assert_called_once_with("m1", "CCO")

    @patch(
        "grimperium.worker.runner._pm7result_to_update",
        return_value={"H298_pm7": -42.0},
    )
    def test_syncs_success_to_server(self, _mock_update: MagicMock) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = _mock_pipeline(success=True)
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        runner.run_one()
        client.sync_results.assert_called()
        client.report_success.assert_not_called()
        client.report_failure.assert_not_called()

    def test_syncs_failure_when_pipeline_result_not_success(self) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = _mock_pipeline(success=False, error_msg="CREST timeout")
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        runner.run_one()
        client.sync_results.assert_called()
        payload = client.sync_results.call_args[0][0][0]
        assert payload["success"] is False
        assert payload["error"] == "CREST timeout"
        client.report_failure.assert_not_called()
        client.report_success.assert_not_called()

    def test_syncs_failure_with_fallback_message_when_no_error(self) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = _mock_pipeline(success=False, error_msg=None)
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        runner.run_one()
        payload = client.sync_results.call_args[0][0][0]
        assert payload["mol_id"] == "m1"
        assert isinstance(payload["error"], str) and len(payload["error"]) > 0

    def test_syncs_failure_when_pipeline_raises(self) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = MagicMock()
        pipeline.process_molecule.side_effect = RuntimeError("MOPAC crashed")
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        runner.run_one()
        client.sync_results.assert_called()
        payload = client.sync_results.call_args[0][0][0]
        assert payload["mol_id"] == "m1"
        assert "MOPAC crashed" in str(payload["error"])

    @patch("grimperium.worker.runner._pm7result_to_update", return_value={})
    def test_store_cleared_after_success(self, _mock_update: MagicMock) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = _mock_pipeline(success=True)
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        runner.run_one()
        assert len(runner._store) == 0

    def test_store_cleared_after_failure(self) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = _mock_pipeline(success=False, error_msg="err")
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        runner.run_one()
        assert len(runner._store) == 0

    def test_update_csv_no_op_when_no_csv_path(self) -> None:
        runner = WorkerRunner(
            _make_config(csv_path=None),
            pipeline=_mock_pipeline(),
            client=_mock_client(),
        )
        runner._update_csv("m1", {"status": "Running"})

    def test_update_csv_writes_to_file(self, tmp_path: Path) -> None:
        csv_path = tmp_path / "batch.csv"
        _write_minimal_batch_csv(csv_path)
        runner = WorkerRunner(
            _make_config(csv_path=str(csv_path)),
            pipeline=_mock_pipeline(),
            client=_mock_client(),
        )

        runner._update_csv("m1", {"status": "Running"})

        df = pd.read_csv(csv_path)
        status = df.loc[df["mol_id"] == "m1", "status"].iloc[0]
        assert status == "Running"

    @patch(
        "grimperium.worker.runner._pm7result_to_update",
        return_value={"H298_pm7": -42.0},
    )
    def test_run_one_updates_csv_running_before_pipeline(
        self, _mock_update: MagicMock
    ) -> None:
        events: list[str] = []
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = _mock_pipeline(success=True)

        def process_molecule(_mol_id: str, _smiles: str) -> MagicMock:
            events.append("pipeline")
            return pipeline.process_molecule.return_value

        pipeline.process_molecule.side_effect = process_molecule
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)

        def update_csv(_mol_id: str, updates: dict[str, Any]) -> None:
            if updates.get("status") == "Running":
                events.append("csv_running")

        runner._update_csv = MagicMock(side_effect=update_csv)

        runner.run_one()

        assert events[:2] == ["csv_running", "pipeline"]


# ── run ───────────────────────────────────────────────────────────────────────


class TestRun:
    @patch("grimperium.worker.runner._pm7result_to_update", return_value={})
    def test_run_returns_count_of_processed(self, _mock_update: MagicMock) -> None:
        client = MagicMock(spec=WorkerClient)
        client.sync_results.side_effect = _echo_sync_results
        client.claim.side_effect = [
            ("m1", "CCO", "att-1"),
            ("m2", "CCC", "att-2"),
            None,
            None,
            None,  # 3 empty to trigger idle stop
        ]
        pipeline = _mock_pipeline(success=True)
        runner = WorkerRunner(
            _make_config(max_idle_polls=3), pipeline=pipeline, client=client
        )
        count = runner.run()
        assert count == 2

    @patch("grimperium.worker.runner._pm7result_to_update", return_value={})
    def test_run_respects_max_molecules(self, _mock_update: MagicMock) -> None:
        client = MagicMock(spec=WorkerClient)
        client.sync_results.side_effect = _echo_sync_results
        client.claim.return_value = ("m1", "CCO", "att-1")
        pipeline = _mock_pipeline(success=True)
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)
        count = runner.run(max_molecules=2)
        assert count == 2
        assert client.claim.call_count == 2

    def test_run_stops_after_max_idle_polls(self) -> None:
        client = MagicMock(spec=WorkerClient)
        client.sync_results.side_effect = _echo_sync_results
        client.claim.return_value = None
        pipeline = _mock_pipeline()
        runner = WorkerRunner(
            _make_config(max_idle_polls=3), pipeline=pipeline, client=client
        )
        count = runner.run()
        assert count == 0
        assert client.claim.call_count == 3

    def test_run_calls_register_before_first_claim(self) -> None:
        client = _mock_client(claim_returns=None)
        pipeline = _mock_pipeline()
        runner = WorkerRunner(
            _make_config(max_idle_polls=1), pipeline=pipeline, client=client
        )
        runner.run()
        client.register.assert_called_once()
        methods = [call[0] for call in client.method_calls]
        assert methods.index("register") < methods.index("claim")

    @patch("grimperium.worker.runner._pm7result_to_update", return_value={})
    def test_run_stops_via_stop_method(self, _mock_update: MagicMock) -> None:
        """stop() between run_one calls should terminate the loop."""
        config = _make_config(max_idle_polls=100)
        client = MagicMock(spec=WorkerClient)
        client.sync_results.side_effect = _echo_sync_results
        pipeline = _mock_pipeline(success=True)
        runner = WorkerRunner(config, pipeline=pipeline, client=client)

        call_count = 0

        def side_effect() -> tuple[str, str, str | None] | None:
            nonlocal call_count
            call_count += 1
            if call_count == 2:
                runner.stop()
            return ("m1", "CCO", "att-1")

        client.claim.side_effect = side_effect
        count = runner.run()
        assert count >= 1


# ── consecutive failure stop ────────────────────────────────────────────────


class TestConsecutiveFailureStop:
    @patch("grimperium.worker.runner._pm7result_to_update", return_value={})
    def test_counter_resets_after_success(self, _mock_update: MagicMock) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        fail_pipeline = _mock_pipeline(success=False, error_msg="err")
        ok_pipeline = _mock_pipeline(success=True)
        runner = WorkerRunner(_make_config(), pipeline=fail_pipeline, client=client)

        runner.run_one()
        runner.run_one()
        runner._pipeline = ok_pipeline
        runner.run_one()

        assert runner._consecutive_failures == 0

    def test_counter_increments_after_pipeline_failure(self) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = _mock_pipeline(success=False, error_msg="err")
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)

        runner.run_one()
        runner.run_one()
        runner.run_one()

        assert runner._consecutive_failures == 3

    def test_counter_increments_after_exception(self) -> None:
        client = _mock_client(claim_returns=("m1", "CCO", "att-1"))
        pipeline = MagicMock()
        pipeline.process_molecule.side_effect = RuntimeError("MOPAC crashed")
        runner = WorkerRunner(_make_config(), pipeline=pipeline, client=client)

        runner.run_one()
        runner.run_one()

        assert runner._consecutive_failures == 2

    def test_empty_queue_does_not_increment_counter(self) -> None:
        client = MagicMock(spec=WorkerClient)
        client.sync_results.side_effect = _echo_sync_results
        client.claim.return_value = None
        pipeline = _mock_pipeline()
        runner = WorkerRunner(
            _make_config(max_idle_polls=5), pipeline=pipeline, client=client
        )

        runner.run()

        assert runner._consecutive_failures == 0

    def test_run_stops_after_max_consecutive_failures(self) -> None:
        client = MagicMock(spec=WorkerClient)
        client.sync_results.side_effect = _echo_sync_results
        client.claim.return_value = ("m1", "CCO", "att-1")
        pipeline = _mock_pipeline(success=False, error_msg="boom")
        cfg = _make_config(
            consecutive_failure_stop=True,
            max_consecutive_failures=3,
        )
        runner = WorkerRunner(cfg, pipeline=pipeline, client=client)

        count = runner.run()

        assert client.claim.call_count == 3
        assert count == 0

    def test_consecutive_failure_stop_false_disables_stop(self) -> None:
        client = MagicMock(spec=WorkerClient)
        client.sync_results.side_effect = _echo_sync_results
        client.claim.return_value = ("m1", "CCO", "att-1")
        pipeline = _mock_pipeline(success=False, error_msg="boom")
        cfg = _make_config(
            consecutive_failure_stop=False,
            max_consecutive_failures=3,
        )
        runner = WorkerRunner(cfg, pipeline=pipeline, client=client)

        count = runner.run(max_molecules=20)

        assert client.claim.call_count == 20
        assert count == 0

    @patch("grimperium.worker.runner._pm7result_to_update", return_value={})
    def test_success_between_failures_resets_stop_threshold(
        self, _mock_update: MagicMock
    ) -> None:
        client = MagicMock(spec=WorkerClient)
        client.sync_results.side_effect = _echo_sync_results
        client.claim.return_value = ("m1", "CCO", "att-1")
        cfg = _make_config(max_consecutive_failures=3)
        runner = WorkerRunner(cfg, pipeline=_mock_pipeline(), client=client)
        outcomes = [False, True, False, False, False]

        def process_molecule(_mol_id: str, _smiles: str) -> MagicMock:
            success = outcomes.pop(0)
            result = MagicMock()
            result.mol_id = "m1"
            result.success = success
            result.error_message = None if success else "boom"
            result.most_stable_hof = -42.0 if success else None
            return result

        runner._pipeline.process_molecule.side_effect = process_molecule

        count = runner.run()

        assert client.claim.call_count == 5
        assert count == 1
        assert runner._consecutive_failures == 3
        assert runner._last_run_succeeded is False

    @patch("grimperium.worker.runner._pm7result_to_update", return_value={})
    def test_processed_counts_all_successes_after_failure(
        self, _mock_update: MagicMock
    ) -> None:
        """F, S, S, S -> processed deve ser 3, nao 1."""
        client = MagicMock(spec=WorkerClient)
        client.sync_results.side_effect = _echo_sync_results
        client.claim.return_value = ("m1", "CCO", "att-1")
        cfg = _make_config(max_consecutive_failures=10)
        runner = WorkerRunner(cfg, pipeline=_mock_pipeline(), client=client)
        outcomes = [False, True, True, True]

        def process_molecule(_mol_id: str, _smiles: str) -> MagicMock:
            success = outcomes.pop(0)
            result = MagicMock()
            result.mol_id = "m1"
            result.success = success
            result.error_message = None if success else "boom"
            result.most_stable_hof = -42.0 if success else None
            return result

        runner._pipeline.process_molecule.side_effect = process_molecule

        count = runner.run(max_molecules=4)

        assert count == 3
        assert runner._consecutive_failures == 0
        assert runner._last_run_succeeded is True


# ── WorkerConfig ─────────────────────────────────────────────────────────────


class TestWorkerConfig:
    def test_heartbeat_interval_default_is_30s(self) -> None:
        cfg = WorkerConfig(server_url="http://x", worker_id="w1")
        assert cfg.heartbeat_interval_s == 30.0

    def test_batch_size_default(self) -> None:
        cfg = WorkerConfig(server_url="http://x", worker_id="w1")
        assert cfg.batch_size == 10

    def test_batch_size_custom(self) -> None:
        cfg = WorkerConfig(server_url="http://x", worker_id="w1", batch_size=5)
        assert cfg.batch_size == 5

    def test_csv_path_default(self) -> None:
        cfg = WorkerConfig(server_url="http://x", worker_id="w1")
        assert cfg.csv_path is None


class TestWorkerRunnerReconfigure:
    def test_reconfigure_empty_dict_is_noop(self) -> None:
        cfg = _make_config(crest_timeout_minutes=30, mopac_timeout_minutes=10)
        runner = WorkerRunner(cfg, pipeline=_mock_pipeline(), client=_mock_client())
        runner.reconfigure({})
        assert runner._config.crest_timeout_minutes == 30

    def test_reconfigure_updates_batch_size(self) -> None:
        cfg = _make_config()
        runner = WorkerRunner(cfg, pipeline=_mock_pipeline(), client=_mock_client())
        runner.reconfigure({"batch_size": 7})
        assert runner._config.batch_size == 7

    def test_reconfigure_updates_crest_timeout(self) -> None:
        cfg = _make_config(crest_timeout_minutes=30)
        runner = WorkerRunner(cfg, pipeline=_mock_pipeline(), client=_mock_client())
        runner.reconfigure({"crest_timeout_minutes": 90})
        assert runner._config.crest_timeout_minutes == 90

    def test_reconfigure_updates_mopac_timeout(self) -> None:
        cfg = _make_config(mopac_timeout_minutes=10)
        runner = WorkerRunner(cfg, pipeline=_mock_pipeline(), client=_mock_client())
        runner.reconfigure({"mopac_timeout_minutes": 45})
        assert runner._config.mopac_timeout_minutes == 45

    def test_reconfigure_rebuilds_pipeline_when_timeout_changes(self) -> None:
        cfg = _make_config(crest_timeout_minutes=30)
        old_pipeline = _mock_pipeline()
        runner = WorkerRunner(cfg, pipeline=old_pipeline, client=_mock_client())
        runner.reconfigure({"crest_timeout_minutes": 90})
        # Pipeline was replaced (not the same mock object)
        assert runner._pipeline is not old_pipeline

    def test_reconfigure_no_rebuild_when_same_timeout(self) -> None:
        cfg = _make_config(crest_timeout_minutes=30)
        old_pipeline = _mock_pipeline()
        runner = WorkerRunner(cfg, pipeline=old_pipeline, client=_mock_client())
        runner.reconfigure({"crest_timeout_minutes": 30})
        # Same value → no rebuild
        assert runner._pipeline is old_pipeline

    def test_reconfigure_no_rebuild_for_batch_size_only(self) -> None:
        cfg = _make_config()
        old_pipeline = _mock_pipeline()
        runner = WorkerRunner(cfg, pipeline=old_pipeline, client=_mock_client())
        runner.reconfigure({"batch_size": 5})
        assert runner._pipeline is old_pipeline


# ── _pm7result_to_update ──────────────────────────────────────────────────────


class TestPm7ResultToUpdate:
    def test_pm7result_to_update_returns_dict(self) -> None:
        from grimperium.crest_pm7 import PM7Result
        from grimperium.crest_pm7.config import CRESTStatus
        from grimperium.worker.runner import _pm7result_to_update

        result = PM7Result(mol_id="m1", smiles="CCO")
        result.crest_status = CRESTStatus.NOT_ATTEMPTED

        update = _pm7result_to_update(result, "w1", 30.0, 10.0)
        assert isinstance(update, dict)

    def test_pm7result_to_update_includes_timeouts(self) -> None:
        from grimperium.crest_pm7 import PM7Result
        from grimperium.worker.runner import _pm7result_to_update

        result = PM7Result(mol_id="m1", smiles="CCO")

        update = _pm7result_to_update(result, "w1", 45.0, 15.0)
        assert update.get("assigned_crest_timeout") == 45.0
        assert update.get("assigned_mopac_timeout") == 15.0


# ── Dead-letter (conflict / stale_attempt) ───────────────────────────────────


def _read_dead_letter(queue_path: Path) -> list[dict[str, Any]]:
    path = dead_letter_path_for(queue_path)
    if not path.exists():
        return []
    records: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line:
            records.append(json.loads(line))
    return records


class TestDeadLetterFlush:
    def test_conflict_writes_dead_letter_with_full_payload(
        self, tmp_path: Path
    ) -> None:
        queue_path = tmp_path / "w1_offline_results.jsonl"
        client = MagicMock(spec=WorkerClient)
        runner = WorkerRunner(
            _make_config(offline_queue_path=str(queue_path), worker_id="w1"),
            pipeline=_mock_pipeline(),
            client=client,
        )
        entry = runner._offline_queue.enqueue(
            mol_id="m1",
            success=True,
            result_update={"H298_pm7": -1.0},
            result_id="rid-conflict",
            attempt_id="att-1",
            completed_at="2026-07-12T10:00:00Z",
        )
        original = entry.to_sync_dict()

        def sync_side_effect(results: list[dict[str, Any]]) -> dict[str, Any]:
            return {
                "accepted": 0,
                "rejected": 1,
                "items": [
                    {
                        "result_id": "rid-conflict",
                        "mol_id": "m1",
                        "status": "conflict",
                        "detail": "fingerprint mismatch",
                    }
                ],
            }

        client.sync_results.side_effect = sync_side_effect
        runner.flush_offline_queue()

        assert runner._offline_queue.pending() == []
        dl = _read_dead_letter(queue_path)
        assert len(dl) == 1
        rec = dl[0]
        assert rec["result_id"] == "rid-conflict"
        assert rec["mol_id"] == "m1"
        assert rec["attempt_id"] == "att-1"
        assert rec["original_payload"] == original
        assert rec["returned_status"] == "conflict"
        assert rec["detail"] == "fingerprint mismatch"
        assert rec["worker_id"] == "w1"
        assert rec["rejected_at"]
        assert rec["rejection_origin"]

    def test_stale_attempt_writes_dead_letter(self, tmp_path: Path) -> None:
        queue_path = tmp_path / "w1_offline_results.jsonl"
        client = MagicMock(spec=WorkerClient)
        runner = WorkerRunner(
            _make_config(offline_queue_path=str(queue_path), worker_id="w1"),
            pipeline=_mock_pipeline(),
            client=client,
        )
        runner._offline_queue.enqueue(
            mol_id="m1",
            success=False,
            error="late",
            result_id="rid-stale",
            attempt_id="att-old",
            completed_at="2026-07-12T10:00:00Z",
        )
        client.sync_results.side_effect = lambda results: {
            "accepted": 0,
            "rejected": 1,
            "items": [
                {
                    "result_id": "rid-stale",
                    "mol_id": "m1",
                    "status": "stale_attempt",
                    "detail": "lease reclaimed",
                }
            ],
        }
        runner.flush_offline_queue()
        assert runner._offline_queue.pending() == []
        dl = _read_dead_letter(queue_path)
        assert len(dl) == 1
        assert dl[0]["returned_status"] == "stale_attempt"
        assert dl[0]["result_id"] == "rid-stale"

    def test_restart_keeps_dead_letter(self, tmp_path: Path) -> None:
        queue_path = tmp_path / "w1_offline_results.jsonl"
        client = MagicMock(spec=WorkerClient)
        runner = WorkerRunner(
            _make_config(offline_queue_path=str(queue_path), worker_id="w1"),
            pipeline=_mock_pipeline(),
            client=client,
        )
        runner._offline_queue.enqueue(
            mol_id="m1",
            success=True,
            result_update={},
            result_id="rid-1",
            attempt_id="att-1",
            completed_at="2026-07-12T10:00:00Z",
        )
        client.sync_results.side_effect = lambda results: {
            "accepted": 0,
            "rejected": 1,
            "items": [
                {
                    "result_id": "rid-1",
                    "mol_id": "m1",
                    "status": "conflict",
                    "detail": "x",
                }
            ],
        }
        runner.flush_offline_queue()
        assert len(_read_dead_letter(queue_path)) == 1

        restarted = WorkerRunner(
            _make_config(offline_queue_path=str(queue_path), worker_id="w1"),
            pipeline=_mock_pipeline(),
            client=client,
        )
        assert restarted._offline_queue.pending() == []
        assert len(_read_dead_letter(queue_path)) == 1

    def test_dead_letter_write_failure_keeps_item_in_queue(
        self, tmp_path: Path
    ) -> None:
        queue_path = tmp_path / "w1_offline_results.jsonl"
        client = MagicMock(spec=WorkerClient)
        runner = WorkerRunner(
            _make_config(offline_queue_path=str(queue_path), worker_id="w1"),
            pipeline=_mock_pipeline(),
            client=client,
        )
        runner._offline_queue.enqueue(
            mol_id="m1",
            success=True,
            result_update={},
            result_id="rid-fail-dl",
            attempt_id="att-1",
            completed_at="2026-07-12T10:00:00Z",
        )
        client.sync_results.side_effect = lambda results: {
            "accepted": 0,
            "rejected": 1,
            "items": [
                {
                    "result_id": "rid-fail-dl",
                    "mol_id": "m1",
                    "status": "conflict",
                    "detail": "x",
                }
            ],
        }
        with patch.object(
            runner._dead_letter,
            "append",
            side_effect=OSError("disk full"),
        ):
            runner.flush_offline_queue()
        pending = runner._offline_queue.pending()
        assert len(pending) == 1
        assert pending[0].result_id == "rid-fail-dl"
        assert _read_dead_letter(queue_path) == []

    def test_applied_and_duplicate_not_in_dead_letter(self, tmp_path: Path) -> None:
        queue_path = tmp_path / "w1_offline_results.jsonl"
        client = MagicMock(spec=WorkerClient)
        runner = WorkerRunner(
            _make_config(offline_queue_path=str(queue_path), worker_id="w1"),
            pipeline=_mock_pipeline(),
            client=client,
        )
        runner._offline_queue.enqueue(
            mol_id="m1",
            success=True,
            result_update={},
            result_id="rid-applied",
            completed_at="2026-07-12T10:00:00Z",
        )
        runner._offline_queue.enqueue(
            mol_id="m2",
            success=True,
            result_update={},
            result_id="rid-dup",
            completed_at="2026-07-12T10:00:00Z",
        )
        client.sync_results.side_effect = lambda results: {
            "accepted": 2,
            "rejected": 0,
            "items": [
                {"result_id": "rid-applied", "mol_id": "m1", "status": "applied"},
                {"result_id": "rid-dup", "mol_id": "m2", "status": "duplicate"},
            ],
        }
        runner.flush_offline_queue()
        assert runner._offline_queue.pending() == []
        assert _read_dead_letter(queue_path) == []

    def test_poison_does_not_block_others(self, tmp_path: Path) -> None:
        queue_path = tmp_path / "w1_offline_results.jsonl"
        client = MagicMock(spec=WorkerClient)
        runner = WorkerRunner(
            _make_config(offline_queue_path=str(queue_path), worker_id="w1"),
            pipeline=_mock_pipeline(),
            client=client,
        )
        runner._offline_queue.enqueue(
            mol_id="m-poison",
            success=True,
            result_update={"x": 1},
            result_id="rid-poison",
            attempt_id="att-p",
            completed_at="2026-07-12T10:00:00Z",
        )
        runner._offline_queue.enqueue(
            mol_id="m-ok",
            success=True,
            result_update={"x": 2},
            result_id="rid-ok",
            attempt_id="att-ok",
            completed_at="2026-07-12T10:00:01Z",
        )
        client.sync_results.side_effect = lambda results: {
            "accepted": 1,
            "rejected": 1,
            "items": [
                {
                    "result_id": "rid-poison",
                    "mol_id": "m-poison",
                    "status": "conflict",
                    "detail": "poison",
                },
                {"result_id": "rid-ok", "mol_id": "m-ok", "status": "applied"},
            ],
        }
        runner.flush_offline_queue()
        assert runner._offline_queue.pending() == []
        dl = _read_dead_letter(queue_path)
        assert len(dl) == 1
        assert dl[0]["result_id"] == "rid-poison"

    def test_temporary_rejected_keeps_item_in_queue(self, tmp_path: Path) -> None:
        queue_path = tmp_path / "w1_offline_results.jsonl"
        client = MagicMock(spec=WorkerClient)
        runner = WorkerRunner(
            _make_config(offline_queue_path=str(queue_path), worker_id="w1"),
            pipeline=_mock_pipeline(),
            client=client,
        )
        runner._offline_queue.enqueue(
            mol_id="m1",
            success=True,
            result_update={},
            result_id="rid-rej",
            completed_at="2026-07-12T10:00:00Z",
        )
        client.sync_results.side_effect = lambda results: {
            "accepted": 0,
            "rejected": 1,
            "items": [
                {
                    "result_id": "rid-rej",
                    "mol_id": "m1",
                    "status": "rejected",
                    "detail": "temporary",
                }
            ],
        }
        runner.flush_offline_queue()
        pending = runner._offline_queue.pending()
        assert len(pending) == 1
        assert pending[0].result_id == "rid-rej"
        assert _read_dead_letter(queue_path) == []


class TestHeartbeatAttemptWiring:
    @patch(
        "grimperium.worker.runner._pm7result_to_update",
        return_value={"H298_pm7": -42.0},
    )
    def test_runner_heartbeat_passes_attempt_id(self, _mock_update: MagicMock) -> None:
        import threading

        client = _mock_client(claim_returns=("m1", "CCO", "att-hb-1"))
        pipeline = _mock_pipeline(success=True)
        hb_seen: list[tuple[Any, ...]] = []
        release = threading.Event()

        def slow_process(mol_id: str, smiles: str) -> MagicMock:
            release.wait(timeout=2.0)
            return pipeline.process_molecule.return_value

        pipeline.process_molecule.side_effect = slow_process

        def capture_hb(mol_id: str, attempt_id: str | None = None) -> None:
            hb_seen.append((mol_id, attempt_id))

        client.heartbeat.side_effect = capture_hb
        runner = WorkerRunner(
            _make_config(heartbeat_interval_s=0.05),
            pipeline=pipeline,
            client=client,
        )
        thread = threading.Thread(target=runner.run_one, daemon=True)
        thread.start()
        for _ in range(50):
            if hb_seen:
                break
            time.sleep(0.05)
        release.set()
        thread.join(timeout=5.0)
        assert hb_seen, "expected at least one heartbeat"
        assert hb_seen[0] == ("m1", "att-hb-1")
