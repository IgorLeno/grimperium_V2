"""Worker main processing loop — claim, process, report."""

import logging
import threading
import time
from dataclasses import dataclass
from typing import Any

from grimperium.crest_pm7 import CRESTPM7Pipeline, PM7Config, PM7Result
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.worker.client import WorkerClient, WorkerClientConfig
from grimperium.worker.local_store import LocalStore

LOG = logging.getLogger(__name__)

# Single shared instance — pm7result_to_csv_update uses only class-level data.
_RESULT_MAPPER = BatchCSVManager(csv_path=None)


def _pm7result_to_update(
    result: PM7Result,
    worker_id: str,
    crest_timeout_minutes: float,
    mopac_timeout_minutes: float,
) -> dict[str, Any]:
    return _RESULT_MAPPER.pm7result_to_csv_update(
        mol_id=result.mol_id,
        result=result,
        batch_id=worker_id,
        batch_order=0,
        crest_timeout_used=crest_timeout_minutes,
        mopac_timeout_used=mopac_timeout_minutes,
    )


def _run_heartbeat(
    mol_id: str,
    stop_event: threading.Event,
    client: WorkerClient,
    interval_s: float,
) -> None:
    """Send periodic heartbeats until stop_event is set."""
    while not stop_event.wait(timeout=interval_s):
        try:
            client.heartbeat(mol_id)
            LOG.debug("Heartbeat sent for %s", mol_id)
        except Exception:
            LOG.warning("Heartbeat failed for %s", mol_id)


@dataclass
class WorkerConfig:
    """Configuration for a single worker process."""

    server_url: str
    worker_id: str
    api_token: str = ""
    heartbeat_interval_s: float = 30.0
    poll_interval_s: float = 5.0
    max_idle_polls: int = 12
    crest_timeout_minutes: int = 30
    mopac_timeout_minutes: int = 10
    batch_size: int = 10
    batch_id: str = "worker"
    consecutive_failure_stop: bool = True
    max_consecutive_failures: int = 10


class WorkerRunner:
    """Single-worker processing loop.

    Typical usage::

        runner = WorkerRunner(config)
        runner.run()          # blocks until queue empty or stopped

    Inject ``pipeline`` and ``client`` for testing.
    """

    def __init__(
        self,
        config: WorkerConfig,
        pipeline: CRESTPM7Pipeline | None = None,
        client: WorkerClient | None = None,
    ) -> None:
        self._config = config
        self._pipeline = pipeline or CRESTPM7Pipeline(
            PM7Config(
                crest_timeout=float(config.crest_timeout_minutes) * 60.0,
                mopac_timeout_base=float(config.mopac_timeout_minutes) * 60.0,
            )
        )
        client_cfg = WorkerClientConfig(
            server_url=config.server_url,
            worker_id=config.worker_id,
            api_token=config.api_token,
        )
        self._client = client or WorkerClient(client_cfg)
        self._store = LocalStore()
        self._stop_event = threading.Event()
        self._consecutive_failures: int = 0

    def stop(self) -> None:
        """Signal run() to exit after the current molecule finishes."""
        self._stop_event.set()

    def reconfigure(self, new_config: dict[str, Any]) -> None:
        """Apply a config override received from the server.

        Only keys present in new_config are updated. If crest_timeout_minutes
        or mopac_timeout_minutes change, the pipeline is rebuilt with the new
        values.

        Args:
            new_config: Partial config dict, e.g. from
                GET /configure/{worker_id}.
        """
        if not new_config:
            return

        timeouts_changed = False
        if "crest_timeout_minutes" in new_config:
            val = int(new_config["crest_timeout_minutes"])
            if val != self._config.crest_timeout_minutes:
                self._config.crest_timeout_minutes = val
                timeouts_changed = True
        if "mopac_timeout_minutes" in new_config:
            val = int(new_config["mopac_timeout_minutes"])
            if val != self._config.mopac_timeout_minutes:
                self._config.mopac_timeout_minutes = val
                timeouts_changed = True
        if "batch_size" in new_config:
            self._config.batch_size = int(new_config["batch_size"])

        if timeouts_changed:
            self._pipeline = CRESTPM7Pipeline(
                PM7Config(
                    crest_timeout=(
                        float(self._config.crest_timeout_minutes) * 60.0
                    ),
                    mopac_timeout_base=(
                        float(self._config.mopac_timeout_minutes) * 60.0
                    ),
                )
            )
            LOG.info(
                "Pipeline rebuilt with crest_timeout=%dm mopac_timeout=%dm",
                self._config.crest_timeout_minutes,
                self._config.mopac_timeout_minutes,
            )

    def run_one(self) -> bool:
        """Claim and process one molecule.

        Returns True if a molecule was processed, False if queue was empty.
        """
        claimed = self._client.claim()
        if claimed is None:
            return False

        mol_id, smiles = claimed
        self._store.add(mol_id, smiles)

        _stop_hb = threading.Event()
        hb_thread = threading.Thread(
            target=_run_heartbeat,
            args=(
                mol_id,
                _stop_hb,
                self._client,
                self._config.heartbeat_interval_s,
            ),
            daemon=True,
        )
        hb_thread.start()

        try:
            result = self._pipeline.process_molecule(mol_id, smiles)
            if result.success:
                update = _pm7result_to_update(
                    result,
                    self._config.worker_id,
                    float(self._config.crest_timeout_minutes),
                    float(self._config.mopac_timeout_minutes),
                )
                self._store.mark_success(mol_id, update)
                self._client.report_success(mol_id, update)
                self._consecutive_failures = 0
            else:
                error = result.error_message or "pipeline failed"
                self._store.mark_failure(mol_id, error)
                self._client.report_failure(mol_id, error)
                self._consecutive_failures += 1
        except Exception as exc:
            LOG.exception("Unhandled error processing %s", mol_id)
            error_str = str(exc)
            self._store.mark_failure(mol_id, error_str)
            self._consecutive_failures += 1
            try:
                self._client.report_failure(mol_id, error_str)
            except Exception:
                LOG.exception("Failed to report failure for %s", mol_id)
        finally:
            _stop_hb.set()
            hb_thread.join(timeout=5.0)

        self._store.clear_completed()
        return True

    def run(self, max_molecules: int | None = None) -> int:
        """Process molecules until stopped, queue exhausted, or max reached.

        Args:
            max_molecules: Stop after processing this many. None means
                unlimited.

        Returns:
            Number of molecules successfully processed.
        """
        self._client.register()
        processed = 0
        attempted = 0
        idle_count = 0

        while not self._stop_event.is_set():
            if max_molecules is not None and attempted >= max_molecules:
                break
            did_work = self.run_one()
            if did_work:
                attempted += 1
                if self._consecutive_failures == 0:
                    processed += 1
                idle_count = 0
                if (
                    self._config.consecutive_failure_stop
                    and self._consecutive_failures
                    >= self._config.max_consecutive_failures
                ):
                    LOG.warning(
                        "Stopping: %d consecutive failures (max=%d)",
                        self._consecutive_failures,
                        self._config.max_consecutive_failures,
                    )
                    break
            else:
                idle_count += 1
                if idle_count >= self._config.max_idle_polls:
                    LOG.info(
                        "No molecules for %d consecutive polls — stopping",
                        idle_count,
                    )
                    break
                time.sleep(self._config.poll_interval_s)

        return processed
