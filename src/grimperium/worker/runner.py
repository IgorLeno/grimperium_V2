"""Worker main processing loop — claim, process, report."""

import logging
import threading
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from grimperium.crest_pm7 import CRESTPM7Pipeline, PM7Config, PM7Result
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.progress import (
    CREST_STATUS_NOT_ATTEMPTED,
    MOPAC_STATUS_NOT_ATTEMPTED,
)
from grimperium.worker.client import (
    LeaseLostError,
    ServerError,
    WorkerClient,
    WorkerClientConfig,
)
from grimperium.worker.local_store import LocalStore
from grimperium.worker.offline_queue import (
    DeadLetterQueue,
    DeadLetterRecord,
    OfflineResult,
    OfflineResultQueue,
    PendingAbortQueue,
    compute_dead_letter_id,
    dead_letter_path_for,
    new_result_id,
    pending_abort_path_for,
)

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
    attempt_id: str | None = None,
    lease_lost: threading.Event | None = None,
) -> None:
    """Send periodic heartbeats until stop_event is set.

    Erros transitórios (timeout/5xx) apenas logam warning.
    Perda definitiva de lease (409/404) seta ``lease_lost`` e encerra o loop.
    """
    while not stop_event.wait(timeout=interval_s):
        try:
            client.heartbeat(mol_id, attempt_id=attempt_id)
            LOG.debug("Heartbeat sent for %s", mol_id)
        except LeaseLostError as exc:
            LOG.error("Lease lost for %s: %s", mol_id, exc)
            if lease_lost is not None:
                lease_lost.set()
            return
        except Exception:
            LOG.warning("Heartbeat failed for %s (transient)", mol_id)


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
    csv_path: str | None = None
    offline_queue_path: str | None = None


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
        queue_path = (
            Path(config.offline_queue_path)
            if config.offline_queue_path
            else Path.home()
            / ".grimperium"
            / "worker"
            / f"{config.worker_id}_offline_results.jsonl"
        )
        self._offline_queue = OfflineResultQueue(queue_path)
        self._dead_letter = DeadLetterQueue(dead_letter_path_for(queue_path))
        self._pending_aborts = PendingAbortQueue(pending_abort_path_for(queue_path))
        self._stop_event = threading.Event()
        self._consecutive_failures: int = 0
        self._last_run_succeeded: bool = False
        self._csv_manager: BatchCSVManager | None = None
        if config.csv_path:
            self._csv_manager = BatchCSVManager(Path(config.csv_path))

    def _update_csv(self, mol_id: str, updates: dict[str, Any]) -> None:
        """Best-effort local CSV write for progress UI refreshes."""
        if self._csv_manager is None:
            return

        try:
            df = self._csv_manager.load_csv()
            mask = df["mol_id"] == mol_id
            if not mask.any():
                raise KeyError(f"mol_id not found: {mol_id}")
            idx = int(df[mask].index[0])
            for column, value in updates.items():
                if column in df.columns:
                    df.at[idx, column] = value
            self._csv_manager.save_csv()
        except Exception as exc:
            LOG.debug("CSV update skipped for %s: %s", mol_id, exc)

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
                    crest_timeout=(float(self._config.crest_timeout_minutes) * 60.0),
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

        mol_id, smiles, attempt_id = claimed
        self._store.add(mol_id, smiles, attempt_id=attempt_id)
        self._update_csv(
            mol_id,
            {
                "status": MoleculeStatus.RUNNING.value,
                "crest_status": CREST_STATUS_NOT_ATTEMPTED,
                "mopac_status": MOPAC_STATUS_NOT_ATTEMPTED,
            },
        )

        _stop_hb = threading.Event()
        lease_lost = threading.Event()
        hb_thread = threading.Thread(
            target=_run_heartbeat,
            args=(
                mol_id,
                _stop_hb,
                self._client,
                self._config.heartbeat_interval_s,
                attempt_id,
                lease_lost,
            ),
            daemon=True,
        )
        hb_thread.start()

        try:
            result = self._pipeline.process_molecule(mol_id, smiles)
            if lease_lost.is_set():
                # Pipeline químico pode não ser cancelável; não publicar como válido.
                self._archive_aborted_lease_loss(
                    mol_id=mol_id,
                    attempt_id=attempt_id,
                    error="lease lost during processing",
                )
                self._last_run_succeeded = False
                self._consecutive_failures += 1
            elif result.success:
                update = _pm7result_to_update(
                    result,
                    self._config.worker_id,
                    float(self._config.crest_timeout_minutes),
                    float(self._config.mopac_timeout_minutes),
                )
                self._update_csv(
                    mol_id,
                    {
                        **update,
                        "status": MoleculeStatus.OK.value,
                    },
                )
                self._store.mark_success(mol_id, update)
                self._report_or_enqueue(
                    mol_id=mol_id,
                    success=True,
                    result_update=update,
                    error=None,
                )
                self._last_run_succeeded = True
                self._consecutive_failures = 0
            else:
                error = result.error_message or "pipeline failed"
                failure_update: dict[str, Any] = {}
                try:
                    failure_update = _pm7result_to_update(
                        result,
                        self._config.worker_id,
                        float(self._config.crest_timeout_minutes),
                        float(self._config.mopac_timeout_minutes),
                    )
                except Exception as exc:
                    LOG.debug("CSV result mapping skipped for %s: %s", mol_id, exc)
                self._update_csv(
                    mol_id,
                    {
                        **failure_update,
                        "status": MoleculeStatus.RERUN.value,
                    },
                )
                self._store.mark_failure(mol_id, error)
                self._report_or_enqueue(
                    mol_id=mol_id,
                    success=False,
                    result_update=failure_update or None,
                    error=error,
                )
                self._last_run_succeeded = False
                self._consecutive_failures += 1
        except Exception as exc:
            LOG.exception("Unhandled error processing %s", mol_id)
            error_str = str(exc)
            if lease_lost.is_set():
                self._archive_aborted_lease_loss(
                    mol_id=mol_id,
                    attempt_id=attempt_id,
                    error=f"lease lost; also: {error_str}",
                )
            else:
                self._update_csv(mol_id, {"status": MoleculeStatus.RERUN.value})
                self._store.mark_failure(mol_id, error_str)
                try:
                    self._report_or_enqueue(
                        mol_id=mol_id,
                        success=False,
                        result_update=None,
                        error=error_str,
                    )
                except Exception:
                    LOG.exception("Failed to report failure for %s", mol_id)
            self._last_run_succeeded = False
            self._consecutive_failures += 1
        finally:
            _stop_hb.set()
            hb_thread.join(timeout=5.0)

        self._flush_pending_aborts()
        self._store.clear_completed()
        self.flush_offline_queue()
        return True

    def _flush_pending_aborts(self) -> None:
        """Retry durable archive of lease-loss aborts still pending dead-letter."""
        for record in list(self._pending_aborts.entries()):
            dead_letter_id = record.dead_letter_id or compute_dead_letter_id(
                result_id=record.result_id,
                returned_status=record.returned_status,
                original_payload=record.original_payload,
            )
            try:
                self._dead_letter.append(record)
            except Exception:
                LOG.exception(
                    "Pending lease-loss abort still not archived "
                    "(dead_letter_id=%s mol_id=%s)",
                    dead_letter_id,
                    record.mol_id,
                )
                continue
            try:
                self._pending_aborts.remove(dead_letter_id)
            except Exception:
                LOG.exception(
                    "Dead-letter archived but pending remove failed "
                    "(dead_letter_id=%s)",
                    dead_letter_id,
                )
            self._apply_local_lease_loss_side_effects(
                mol_id=record.mol_id,
                result_id=record.result_id,
                error=record.detail or "lease lost",
            )

    def _apply_local_lease_loss_side_effects(
        self,
        *,
        mol_id: str,
        result_id: str,
        error: str,
    ) -> None:
        """Atualizar CSV/local store após evidência de lease-loss estar durável."""
        self._update_csv(mol_id, {"status": MoleculeStatus.RERUN.value})
        record = self._store.get(mol_id)
        if record is not None:
            record.result_id = result_id
            if not record.completed:
                self._store.mark_failure(mol_id, error)

    def _archive_aborted_lease_loss(
        self,
        *,
        mol_id: str,
        attempt_id: str | None,
        error: str,
    ) -> bool:
        """Arquivar processamento abortado por perda de lease (não sync como válido)."""
        result_id = new_result_id()
        payload = {
            "result_id": result_id,
            "mol_id": mol_id,
            "success": False,
            "result_update": None,
            "error": error,
            "completed_at": datetime.now(timezone.utc)
            .isoformat()
            .replace("+00:00", "Z"),
            "attempt_id": attempt_id,
        }
        dl_record = DeadLetterRecord(
            result_id=result_id,
            mol_id=mol_id,
            attempt_id=attempt_id,
            original_payload=payload,
            returned_status="stale_attempt",
            detail=error,
            worker_id=self._config.worker_id,
            rejected_at=datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            rejection_origin="lease_lost",
        )
        dead_letter_id = compute_dead_letter_id(
            result_id=dl_record.result_id,
            returned_status=dl_record.returned_status,
            original_payload=dl_record.original_payload,
        )
        try:
            self._pending_aborts.append(dl_record)
        except Exception:
            LOG.exception("Failed to persist pending lease-loss abort for %s", mol_id)
            return False
        # Evidência pendente durável: refletir Rerun local sem marcar completed
        # (clear_completed não apaga até o dead-letter confirmar).
        self._update_csv(mol_id, {"status": MoleculeStatus.RERUN.value})
        local = self._store.get(mol_id)
        if local is not None:
            local.result_id = result_id
            local.error = error
        try:
            self._dead_letter.append(dl_record)
        except Exception:
            LOG.exception("Failed to archive lease-lost abort for %s", mol_id)
            return False
        try:
            self._pending_aborts.remove(dead_letter_id)
        except Exception:
            LOG.exception(
                "Lease-loss archived but pending remove failed for %s", mol_id
            )
        self._apply_local_lease_loss_side_effects(
            mol_id=mol_id, result_id=result_id, error=error
        )
        return True

    def _report_or_enqueue(
        self,
        *,
        mol_id: str,
        success: bool,
        result_update: dict[str, Any] | None,
        error: str | None,
    ) -> None:
        """Persistir na fila e tentar esvaziar imediatamente via /sync_results."""
        record = self._store.get(mol_id)
        result_id = record.result_id if record is not None else None
        attempt_id = record.attempt_id if record is not None else None
        entry = self._offline_queue.enqueue(
            mol_id=mol_id,
            success=success,
            result_update=result_update,
            error=error,
            result_id=result_id,
            attempt_id=attempt_id,
            completed_at=(
                record.completed_at.isoformat().replace("+00:00", "Z")
                if record is not None and record.completed_at is not None
                else None
            ),
        )
        if record is not None:
            record.result_id = entry.result_id
        # Tentativa online imediata = mesmo protocolo offline (fila persistida).
        self._flush_entries([entry])

    def flush_offline_queue(self) -> tuple[int, int]:
        """Reenviar resultados pendentes via /sync_results com o mesmo result_id."""
        pending = self._offline_queue.pending()
        if not pending:
            return 0, 0
        return self._flush_entries(pending)

    def _flush_entries(self, pending: list[OfflineResult]) -> tuple[int, int]:
        """Confirmar apenas itens applied/duplicate; conflict/stale → dead-letter."""
        try:
            response = self._client.sync_results(
                [entry.to_sync_dict() for entry in pending]
            )
        except ServerError as exc:
            LOG.warning("Offline queue flush failed: %s", exc)
            return 0, len(pending)

        items = response.get("items")
        confirmed: set[str] = set()
        by_result_id = {entry.result_id: entry for entry in pending}
        dead_letter_statuses = {"conflict", "stale_attempt"}
        confirm_statuses = {"applied", "duplicate"} | dead_letter_statuses
        if isinstance(items, list) and items:
            for item in items:
                if not isinstance(item, dict):
                    continue
                status = str(item.get("status", ""))
                result_id = str(item.get("result_id", ""))
                if not result_id or status not in confirm_statuses:
                    continue
                if status in dead_letter_statuses:
                    entry = by_result_id.get(result_id)
                    if entry is None:
                        LOG.error(
                            "Terminal sync rejection for unknown result_id=%s "
                            "status=%s — keeping offline queue unchanged",
                            result_id,
                            status,
                        )
                        continue
                    try:
                        self._dead_letter.append(
                            DeadLetterRecord(
                                result_id=result_id,
                                mol_id=entry.mol_id,
                                attempt_id=entry.attempt_id,
                                original_payload=entry.to_sync_dict(),
                                returned_status=status,
                                detail=(
                                    str(item["detail"])
                                    if item.get("detail") is not None
                                    else None
                                ),
                                worker_id=self._config.worker_id,
                                rejected_at=datetime.now(timezone.utc)
                                .isoformat()
                                .replace("+00:00", "Z"),
                                rejection_origin="sync_results",
                            )
                        )
                    except Exception:
                        LOG.exception(
                            "Dead-letter write failed for result_id=%s status=%s — "
                            "keeping item in offline queue",
                            result_id,
                            status,
                        )
                        continue
                    LOG.error(
                        "Terminal sync rejection for result_id=%s status=%s — "
                        "moved to dead-letter",
                        result_id,
                        status,
                    )
                confirmed.add(result_id)
                self._offline_queue.confirm(result_id)
        else:
            # Compat com servidores sem items: só confirmar se rejected==0.
            rejected = int(response.get("rejected", 0))
            if rejected == 0:
                for entry in pending:
                    self._offline_queue.confirm(entry.result_id)
                    confirmed.add(entry.result_id)

        accepted = int(response.get("accepted", len(confirmed)))
        rejected = int(response.get("rejected", max(0, len(pending) - len(confirmed))))
        return accepted, rejected

    def run(self, max_molecules: int | None = None) -> int:
        """Process molecules until stopped, queue exhausted, or max reached.

        Args:
            max_molecules: Stop after processing this many. None means
                unlimited.

        Returns:
            Number of molecules successfully processed.
        """
        self._client.register()
        self._flush_pending_aborts()
        self.flush_offline_queue()
        processed = 0
        attempted = 0
        idle_count = 0

        while not self._stop_event.is_set():
            if max_molecules is not None and attempted >= max_molecules:
                break
            did_work = self.run_one()
            if did_work:
                attempted += 1
                if self._last_run_succeeded:
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
