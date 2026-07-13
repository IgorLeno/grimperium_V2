"""Serviço transacional único para entrega de resultados do worker.

Online e offline compartilham o mesmo fluxo:
check → prepare → dual-write → commit → métricas → cleanup.

Métricas do ``WorkerRegistry`` são best-effort / observabilidade em memória:
não entram no journal e um retry ``duplicate`` não as reconcilia. O estado
científico e operacional (CSV + ledger) permanece a autoridade transacional.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from grimperium.crest_pm7.batch.result_applier import BatchResultApplier
from grimperium.crest_pm7.batch.result_ledger import (
    JournalEntry,
    JournalTxnStatus,
    LedgerStatus,
    ResultLedger,
    build_legacy_result_id,
    build_result_fingerprint,
)
from grimperium.crest_pm7.batch.state_manager import BatchStateManager
from grimperium.server.models import SyncItemOutcome, SyncItemResult, SyncResult
from grimperium.server.worker_registry import WorkerRegistry

LOG = logging.getLogger(__name__)


class SyncConflictError(Exception):
    """result_id reutilizado com fingerprint diferente."""

    def __init__(self, result_id: str) -> None:
        self.result_id = result_id
        super().__init__(f"Conflicting sync result_id: {result_id}")


@dataclass(frozen=True)
class SyncApplyOutcome:
    """Resultado de processar um único SyncResult."""

    item: SyncItemResult
    final_status: str | None = None


class SyncResultApplicationService:
    """Centraliza idempotência, journal e dual-write de resultados."""

    def __init__(
        self,
        *,
        csv_manager: BatchCSVManager,
        state_manager: BatchStateManager,
        ledger: ResultLedger,
        worker_registry: WorkerRegistry,
        running_molecules: dict[str, str],
    ) -> None:
        self._csv_manager = csv_manager
        self._state_manager = state_manager
        self._ledger = ledger
        self._worker_registry = worker_registry
        self._running_molecules = running_molecules

    def apply_one(self, worker_id: str, result: SyncResult) -> SyncApplyOutcome:
        """Aplicar um resultado com efeito operacional no máximo uma vez."""
        result_id = self.resolve_result_id(result)
        payload = sync_result_payload(result)
        fingerprint = build_result_fingerprint(payload)

        decision = self._ledger.check(result_id, fingerprint)
        if decision.status is LedgerStatus.DUPLICATE:
            self._cleanup_assignment(worker_id, result.mol_id, result.attempt_id)
            return SyncApplyOutcome(
                item=SyncItemResult(
                    result_id=result_id,
                    mol_id=result.mol_id,
                    status=SyncItemOutcome.DUPLICATE,
                )
            )
        if decision.status is LedgerStatus.CONFLICT:
            raise SyncConflictError(result_id)

        incomplete = {entry.result_id: entry for entry in self._ledger.get_incomplete()}
        prepared = incomplete.get(result_id)
        if prepared is not None and prepared.fingerprint == fingerprint:
            if self.verify_applied(prepared):
                self._ledger.commit(
                    result_id,
                    final_status=prepared.expected_final_status
                    or self._state_manager.get_status(result.mol_id),
                )
                self._cleanup_assignment(worker_id, result.mol_id, result.attempt_id)
                return SyncApplyOutcome(
                    item=SyncItemResult(
                        result_id=result_id,
                        mol_id=result.mol_id,
                        status=SyncItemOutcome.DUPLICATE,
                        detail="recovered prepared transaction",
                    )
                )
            current_status = self._state_manager.get_status(result.mol_id)
            if str(current_status).lower() in {"running", "claimed", "selected"}:
                return self._resume_prepared(
                    worker_id=worker_id,
                    result=result,
                    result_id=result_id,
                )
            # Prepare sem apply + reclaim para Pending: retomar com segurança.
            prev = (prepared.previous_status or "").lower()
            if str(current_status).lower() == "pending" and prev in {
                "running",
                "claimed",
                "selected",
            }:
                return self._resume_prepared(
                    worker_id=worker_id,
                    result=result,
                    result_id=result_id,
                )
            LOG.warning(
                "Prepared journal without exact proof; refusing blind reapply "
                "(result_id=%s mol_id=%s status=%s)",
                result_id,
                result.mol_id,
                current_status,
            )
            return SyncApplyOutcome(
                item=SyncItemResult(
                    result_id=result_id,
                    mol_id=result.mol_id,
                    status=SyncItemOutcome.REJECTED,
                    detail="prepared without verifiable effect",
                )
            )

        stale = self._stale_attempt_outcome(result_id, result)
        if stale is not None:
            return stale

        previous_status = self._state_manager.get_status(result.mol_id)
        previous_reruns = self._state_manager.get_reruns(result.mol_id)
        expected = self._expected_effect(
            success=result.success,
            previous_status=previous_status,
            previous_reruns=previous_reruns,
            force_skip=False,
            result_update=result.result_update,
        )
        self._ledger.prepare(
            result_id=result_id,
            mol_id=result.mol_id,
            fingerprint=fingerprint,
            desired_success=result.success,
            previous_status=previous_status,
            previous_reruns=previous_reruns,
            expected_final_status=expected.final_status,
            expected_reruns=expected.reruns,
            expected_science_hash=expected.science_hash,
            worker_id=worker_id,
            attempt_id=result.attempt_id,
        )
        try:
            final_status = apply_worker_result(
                self._csv_manager,
                self._state_manager,
                mol_id=result.mol_id,
                success=result.success,
                error=result.error,
                result_update=result.result_update,
            )
        except Exception as apply_exc:
            self._ledger.mark_failed(result_id, error=str(apply_exc))
            raise

        self._ledger.commit(result_id, final_status=final_status)
        self._record_metrics(worker_id, final_status=final_status)
        self._cleanup_assignment(worker_id, result.mol_id, result.attempt_id)
        return SyncApplyOutcome(
            item=SyncItemResult(
                result_id=result_id,
                mol_id=result.mol_id,
                status=SyncItemOutcome.APPLIED,
            ),
            final_status=final_status,
        )

    def _stale_attempt_outcome(
        self, result_id: str, result: SyncResult
    ) -> SyncApplyOutcome | None:
        """Reject results that do not match the current assignment lease."""
        current_attempt = self._state_manager.get_attempt_id(result.mol_id)
        if result.attempt_id is None:
            if current_attempt:
                return SyncApplyOutcome(
                    item=SyncItemResult(
                        result_id=result_id,
                        mol_id=result.mol_id,
                        status=SyncItemOutcome.STALE_ATTEMPT,
                        detail="missing attempt_id for active assignment",
                    )
                )
            # CSV/payload legado sem lease: aceitar como antes.
            return None
        if current_attempt is None:
            try:
                status = str(self._state_manager.get_status(result.mol_id)).lower()
            except Exception:
                status = ""
            # Molécula terminal: deferir ao ledger (duplicate/recovery).
            if status in {"ok", "skip"}:
                return None
            return SyncApplyOutcome(
                item=SyncItemResult(
                    result_id=result_id,
                    mol_id=result.mol_id,
                    status=SyncItemOutcome.STALE_ATTEMPT,
                    detail=(
                        f"attempt_id={result.attempt_id!r} no longer assigned "
                        f"(status={status!r})"
                    ),
                )
            )
        if result.attempt_id != current_attempt:
            return SyncApplyOutcome(
                item=SyncItemResult(
                    result_id=result_id,
                    mol_id=result.mol_id,
                    status=SyncItemOutcome.STALE_ATTEMPT,
                    detail=(
                        f"attempt_id mismatch: got={result.attempt_id!r} "
                        f"current={current_attempt!r}"
                    ),
                )
            )
        return None

    def _resume_prepared(
        self,
        *,
        worker_id: str,
        result: SyncResult,
        result_id: str,
    ) -> SyncApplyOutcome:
        """Retomar dual-write de uma transação prepared ainda em voo."""
        try:
            final_status = apply_worker_result(
                self._csv_manager,
                self._state_manager,
                mol_id=result.mol_id,
                success=result.success,
                error=result.error,
                result_update=result.result_update,
            )
        except Exception as apply_exc:
            self._ledger.mark_failed(result_id, error=str(apply_exc))
            raise
        self._ledger.commit(result_id, final_status=final_status)
        self._record_metrics(worker_id, final_status=final_status)
        self._cleanup_assignment(worker_id, result.mol_id, result.attempt_id)
        return SyncApplyOutcome(
            item=SyncItemResult(
                result_id=result_id,
                mol_id=result.mol_id,
                status=SyncItemOutcome.APPLIED,
                detail="resumed prepared transaction",
            ),
            final_status=final_status,
        )

    def apply_force_skip(
        self, worker_id: str, result: SyncResult, error: str
    ) -> SyncApplyOutcome:
        """Aplicar force_skip legado pelo mesmo contrato de ledger."""
        result_id = self.resolve_result_id(result)
        payload = sync_result_payload(result)
        fingerprint = build_result_fingerprint(payload)
        decision = self._ledger.check(result_id, fingerprint)
        if decision.status is LedgerStatus.DUPLICATE:
            self._cleanup_assignment(worker_id, result.mol_id, result.attempt_id)
            return SyncApplyOutcome(
                item=SyncItemResult(
                    result_id=result_id,
                    mol_id=result.mol_id,
                    status=SyncItemOutcome.DUPLICATE,
                )
            )
        if decision.status is LedgerStatus.CONFLICT:
            raise SyncConflictError(result_id)

        stale = self._stale_attempt_outcome(result_id, result)
        if stale is not None:
            return stale

        previous_status = self._state_manager.get_status(result.mol_id)
        previous_reruns = self._state_manager.get_reruns(result.mol_id)
        self._ledger.prepare(
            result_id=result_id,
            mol_id=result.mol_id,
            fingerprint=fingerprint,
            desired_success=False,
            previous_status=previous_status,
            previous_reruns=previous_reruns,
            expected_final_status=MoleculeStatus.SKIP.value,
            expected_reruns=previous_reruns + 1,
            expected_science_hash=None,
            worker_id=worker_id,
            attempt_id=result.attempt_id,
        )
        try:
            final_status = apply_worker_result(
                self._csv_manager,
                self._state_manager,
                mol_id=result.mol_id,
                success=False,
                error=error,
                force_skip=True,
            )
        except Exception as apply_exc:
            self._ledger.mark_failed(result_id, error=str(apply_exc))
            raise
        self._ledger.commit(result_id, final_status=final_status)
        self._record_metrics(worker_id, final_status=final_status)
        self._cleanup_assignment(worker_id, result.mol_id, result.attempt_id)
        return SyncApplyOutcome(
            item=SyncItemResult(
                result_id=result_id,
                mol_id=result.mol_id,
                status=SyncItemOutcome.APPLIED,
            ),
            final_status=final_status,
        )

    def resolve_result_id(self, result: SyncResult) -> str:
        """Resolver result_id explícito ou fallback legado estável."""
        if result.result_id:
            return result.result_id
        return build_legacy_result_id(sync_result_payload(result))

    def verify_applied(self, entry: JournalEntry) -> bool:
        """Provar que o efeito operacional/científico já está no estado persistido."""
        try:
            current_status = self._state_manager.get_status(entry.mol_id)
            current_reruns = self._state_manager.get_reruns(entry.mol_id)
        except Exception:
            return False

        if str(current_status).lower() in {"running", "claimed", "selected", "pending"}:
            return False

        expected_status = entry.expected_final_status
        if expected_status is not None:
            if current_status != expected_status:
                return False
            if (
                entry.expected_reruns is not None
                and current_reruns != entry.expected_reruns
            ):
                return False
            if entry.desired_success:
                return self._verify_success_science(entry)
            return True

        # Journal legado sem expectativa explícita: prova mínima determinística.
        if entry.desired_success:
            if current_status != MoleculeStatus.OK.value:
                return False
            return self._verify_success_science(entry)

        if current_status not in {
            MoleculeStatus.RERUN.value,
            MoleculeStatus.SKIP.value,
        }:
            return False
        if entry.previous_reruns is None:
            return False
        return current_reruns == entry.previous_reruns + 1

    def _verify_success_science(self, entry: JournalEntry) -> bool:
        """Confirmar campos científicos essenciais quando o journal registra hash."""
        if not entry.expected_science_hash:
            return True
        if self._csv_manager.df is None:
            return False
        try:
            row = self._csv_manager.df.loc[
                self._csv_manager.df["mol_id"] == entry.mol_id
            ].iloc[0]
        except (IndexError, KeyError):
            return False
        science = {
            "H298_pm7": _nullable_float(row.get("H298_pm7")),
            "G298_pm7": _nullable_float(row.get("G298_pm7")),
            "gap": _nullable_float(row.get("gap")),
        }
        return build_result_fingerprint(science) == entry.expected_science_hash

    def _expected_effect(
        self,
        *,
        success: bool,
        previous_status: str,
        previous_reruns: int,
        force_skip: bool,
        result_update: dict[str, Any] | None,
    ) -> _ExpectedEffect:
        del previous_status  # usado apenas para auditoria no journal
        if success:
            science = {
                "H298_pm7": (result_update or {}).get("H298_pm7"),
                "G298_pm7": (result_update or {}).get("G298_pm7"),
                "gap": (result_update or {}).get("gap"),
            }
            return _ExpectedEffect(
                final_status=MoleculeStatus.OK.value,
                reruns=previous_reruns,
                science_hash=build_result_fingerprint(science),
            )
        if force_skip:
            return _ExpectedEffect(
                final_status=MoleculeStatus.SKIP.value,
                reruns=previous_reruns + 1,
                science_hash=None,
            )
        max_reruns = self._state_manager.max_reruns
        next_reruns = previous_reruns + 1
        if next_reruns >= max_reruns:
            return _ExpectedEffect(
                final_status=MoleculeStatus.SKIP.value,
                reruns=next_reruns,
                science_hash=None,
            )
        return _ExpectedEffect(
            final_status=MoleculeStatus.RERUN.value,
            reruns=next_reruns,
            science_hash=None,
        )

    def _cleanup_assignment(
        self,
        worker_id: str,
        mol_id: str,
        attempt_id: str | None,
    ) -> None:
        """Clear assignment only when the result still matches the current lease."""
        current_attempt = self._state_manager.get_attempt_id(mol_id)
        if current_attempt and attempt_id and attempt_id != current_attempt:
            return
        self._running_molecules.pop(mol_id, None)
        self._worker_registry.clear_current_mol(worker_id)

    def _record_metrics(self, worker_id: str, *, final_status: str) -> None:
        # Best-effort: registry is in-memory observability, not transactional.
        if final_status == MoleculeStatus.OK.value:
            self._worker_registry.record_success(worker_id)
        elif final_status == MoleculeStatus.SKIP.value:
            self._worker_registry.record_skip(worker_id)
        else:
            self._worker_registry.record_failure(worker_id)


@dataclass(frozen=True)
class _ExpectedEffect:
    final_status: str
    reruns: int
    science_hash: str | None


def sync_result_payload(result: SyncResult) -> dict[str, Any]:
    """Payload estável usado no fingerprint do ledger."""
    payload: dict[str, Any] = {
        "mol_id": result.mol_id,
        "success": result.success,
        "result_update": result.result_update,
        "error": result.error,
        "completed_at": result.completed_at,
    }
    if result.attempt_id is not None:
        payload["attempt_id"] = result.attempt_id
    return payload


def apply_worker_result(
    csv_manager: BatchCSVManager,
    state_manager: BatchStateManager,
    *,
    mol_id: str,
    success: bool,
    error: str | None = None,
    force_skip: bool = False,
    result_update: dict[str, Any] | None = None,
) -> str:
    """Aplicar um resultado via dual-write compartilhado."""
    applier = BatchResultApplier(state_manager=state_manager, csv_manager=csv_manager)
    if success:
        return applier.apply_success(mol_id, result_update or {}).final_status
    return applier.apply_failure(
        mol_id,
        error or "failure",
        force_skip=force_skip,
        result_update=result_update,
    ).final_status


def recover_sync_journal(
    ledger: ResultLedger,
    service: SyncResultApplicationService,
) -> list[JournalEntry]:
    """Commit apenas quando a prova de efeito for exata."""

    def _already_applied(entry: JournalEntry) -> bool:
        return service.verify_applied(entry)

    recovered = ledger.recover_incomplete(is_already_applied=_already_applied)
    committed = [
        entry for entry in recovered if entry.txn_status is JournalTxnStatus.COMMITTED
    ]
    if committed:
        LOG.info(
            "Recovered %d incomplete sync journal entr(y/ies)",
            len(committed),
        )
    stuck = [
        entry for entry in recovered if entry.txn_status is JournalTxnStatus.PREPARED
    ]
    for entry in stuck:
        LOG.warning(
            "Leaving prepared journal entry without exact proof: result_id=%s mol_id=%s",
            entry.result_id,
            entry.mol_id,
        )
    return recovered


def _nullable_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        if str(value).strip() == "" or str(value).lower() == "nan":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None
