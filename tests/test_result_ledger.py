import json
from pathlib import Path

import pytest

from grimperium.crest_pm7.batch.result_ledger import (
    JournalTxnStatus,
    LedgerStatus,
    OperationKind,
    ResultLedger,
    build_result_fingerprint,
    resolve_compatible_fingerprint,
    with_operation_kind,
)


def test_result_ledger_records_first_result_as_applied(tmp_path: Path) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")

    decision = ledger.check_and_record("worker-1|mol_001|0", "mol_001", "fingerprint-a")

    assert decision.status is LedgerStatus.APPLIED
    assert decision.duplicate is False
    assert (tmp_path / "result_ledger.jsonl").exists()


def test_result_ledger_reports_duplicate_for_same_fingerprint(tmp_path: Path) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    ledger.check_and_record("result-1", "mol_001", "fingerprint-a")

    decision = ledger.check_and_record("result-1", "mol_001", "fingerprint-a")

    assert decision.status is LedgerStatus.DUPLICATE
    assert decision.duplicate is True


def test_result_ledger_reports_conflict_for_reused_id_with_new_fingerprint(
    tmp_path: Path,
) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    ledger.check_and_record("result-1", "mol_001", "fingerprint-a")

    decision = ledger.check_and_record("result-1", "mol_001", "fingerprint-b")

    assert decision.status is LedgerStatus.CONFLICT
    assert decision.conflict is True


def test_result_ledger_rejects_blank_result_id(tmp_path: Path) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")

    with pytest.raises(ValueError, match="result_id"):
        ledger.check_and_record("", "mol_001", "fingerprint-a")


def test_check_conflicts_when_prepared_fingerprint_differs(tmp_path: Path) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    ledger.prepare(
        result_id="result-1",
        mol_id="mol_001",
        fingerprint="fingerprint-a",
        desired_success=True,
    )

    decision = ledger.check("result-1", "fingerprint-b")

    assert decision.status is LedgerStatus.CONFLICT
    assert decision.conflict is True


def test_check_allows_retry_when_failed_fingerprint_matches(tmp_path: Path) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    ledger.prepare(
        result_id="result-1",
        mol_id="mol_001",
        fingerprint="fingerprint-a",
        desired_success=True,
    )
    ledger.mark_failed("result-1", error="dual-write boom")

    decision = ledger.check("result-1", "fingerprint-a")
    assert decision.status is LedgerStatus.APPLIED

    retried = ledger.prepare(
        result_id="result-1",
        mol_id="mol_001",
        fingerprint="fingerprint-a",
        desired_success=True,
    )
    assert retried.txn_status.value == "prepared"


def test_check_conflicts_when_failed_fingerprint_differs(tmp_path: Path) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    ledger.prepare(
        result_id="result-1",
        mol_id="mol_001",
        fingerprint="fingerprint-a",
        desired_success=True,
    )
    ledger.mark_failed("result-1", error="dual-write boom")

    decision = ledger.check("result-1", "fingerprint-b")

    assert decision.status is LedgerStatus.CONFLICT


def test_prepare_rejects_fingerprint_mismatch_while_prepared(tmp_path: Path) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    ledger.prepare(
        result_id="result-1",
        mol_id="mol_001",
        fingerprint="fingerprint-a",
        desired_success=True,
    )

    with pytest.raises(ValueError, match="fingerprint mismatch"):
        ledger.prepare(
            result_id="result-1",
            mol_id="mol_001",
            fingerprint="fingerprint-b",
            desired_success=True,
        )


def test_find_committed_by_attempt_returns_committed_entry(tmp_path: Path) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    ledger.prepare(
        result_id="result-1",
        mol_id="mol_001",
        fingerprint="fingerprint-a",
        desired_success=True,
        attempt_id="attempt-a",
    )
    ledger.commit("result-1", final_status="OK")

    found = ledger.find_committed_by_attempt("attempt-a")

    assert found is not None
    assert found.result_id == "result-1"
    assert found.txn_status is JournalTxnStatus.COMMITTED
    assert found.attempt_id == "attempt-a"


def test_find_committed_by_attempt_ignores_prepared_and_other_attempts(
    tmp_path: Path,
) -> None:
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    ledger.prepare(
        result_id="result-prepared",
        mol_id="mol_001",
        fingerprint="fingerprint-a",
        desired_success=True,
        attempt_id="attempt-a",
    )
    ledger.prepare(
        result_id="result-other",
        mol_id="mol_002",
        fingerprint="fingerprint-b",
        desired_success=True,
        attempt_id="attempt-b",
    )
    ledger.commit("result-other", final_status="OK")

    assert ledger.find_committed_by_attempt("attempt-a") is None
    assert ledger.find_committed_by_attempt("attempt-b") is not None
    assert ledger.find_committed_by_attempt("") is None


def test_operation_kind_changes_fingerprint(tmp_path: Path) -> None:
    base = {"mol_id": "m1", "success": False, "error": "x"}
    normal = build_result_fingerprint(
        with_operation_kind(base, OperationKind.NORMAL_RESULT)
    )
    force = build_result_fingerprint(
        with_operation_kind(base, OperationKind.FORCE_SKIP)
    )
    assert normal != force


def test_journal_preserves_operation_kind_and_legacy_default(tmp_path: Path) -> None:
    path = tmp_path / "result_ledger.jsonl"
    ledger = ResultLedger(path)
    ledger.prepare(
        result_id="r-kind",
        mol_id="mol_001",
        fingerprint="fp-force",
        desired_success=False,
        operation_kind=OperationKind.FORCE_SKIP,
    )
    reloaded = ResultLedger(path)
    entry = reloaded.get_incomplete()[0]
    assert entry.operation_kind == OperationKind.FORCE_SKIP.value

    # Legacy journal line without operation_kind → NORMAL_RESULT
    journal = path.with_name("result_ledger_journal.jsonl")
    legacy_line = {
        "schema_version": 1,
        "result_id": "legacy-r",
        "fingerprint": "fp-legacy",
        "mol_id": "mol_002",
        "txn_status": JournalTxnStatus.PREPARED.value,
        "desired_success": True,
    }
    with journal.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(legacy_line) + "\n")
    reloaded2 = ResultLedger(path)
    legacy = reloaded2.get_incomplete()
    by_id = {e.result_id: e for e in legacy}
    assert by_id["legacy-r"].operation_kind == OperationKind.NORMAL_RESULT.value


def test_resolve_compatible_fingerprint_accepts_pre_upgrade_attempt_id(
    tmp_path: Path,
) -> None:
    """Journals com attempt_id no hash ainda fazem duplicate no payload novo."""
    ledger = ResultLedger(tmp_path / "result_ledger.jsonl")
    payload = {
        "mol_id": "mol_001",
        "success": False,
        "result_update": None,
        "error": "timeout",
        "completed_at": "2026-07-13T00:00:00Z",
    }
    legacy_fp = build_result_fingerprint({**payload, "attempt_id": "att-old"})
    ledger.prepare(
        result_id="compat-r1",
        mol_id="mol_001",
        fingerprint=legacy_fp,
        desired_success=False,
        attempt_id="att-old",
    )
    ledger.commit("compat-r1", final_status="Rerun")

    fingerprint, decision = resolve_compatible_fingerprint(
        ledger,
        result_id="compat-r1",
        payload=payload,
        operation_kind=OperationKind.NORMAL_RESULT,
        attempt_id=None,
    )
    assert decision.status is LedgerStatus.DUPLICATE
    assert fingerprint == legacy_fp
