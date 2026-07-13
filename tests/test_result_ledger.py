from pathlib import Path

import pytest

from grimperium.crest_pm7.batch.result_ledger import (
    LedgerStatus,
    ResultLedger,
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
