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
