"""Tests for the DatabaseAnalyzer module."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from grimperium.crest_pm7.database_analyzer import DatabaseAnalyzer


def _write_csv(path: Path, df: pd.DataFrame) -> Path:
    """Helper to write a DataFrame to a CSV file."""
    df.to_csv(path, index=False)
    return path


@pytest.fixture()
def analyzer() -> DatabaseAnalyzer:
    return DatabaseAnalyzer()


# ── Test 1: status_counts ──


def test_status_counts(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """Mix of OK/PENDING/RUNNING → correct counts."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "PENDING", "RUNNING", "OK"],
            "mol_id": ["m1", "m2", "m3", "m4", "m5"],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.total_rows == 5
    assert report.ok_count == 3
    assert report.status_counts["OK"] == 3
    assert report.status_counts["PENDING"] == 1
    assert report.status_counts["RUNNING"] == 1


# ── Test 2: orphan_running ──


def test_orphan_running(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """RUNNING rows without timestamp → detected as orphan."""
    df = pd.DataFrame(
        {
            "status": ["RUNNING", "RUNNING", "OK"],
            "mol_id": ["m1", "m2", "m3"],
            "timestamp": [None, "2026-01-01", None],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.orphan_running == 1  # m1 is RUNNING without timestamp


# ── Test 3: grade_ab_rate ──


def test_grade_ab_rate(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """Grade distribution from mock JSON files, alert if < 67%."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "OK", "OK"],
            "mol_id": ["m1", "m2", "m3", "m4"],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)

    # Create conformer details with grades
    detail_dir = tmp_path / "details"
    detail_dir.mkdir()
    grades = {"m1": "A", "m2": "B", "m3": "C", "m4": "FAILED"}
    for mol_id, grade in grades.items():
        (detail_dir / f"{mol_id}.json").write_text(
            json.dumps({"mol_id": mol_id, "quality_grade": grade}),
            encoding="utf-8",
        )

    report = analyzer.analyze(csv, detail_dir=detail_dir)

    assert report.grade_distribution == {"A": 1, "B": 1, "C": 1, "FAILED": 1}
    assert report.grade_ab_rate == pytest.approx(0.5)
    # A+B = 50% < 67% threshold → alert
    assert any("below 67%" in a for a in report.alerts)


# ── Test 4: null_h298_ok ──


def test_null_h298_ok(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """OK rows with NaN H298_pm7 → alert generated."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "OK"],
            "mol_id": ["m1", "m2", "m3"],
            "H298_pm7": [-100.0, None, -50.0],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.null_h298_pm7 == 1
    assert any("null H298_pm7" in a for a in report.alerts)


# ── Test 5: target_delta_outliers ──


def test_target_delta_outliers(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """|delta| > 200 and > 500 → counted."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "OK", "OK"],
            "mol_id": ["m1", "m2", "m3", "m4"],
            "target_delta_kcalmol": [10.0, 250.0, -600.0, 5.0],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.target_delta_outliers_200 == 2  # 250, -600
    assert report.target_delta_outliers_500 == 1  # -600


# ── Test 6: homo_suspicious ──


def test_homo_suspicious(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """HOMO = -8696 → flagged (|HOMO| > 100)."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "OK"],
            "mol_id": ["m1", "m2", "m3"],
            "mopac_homo_ev": [-8.5, -8696.0, -9.1],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.homo_suspicious == 1
    assert any("|HOMO| > 100" in a for a in report.alerts)


# ── Test 7: gap_negative ──


def test_gap_negative(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """gap = -5.0 → critical alert."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK"],
            "mol_id": ["m1", "m2"],
            "mopac_gap_ev": [8.5, -5.0],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.gap_negative_count == 1
    assert any("negative HOMO-LUMO gap" in a for a in report.alerts)


# ── Test 8: missing_column ──


def test_missing_column(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """CSV without mopac_homo_ev → None, no crash."""
    df = pd.DataFrame(
        {
            "status": ["OK"],
            "mol_id": ["m1"],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.homo_mean is None
    assert report.homo_suspicious == 0
    assert report.h298_pm7_mean is None
    assert report.target_delta_mean is None


# ── Test 9: empty_csv ──


def test_empty_csv(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """Header-only CSV → zeros, no exception."""
    df = pd.DataFrame(columns=["status", "mol_id", "H298_pm7"])
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.total_rows == 0
    assert report.ok_count == 0
    assert report.h298_pm7_mean is None


# ── Test 10: conformers_zero ──


def test_conformers_zero(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """crest_conformers_generated = 0 → counted."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "OK"],
            "mol_id": ["m1", "m2", "m3"],
            "crest_conformers_generated": [0, 5, 0],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.conformers_zero_count == 2


# ── Test 11: crest_timeout_dynamic ──


def test_crest_timeout_dynamic(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """crest_time_s >= assigned_crest_timeout * 60 → counted as timeout."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "OK"],
            "mol_id": ["m1", "m2", "m3"],
            "crest_time_s": [1800.0, 500.0, 1801.0],  # m1 at limit, m3 over
            "assigned_crest_timeout": [30.0, 30.0, 30.0],  # 30 min = 1800s
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.crest_timeout_count == 2  # m1 (==1800) and m3 (>1800)


# ── Test 12: reruns_distribution ──


def test_reruns_distribution(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """Mix of reruns values → correct distribution (OK rows only), alert for reruns=3."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "OK", "SKIP", "OK"],
            "mol_id": ["m1", "m2", "m3", "m4", "m5"],
            "reruns": [0, 1, 0, 3, 3],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    # Pipeline health is scoped to OK rows only (m1, m2, m3, m5)
    assert report.reruns_distribution == {0: 2, 1: 1, 3: 1}
    assert report.max_reruns_count == 1  # m5 has reruns=3 and status=OK
    assert any("reruns=3" in a for a in report.alerts)


# ── Test 13: h298_extreme_per_nheavy ──


def test_h298_extreme_per_nheavy(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """H298_pm7/nheavy outside [-80, +20] → flagged as extreme."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "OK"],
            "mol_id": ["m1", "m2", "m3"],
            "H298_pm7": [-100.0, -500.0, 50.0],
            "nheavy": [5, 5, 2],
            # per_heavy: -20, -100, +25
            # -20 is within [-80, +20] → OK
            # -100 < -80 → extreme
            # +25 > +20 → extreme
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.h298_pm7_extreme_count == 2


# ── Test 14: crest_status SUCCESS not counted as failed ──


def test_crest_status_success_not_failed(
    analyzer: DatabaseAnalyzer, tmp_path: Path
) -> None:
    """crest_status='SUCCESS' (pipeline value) should not be counted as failed."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK", "OK"],
            "mol_id": ["m1", "m2", "m3"],
            "crest_status": ["SUCCESS", "ok", "success"],
            "mopac_status": ["OK", "OK", "OK"],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.crest_failed_count == 0


# ── Test 15: mopac_status OK uppercase not counted as failed ──


def test_mopac_status_ok_uppercase_not_failed(
    analyzer: DatabaseAnalyzer, tmp_path: Path
) -> None:
    """mopac_status='OK' (uppercase) should not be counted as failed."""
    df = pd.DataFrame(
        {
            "status": ["OK", "OK"],
            "mol_id": ["m1", "m2"],
            "mopac_status": ["OK", "ok"],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.mopac_failed_count == 0


# ── Test 16: cbs_suspect_flag ──


def test_cbs_suspect_flag(analyzer: DatabaseAnalyzer, tmp_path: Path) -> None:
    """Molecules with cbs_quality_flag=SUSPECT are counted separately."""
    df = pd.DataFrame(
        {
            "status": ["OK"] * 10,
            "mol_id": [f"m{i}" for i in range(10)],
            "cbs_quality_flag": ["SUSPECT", "SUSPECT"] + ["OK"] * 8,
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    assert report.cbs_suspect_count == 2
    assert report.cbs_ok_count == 8
    assert any("SUSPECT" in a for a in report.alerts)


def test_target_delta_excludes_cbs_suspect(
    analyzer: DatabaseAnalyzer, tmp_path: Path
) -> None:
    """SUSPECT CBS rows must not pollute Target Delta statistics."""
    df = pd.DataFrame(
        {
            "status": ["OK"] * 5,
            "mol_id": [f"m{i}" for i in range(5)],
            "target_delta_kcalmol": [10.0, 20.0, 15.0, -145000.0, -130000.0],
            "cbs_quality_flag": ["OK", "OK", "OK", "SUSPECT", "SUSPECT"],
        }
    )
    csv = _write_csv(tmp_path / "test.csv", df)
    report = analyzer.analyze(csv)

    # Mean should reflect only the 3 OK rows (10, 20, 15), not the −145k outliers
    assert report.target_delta_mean is not None
    assert abs(report.target_delta_mean - 15.0) < 1.0
    assert report.target_delta_outliers_500 == 0

    # Top delta outliers should not include SUSPECT molecules
    for entry in report.top_delta_outliers:
        assert abs(entry["delta"]) < 1000  # type: ignore[arg-type]
