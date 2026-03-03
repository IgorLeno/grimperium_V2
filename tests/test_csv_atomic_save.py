"""Tests for atomic CSV save, backup, and recovery.

Covers:
- Atomic write prevents corruption on interrupted saves
- Backup (.bak) creation on every save
- Auto-recovery from 0-byte / corrupted CSV
- Orphan .tmp file cleanup on load
- CSVCorruptedError when both CSV and backup are broken
- Data fidelity across multiple save/load cycles
"""

import pandas as pd
import pytest

from grimperium.core.csv_utils import atomic_to_csv
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager, CSVCorruptedError


def _make_df(n: int = 5) -> pd.DataFrame:
    """Create a minimal DataFrame matching BatchCSVManager required columns."""
    return pd.DataFrame(
        {
            "mol_id": [f"mol_{i:03d}" for i in range(n)],
            "smiles": [f"C{'C' * i}" for i in range(n)],
            "nheavy": list(range(1, n + 1)),
            "status": ["PENDING"] * n,
        }
    )


class TestAtomicSaveSurvivesInterruption:
    """Original CSV must remain intact if write fails mid-stream."""

    def test_original_intact_after_write_error(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        csv_path = tmp_path / "thermo.csv"
        original_df = _make_df(3)
        original_df.to_csv(csv_path, index=False)
        original_content = csv_path.read_text()

        # Monkeypatch to_csv to explode during temp-file write
        real_to_csv = pd.DataFrame.to_csv

        def _exploding_to_csv(
            self: pd.DataFrame, *args: object, **kwargs: object
        ) -> None:
            raise OSError("Simulated disk failure")

        pd.DataFrame.to_csv = _exploding_to_csv  # type: ignore[assignment]
        try:
            with pytest.raises(OSError, match="Simulated disk failure"):
                atomic_to_csv(csv_path, _make_df(10))
        finally:
            pd.DataFrame.to_csv = real_to_csv  # type: ignore[assignment]

        # Original must be untouched
        assert csv_path.read_text() == original_content
        # No orphan temp files
        assert list(tmp_path.glob("*.csv.tmp")) == []


class TestAtomicSaveCreatesBackup:
    """Every atomic save of an existing file must create a .bak."""

    def test_backup_created_with_all_columns(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        csv_path = tmp_path / "thermo.csv"
        df = _make_df(5)
        # First save (no existing file -> no .bak yet)
        atomic_to_csv(csv_path, df)
        assert not (tmp_path / "thermo.csv.bak").exists()

        # Second save -> .bak should appear
        df2 = _make_df(8)
        atomic_to_csv(csv_path, df2)
        bak = tmp_path / "thermo.csv.bak"
        assert bak.exists()

        # Backup should contain the FIRST save's data (5 rows)
        bak_df = pd.read_csv(bak)
        assert len(bak_df) == 5

        # Current file should have the new data (8 rows)
        current_df = pd.read_csv(csv_path)
        assert len(current_df) == 8


class TestLoadCSVRecoversFromBackup:
    """0-byte CSV + valid .bak -> auto-recovery."""

    def test_zero_byte_csv_recovered(self, tmp_path: "pytest.TempPathFactory") -> None:
        csv_path = tmp_path / "thermo_pm7.csv"
        bak_path = tmp_path / "thermo_pm7.csv.bak"

        # Create valid backup
        df = _make_df(4)
        df.to_csv(bak_path, index=False)

        # Create 0-byte primary CSV
        csv_path.write_text("")

        mgr = BatchCSVManager(csv_path=csv_path)
        result = mgr.load_csv()

        assert len(result) == 4
        assert list(result["mol_id"]) == ["mol_000", "mol_001", "mol_002", "mol_003"]
        # Primary file should now be restored
        assert csv_path.stat().st_size > 0


class TestOrphanTmpFilesCleaned:
    """Stale .tmp files from interrupted writes are removed on load."""

    def test_tmp_files_removed_on_load(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        csv_path = tmp_path / "thermo_pm7.csv"
        df = _make_df(3)
        df.to_csv(csv_path, index=False)

        # Create orphan temp files
        (tmp_path / ".atomic_abc123.csv.tmp").write_text("garbage")
        (tmp_path / ".atomic_def456.csv.tmp").write_text("more garbage")

        mgr = BatchCSVManager(csv_path=csv_path)
        mgr.load_csv()

        assert list(tmp_path.glob("*.csv.tmp")) == []


class TestBothCSVAndBackupCorrupted:
    """Raises CSVCorruptedError with actionable message."""

    def test_raises_csv_corrupted_error(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        csv_path = tmp_path / "thermo_pm7.csv"
        bak_path = tmp_path / "thermo_pm7.csv.bak"

        # Both empty
        csv_path.write_text("")
        bak_path.write_text("")

        mgr = BatchCSVManager(csv_path=csv_path)
        with pytest.raises(CSVCorruptedError, match="corrupted/missing"):
            mgr.load_csv()

    def test_raises_when_no_backup_exists(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        csv_path = tmp_path / "thermo_pm7.csv"
        csv_path.write_text("")  # 0-byte, no .bak

        mgr = BatchCSVManager(csv_path=csv_path)
        with pytest.raises(CSVCorruptedError, match="corrupted/missing"):
            mgr.load_csv()


class TestDataFidelityAcrossCycles:
    """Float precision and all columns persist across 10 save/load cycles."""

    def test_ten_cycles_preserve_precision(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        csv_path = tmp_path / "thermo_pm7.csv"

        df = _make_df(5)
        df["h298_pm7"] = [
            -123.456789012345,
            0.0,
            1e-10,
            9999999.99999,
            -0.000001,
        ]
        df["dipole_moment"] = [1.23456, 2.34567, 3.45678, 4.56789, 5.67890]

        mgr = BatchCSVManager(csv_path=csv_path)
        mgr.df = df

        for _ in range(10):
            mgr.save_csv()
            mgr.load_csv()

        assert len(mgr.df) == 5
        assert list(mgr.df.columns) == list(df.columns)
        # Float precision: pandas default to_csv uses repr precision
        pd.testing.assert_frame_equal(
            mgr.df[["h298_pm7", "dipole_moment"]].astype(float),
            df[["h298_pm7", "dipole_moment"]].astype(float),
        )


class TestAtomicToCsvEdgeCases:
    """Edge cases for the atomic_to_csv utility itself."""

    def test_none_path_raises_type_error(self) -> None:
        with pytest.raises(TypeError, match="path must not be None"):
            atomic_to_csv(None, _make_df())  # type: ignore[arg-type]

    def test_none_df_raises_type_error(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        with pytest.raises(TypeError, match="df must not be None"):
            atomic_to_csv(tmp_path / "x.csv", None)  # type: ignore[arg-type]

    def test_creates_parent_directories(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        deep = tmp_path / "a" / "b" / "c" / "data.csv"
        atomic_to_csv(deep, _make_df(2))
        assert deep.exists()
        assert len(pd.read_csv(deep)) == 2

    def test_missing_primary_csv_recovered_from_backup(
        self, tmp_path: "pytest.TempPathFactory"
    ) -> None:
        csv_path = tmp_path / "thermo_pm7.csv"
        bak_path = tmp_path / "thermo_pm7.csv.bak"

        # Only backup exists, primary is missing entirely
        df = _make_df(3)
        df.to_csv(bak_path, index=False)

        mgr = BatchCSVManager(csv_path=csv_path)
        result = mgr.load_csv()
        assert len(result) == 3
