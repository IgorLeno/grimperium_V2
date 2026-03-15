"""Tests for grimperium.ml.data_loader — ML-specific data loading."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from grimperium.ml.data_loader import load_ml_data


class TestLoadMLData:
    """Tests for load_ml_data function."""

    def test_returns_df_y_cbs_y_pm7(self, synthetic_csv_path: Path) -> None:
        """load_ml_data returns (DataFrame, y_cbs, y_pm7) tuple."""
        df, y_cbs, y_pm7 = load_ml_data(synthetic_csv_path)
        assert isinstance(df, pd.DataFrame)
        assert isinstance(y_cbs, np.ndarray)
        assert isinstance(y_pm7, np.ndarray)

    def test_filters_only_ok_status(self, synthetic_csv_path: Path) -> None:
        """Only rows with status='OK' are returned."""
        df, y_cbs, y_pm7 = load_ml_data(synthetic_csv_path)
        # All 10 rows in fixture have status=OK
        assert len(df) == 10
        assert len(y_cbs) == 10
        assert len(y_pm7) == 10

    def test_drops_rows_missing_targets(
        self, csv_with_missing_target: Path
    ) -> None:
        """Rows where H298_cbs or H298_pm7 is NaN are dropped."""
        df, y_cbs, y_pm7 = load_ml_data(csv_with_missing_target)
        assert len(df) == 1  # Only first row survives

    def test_raises_on_no_valid_rows(self, csv_no_ok_status: Path) -> None:
        """Raises ValueError if no rows remain after filtering."""
        with pytest.raises(ValueError, match="No valid"):
            load_ml_data(csv_no_ok_status)

    def test_preserves_feature_columns(self, synthetic_csv_path: Path) -> None:
        """Returned DataFrame contains feature columns needed by pipeline."""
        df, _, _ = load_ml_data(synthetic_csv_path)
        for col in [
            "nheavy",
            "rdkit_nrotbonds",
            "mopac_homo_ev",
            "mopac_lumo_ev",
            "mopac_gap_ev",
        ]:
            assert col in df.columns, f"Missing feature column: {col}"
