"""Tests for batch prediction (Phase F).

8 tests covering predict_batch: columns, eligibility, math, I/O, errors, stats.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from grimperium.ml.predictor import predict_batch


class TestPredictBatch:
    """Tests for the predict_batch function."""

    def test_returns_dataframe_with_prediction_columns(
        self,
        csv_with_predictions: Path,
        trained_model_fixture: dict[str, Any],
        tmp_path: Path,
    ) -> None:
        """Result DataFrame has H298_predicted and delta_correction columns."""
        from grimperium.ml.persistence import save_model

        model_path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, model_path)

        df = predict_batch(csv_with_predictions, model_path)
        assert "H298_predicted" in df.columns
        assert "delta_correction" in df.columns

    def test_predictions_only_for_eligible_rows(
        self,
        csv_with_predictions: Path,
        trained_model_fixture: dict[str, Any],
        tmp_path: Path,
    ) -> None:
        """Only status==OK & cbs_quality_flag==OK rows get predictions; others NaN."""
        from grimperium.ml.persistence import save_model

        model_path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, model_path)

        df = predict_batch(csv_with_predictions, model_path)

        # 10 eligible rows (OK/OK), 2 non-eligible (SUSPECT, PENDING)
        eligible = df[(df["status"] == "OK") & (df["cbs_quality_flag"] == "OK")]
        non_eligible = df[~((df["status"] == "OK") & (df["cbs_quality_flag"] == "OK"))]

        assert eligible["H298_predicted"].notna().all()
        assert non_eligible["H298_predicted"].isna().all()

    def test_delta_correction_equals_predicted_minus_pm7(
        self,
        csv_with_predictions: Path,
        trained_model_fixture: dict[str, Any],
        tmp_path: Path,
    ) -> None:
        """delta_correction = H298_predicted - H298_pm7 for eligible rows."""
        from grimperium.ml.persistence import save_model

        model_path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, model_path)

        df = predict_batch(csv_with_predictions, model_path)

        eligible = df[(df["status"] == "OK") & (df["cbs_quality_flag"] == "OK")]
        expected_delta = eligible["H298_predicted"] - eligible["H298_pm7"]
        np.testing.assert_allclose(
            eligible["delta_correction"].to_numpy(),
            expected_delta.to_numpy(),
            rtol=1e-10,
        )

    def test_csv_is_written_with_new_columns(
        self,
        csv_with_predictions: Path,
        trained_model_fixture: dict[str, Any],
        tmp_path: Path,
    ) -> None:
        """CSV on disk is updated with new columns."""
        from grimperium.ml.persistence import save_model

        model_path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, model_path)

        predict_batch(csv_with_predictions, model_path)

        # Re-read from disk
        df_disk = pd.read_csv(csv_with_predictions)
        assert "H298_predicted" in df_disk.columns
        assert "delta_correction" in df_disk.columns

    def test_original_columns_preserved(
        self,
        csv_with_predictions: Path,
        trained_model_fixture: dict[str, Any],
        tmp_path: Path,
    ) -> None:
        """Existing columns are not altered."""
        from grimperium.ml.persistence import save_model

        model_path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, model_path)

        df_before = pd.read_csv(csv_with_predictions)
        original_cols = set(df_before.columns)

        predict_batch(csv_with_predictions, model_path)

        df_after = pd.read_csv(csv_with_predictions)
        assert original_cols.issubset(set(df_after.columns))

        # Verify original data unchanged for a known column
        pd.testing.assert_series_equal(df_before["H298_pm7"], df_after["H298_pm7"])

    def test_raises_on_missing_model(
        self,
        csv_with_predictions: Path,
        tmp_path: Path,
    ) -> None:
        """FileNotFoundError raised when model file does not exist."""
        fake_model = tmp_path / "nonexistent.joblib"
        with pytest.raises(FileNotFoundError, match="Model file not found"):
            predict_batch(csv_with_predictions, fake_model)

    def test_raises_on_missing_csv(
        self,
        trained_model_fixture: dict[str, Any],
        tmp_path: Path,
    ) -> None:
        """FileNotFoundError raised when CSV file does not exist."""
        from grimperium.ml.persistence import save_model

        model_path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, model_path)

        fake_csv = tmp_path / "nonexistent.csv"
        with pytest.raises(FileNotFoundError, match="CSV file not found"):
            predict_batch(fake_csv, model_path)

    def test_returns_summary_stats(
        self,
        csv_with_predictions: Path,
        trained_model_fixture: dict[str, Any],
        tmp_path: Path,
    ) -> None:
        """return_stats=True gives (df, stats_dict) with expected keys."""
        from grimperium.ml.persistence import save_model

        model_path = tmp_path / "model.joblib"
        save_model(trained_model_fixture, model_path)

        result = predict_batch(csv_with_predictions, model_path, return_stats=True)
        assert isinstance(result, tuple)
        assert len(result) == 2

        df, stats = result
        assert isinstance(df, pd.DataFrame)
        assert isinstance(stats, dict)
        assert "n_predicted" in stats
        assert "mean_delta" in stats
        assert "std_delta" in stats
        assert stats["n_predicted"] == 10  # 10 OK/OK rows
