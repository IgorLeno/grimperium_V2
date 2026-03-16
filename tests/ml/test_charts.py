"""Tests for chart generation (Phase G)."""

from __future__ import annotations

from pathlib import Path

import pytest

from grimperium.ml.charts import ChartGenerationResult, generate_charts


class TestGenerateCharts:
    """Tests for generate_charts function."""

    def test_returns_chart_result_object(
        self, tmp_path: Path, csv_with_predictions_for_charts: Path
    ) -> None:
        """generate_charts returns a ChartGenerationResult dataclass."""
        output_dir = tmp_path / "charts_out"
        result = generate_charts(csv_with_predictions_for_charts, output_dir)
        assert isinstance(result, ChartGenerationResult)

    def test_parity_plot_created(
        self, tmp_path: Path, csv_with_predictions_for_charts: Path
    ) -> None:
        """Parity plot PNG is created on disk."""
        output_dir = tmp_path / "charts_out"
        result = generate_charts(csv_with_predictions_for_charts, output_dir)
        assert result.parity_plot.exists()
        assert result.parity_plot.suffix == ".png"

    def test_delta_histogram_created(
        self, tmp_path: Path, csv_with_predictions_for_charts: Path
    ) -> None:
        """Delta histogram PNG is created on disk."""
        output_dir = tmp_path / "charts_out"
        result = generate_charts(csv_with_predictions_for_charts, output_dir)
        assert result.delta_histogram.exists()
        assert result.delta_histogram.suffix == ".png"

    def test_residuals_plot_created(
        self, tmp_path: Path, csv_with_predictions_for_charts: Path
    ) -> None:
        """Residuals plot PNG is created on disk."""
        output_dir = tmp_path / "charts_out"
        result = generate_charts(csv_with_predictions_for_charts, output_dir)
        assert result.residuals_plot.exists()
        assert result.residuals_plot.suffix == ".png"

    def test_output_dir_created_if_missing(
        self, tmp_path: Path, csv_with_predictions_for_charts: Path
    ) -> None:
        """Output directory is created automatically if it doesn't exist."""
        output_dir = tmp_path / "nested" / "charts"
        assert not output_dir.exists()
        generate_charts(csv_with_predictions_for_charts, output_dir)
        assert output_dir.exists()

    def test_raises_on_missing_csv(self, tmp_path: Path) -> None:
        """FileNotFoundError is raised when CSV path doesn't exist."""
        fake_csv = tmp_path / "nonexistent.csv"
        output_dir = tmp_path / "charts_out"
        with pytest.raises(FileNotFoundError):
            generate_charts(fake_csv, output_dir)

    def test_raises_if_no_predictions(
        self, tmp_path: Path, synthetic_csv_path: Path
    ) -> None:
        """ValueError is raised when CSV has no H298_predicted column."""
        output_dir = tmp_path / "charts_out"
        with pytest.raises(ValueError, match="H298_predicted"):
            generate_charts(synthetic_csv_path, output_dir)

    def test_chart_result_has_n_points(
        self, tmp_path: Path, csv_with_predictions_for_charts: Path
    ) -> None:
        """ChartGenerationResult.n_points is a positive integer."""
        output_dir = tmp_path / "charts_out"
        result = generate_charts(csv_with_predictions_for_charts, output_dir)
        assert result.n_points > 0

    def test_chart_result_has_metrics(
        self, tmp_path: Path, csv_with_predictions_for_charts: Path
    ) -> None:
        """ChartGenerationResult has rmse and r2 as floats."""
        output_dir = tmp_path / "charts_out"
        result = generate_charts(csv_with_predictions_for_charts, output_dir)
        assert isinstance(result.rmse, float)
        assert isinstance(result.r2, float)
