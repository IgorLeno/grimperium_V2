"""
Results view for GRIMPERIUM CLI.

Displays performance analytics and divergence analysis using real model data.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from rich.panel import Panel
from rich.table import Table

from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView
from grimperium.core.metrics import mae, max_error, r2_score, rmse
from grimperium.ml.gate import evaluate_gate
from grimperium.ml.persistence import load_model_metadata, save_model

if TYPE_CHECKING:
    from grimperium import DictStrAny

logger = logging.getLogger(__name__)

# Static thresholds for divergence severity (no env dependency)
_DIVERGENCE_THRESHOLDS = [
    ("LOW", 0, 10),
    ("MEDIUM", 10, 25),
    ("HIGH", 25, 50),
    ("CRITICAL", 50, float("inf")),
]


class ResultsView(BaseView):
    """View for performance analytics and results."""

    name = "results"
    title = "Results"
    icon = ICONS["results"]
    color = COLORS["results"]

    def _get_model_path(self) -> Path:
        """Resolve model path from env var (read fresh on each call).

        Returns:
            Path: The resolved model file path.
        """
        return Path(
            os.environ.get("GRIMPERIUM_MODEL_PATH", "models/delta_learner_v1.joblib")
        )

    def _get_csv_path(self) -> Path:
        """Resolve CSV path from env var (read fresh on each call).

        Returns:
            Path: The resolved CSV data file path.
        """
        return Path(os.environ.get("GRIMPERIUM_DATA_PATH", "data/thermo_pm7.csv"))

    def _load_real_model_row(self) -> dict[str, Any] | None:
        """Load model metadata for display. Returns None if not available."""
        try:
            return load_model_metadata(self._get_model_path())
        except (FileNotFoundError, KeyError, Exception):
            return None

    def _compute_divergence_stats(self) -> list[DictStrAny] | None:
        """Compute divergence statistics from predictions vs CBS values.

        Returns None if CSV is missing or has no predictions.
        """
        csv_path = self._get_csv_path()
        if not csv_path.exists():
            return None

        df = pd.read_csv(csv_path)
        if "H298_predicted" not in df.columns or "H298_cbs" not in df.columns:
            return None

        # Only rows with both predicted and CBS values
        valid = df.dropna(subset=["H298_predicted", "H298_cbs"])
        if len(valid) == 0:
            return None

        # Compute relative divergence: |predicted - cbs| / |cbs| * 100
        cbs_abs = valid["H298_cbs"].abs()
        # Avoid division by zero: use small epsilon
        safe_cbs = cbs_abs.replace(0, np.finfo(float).eps)
        divergence_pct = (
            (valid["H298_predicted"] - valid["H298_cbs"]).abs() / safe_cbs * 100
        )

        total = len(divergence_pct)
        stats: list[DictStrAny] = []
        for severity, lo, hi in _DIVERGENCE_THRESHOLDS:
            count = int(((divergence_pct >= lo) & (divergence_pct < hi)).sum())
            pct = (count / total * 100) if total > 0 else 0.0
            stats.append(
                {
                    "severity": severity,
                    "range_min": lo,
                    "range_max": hi if hi != float("inf") else 100,
                    "count": count,
                    "percentage": pct,
                }
            )
        return stats

    def render(self) -> None:
        """Render the results overview."""
        self.clear_screen()
        self.show_header()

        self._render_model_comparison()
        self._render_divergence_analysis()

    def _render_model_comparison(self) -> None:
        """Render model performance from real metadata."""
        table = Table(
            title="Model Performance",
            show_header=True,
            header_style=f"bold {COLORS['results']}",
            border_style=COLORS["border"],
        )
        table.add_column("Model", style="bold")
        table.add_column("Algorithm")
        table.add_column("MAE (kcal/mol)", justify="right")
        table.add_column("R\u00b2", justify="right")
        table.add_column("Status")

        metadata = self._load_real_model_row()
        if metadata is not None:
            test_metrics = metadata.get("metrics", {}).get("test", {})
            mae_val = test_metrics.get("mae")
            r2_val = test_metrics.get("r2")
            mae_str = f"{mae_val:.3f}" if mae_val is not None else "-"
            r2_str = f"{r2_val:.4f}" if r2_val is not None else "-"
            status = (
                f"[{COLORS['success']}]{ICONS['success']} Ready"
                f"[/{COLORS['success']}]"
            )
            table.add_row(
                "DeltaLearner v1",
                "KRR + XGBoost Ensemble",
                mae_str,
                r2_str,
                status,
            )
        else:
            table.add_row(
                "DeltaLearner v1",
                "KRR + XGBoost Ensemble",
                "-",
                "-",
                f"[{COLORS['in_dev']}]{ICONS['in_dev']} Not Trained"
                f"[/{COLORS['in_dev']}]",
            )

        self.console.print(table)
        self.console.print()

    def _render_divergence_analysis(self) -> None:
        """Render CBS vs Predicted divergence analysis from real data."""
        stats = self._compute_divergence_stats()

        if stats is None:
            self.console.print(
                Panel(
                    f"[{COLORS['muted']}]No predictions available. "
                    f"Run 'Predict Batch' from this menu to generate predictions."
                    f"[/{COLORS['muted']}]",
                    title=f"[bold {COLORS['results']}]Divergence Analysis"
                    f"[/bold {COLORS['results']}]",
                    border_style=COLORS["border"],
                    padding=(1, 2),
                )
            )
            self.console.print()
            return

        # Divergence distribution table
        table = Table(
            title="Predicted vs CBS Divergence Distribution",
            show_header=True,
            header_style=f"bold {COLORS['results']}",
            border_style=COLORS["border"],
        )
        table.add_column("Severity", style="bold")
        table.add_column("Range (%)", justify="center")
        table.add_column("Count", justify="right")
        table.add_column("Percentage", justify="right")
        table.add_column("Bar", min_width=20)

        severity_colors = {
            "LOW": COLORS["success"],
            "MEDIUM": COLORS["warning"],
            "HIGH": COLORS["high"],
            "CRITICAL": COLORS["error"],
        }

        total_molecules = sum(s["count"] for s in stats)

        for stat in stats:
            color = severity_colors.get(stat["severity"], COLORS["muted"])
            pct = max(0, min(stat["percentage"], 100))
            bar_length = max(0, min(int(pct / 5), 20))
            bar = (
                f"[{color}]"
                f"{'█' * bar_length}{'░' * (20 - bar_length)}"
                f"[/{color}]"
            )

            range_max = stat["range_max"]
            range_str = (
                f"{stat['range_min']:.0f}% - {range_max:.0f}%"
                if range_max < 100
                else f"{stat['range_min']:.0f}%+"
            )

            table.add_row(
                f"[{color}]{stat['severity']}[/{color}]",
                range_str,
                f"{stat['count']:,}",
                f"{stat['percentage']:.1f}%",
                bar,
            )

        self.console.print(table)
        self.console.print()

        # Summary panel
        low_medium = sum(
            s["count"] for s in stats if s["severity"] in ["LOW", "MEDIUM"]
        )
        low_medium_pct = (
            (low_medium / total_molecules) * 100 if total_molecules > 0 else 0
        )

        summary = f"""
[bold]Key Findings:[/bold]

\u2022 Total molecules analyzed: {total_molecules:,}
\u2022 Molecules with LOW/MEDIUM divergence: {low_medium:,} ({low_medium_pct:.1f}%)
\u2022 The Delta-Learning approach effectively corrects PM7 predictions

[bold]Interpretation:[/bold]

\u2022 [{COLORS['success']}]LOW (0-10%)[/{COLORS['success']}]: Predictions are accurate, small corrections needed
\u2022 [{COLORS['warning']}]MEDIUM (10-25%)[/{COLORS['warning']}]: Moderate corrections, ML performs well
\u2022 [{COLORS['high']}]HIGH (25-50%)[/{COLORS['high']}]: Significant corrections needed, challenging cases
\u2022 [{COLORS['error']}]CRITICAL (>50%)[/{COLORS['error']}]: Large deviations, may require special handling
"""
        self.console.print(
            Panel(
                summary,
                title=f"[bold {COLORS['results']}]Analysis Summary"
                f"[/bold {COLORS['results']}]",
                border_style=COLORS["border"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the results view."""
        return [
            MenuOption(
                label="Run New Analysis",
                value="run_analysis",
                icon=ICONS.get("models", "\U0001f916"),
            ),
            MenuOption(
                label="Predict Batch",
                value="predict_batch",
                icon=ICONS.get("calc", "\u26a1"),
            ),
            MenuOption(
                label="Detailed Metrics",
                value="detailed",
                icon=ICONS.get("results", "\U0001f4ca"),
            ),
            MenuOption(
                label="Visualization Charts",
                value="charts",
                icon=ICONS.get("results", "📊"),
            ),
        ]

    def handle_action(self, action: str) -> str | None:
        """Handle menu actions."""
        if action == "back":
            return "main"

        if action == "run_analysis":
            self._handle_run_analysis()
            return None

        if action == "predict_batch":
            self._handle_predict_batch()
            return None

        if action == "charts":
            self._handle_charts()
            return None

        if action == "detailed":
            self._handle_detailed_metrics()
            return None

        return None

    def _get_charts_dir(self) -> Path:
        """Resolve charts output directory from env var.

        Returns:
            Path: The resolved charts output directory.
        """
        return Path(os.environ.get("GRIMPERIUM_CHARTS_DIR", "reports/charts"))

    def _handle_detailed_metrics(self) -> None:
        """Compute and display detailed statistical metrics for predictions."""
        csv_path = self._get_csv_path()
        if not csv_path.exists():
            self.show_error(f"Data file not found: {csv_path}")
            self.wait_for_enter()
            return

        df = pd.read_csv(csv_path)
        if "H298_predicted" not in df.columns or "H298_cbs" not in df.columns:
            self.show_error(
                "CSV must contain both 'H298_predicted' and 'H298_cbs' columns."
            )
            self.wait_for_enter()
            return

        valid = df.dropna(subset=["H298_predicted", "H298_cbs"])
        if len(valid) == 0:
            self.show_error("No valid rows with both predicted and CBS values.")
            self.wait_for_enter()
            return

        y_cbs = valid["H298_cbs"].to_numpy()
        y_pred = valid["H298_predicted"].to_numpy()

        # Core metrics from metrics.py
        mae_val = mae(y_cbs, y_pred)
        rmse_val = rmse(y_cbs, y_pred)
        r2_val = r2_score(y_cbs, y_pred)
        max_err = max_error(y_cbs, y_pred)

        # Additional statistics
        errors = y_pred - y_cbs
        abs_errors = np.abs(errors)
        bias = float(np.mean(errors))

        # Pearson r (handle constant arrays gracefully)
        if np.std(y_cbs) == 0 or np.std(y_pred) == 0:
            pearson_r = float("nan")
        else:
            pearson_r = float(np.corrcoef(y_cbs, y_pred)[0, 1])

        # Percentage within thresholds
        n = len(abs_errors)
        pct_1 = float(np.sum(abs_errors <= 1.0) / n * 100)
        pct_2 = float(np.sum(abs_errors <= 2.0) / n * 100)
        pct_5 = float(np.sum(abs_errors <= 5.0) / n * 100)

        # Display table
        table = Table(
            title="Detailed Prediction Metrics",
            show_header=True,
            header_style=f"bold {COLORS['results']}",
            border_style=COLORS["border"],
            show_edge=False,
        )
        table.add_column("Metric", style="bold")
        table.add_column("Value", justify="right")

        table.add_row("MAE (kcal/mol)", f"{mae_val:.4f}")
        table.add_row("RMSE (kcal/mol)", f"{rmse_val:.4f}")
        table.add_row("Max Error (kcal/mol)", f"{max_err:.4f}")
        table.add_row("Bias (kcal/mol)", f"{bias:.4f}")
        table.add_row("R\u00b2", f"{r2_val:.4f}")
        table.add_row(
            "Pearson r",
            f"{pearson_r:.5f}" if not np.isnan(pearson_r) else "N/A",
        )
        table.add_row("", "")
        table.add_row("Within \u00b11 kcal/mol", f"{pct_1:.1f}%")
        table.add_row("Within \u00b12 kcal/mol", f"{pct_2:.1f}%")
        table.add_row("Within \u00b15 kcal/mol", f"{pct_5:.1f}%")
        table.add_row("", "")
        table.add_row("Molecules analyzed", f"{n:,}")

        self.console.print()
        self.console.print(table)
        self.console.print()

        # Interpretation panel
        interpretation = (
            f"[bold]Interpretation:[/bold]\n\n"
            f"\u2022 Mean absolute error: {mae_val:.4f} kcal/mol\n"
            f"\u2022 Systematic bias: {bias:+.4f} kcal/mol "
            f"({'overestimates' if bias > 0 else 'underestimates' if bias < 0 else 'unbiased'})\n"
            f"\u2022 {pct_1:.1f}% of predictions within chemical accuracy (\u00b11 kcal/mol)\n"
            f"\u2022 R\u00b2 = {r2_val:.4f} — "
            f"{'excellent' if r2_val >= 0.999 else 'good' if r2_val >= 0.99 else 'moderate' if r2_val >= 0.95 else 'needs improvement'} fit"
        )
        self.console.print(
            Panel(
                interpretation,
                title=f"[bold {COLORS['results']}]Analysis Summary"
                f"[/bold {COLORS['results']}]",
                border_style=COLORS["border"],
                padding=(1, 2),
            )
        )
        self.console.print()
        self.wait_for_enter()

    def _handle_charts(self) -> None:
        """Generate visualization charts and display results."""
        from grimperium.ml.charts import ChartGenerationResult, generate_charts

        csv_path = self._get_csv_path()
        charts_dir = self._get_charts_dir()

        self.console.print()
        self.console.print(
            f"[bold {COLORS['results']}]Generating charts..."
            f"[/bold {COLORS['results']}]"
        )
        self.console.print()

        try:
            result: ChartGenerationResult = generate_charts(csv_path, charts_dir)

            result_text = f"""
[bold]Charts Generated Successfully![/bold]

[bold]Summary:[/bold]
  Data points plotted:  {result.n_points}
  RMSE:                 {result.rmse:.4f} kcal/mol
  R\u00b2:                   {result.r2:.4f}

[bold]Files:[/bold]
  Parity plot:      {result.parity_plot}
  Delta histogram:  {result.delta_histogram}
  Residuals plot:   {result.residuals_plot}
"""

            self.console.print(
                Panel(
                    result_text,
                    title=f"[bold {COLORS['success']}]Visualization Charts"
                    f"[/bold {COLORS['success']}]",
                    border_style=COLORS["success"],
                    padding=(1, 2),
                )
            )
            self.show_success("Charts saved!")
        except FileNotFoundError as e:
            self.show_error(str(e))
        except ValueError as e:
            self.show_error(str(e))
        except Exception as e:
            self.show_error(f"Chart generation failed: {e}")

        self.wait_for_enter()

    def _handle_predict_batch(self) -> None:
        """Run batch prediction and display results."""
        from grimperium.ml.predictor import predict_batch

        model_path = self._get_model_path()
        csv_path = self._get_csv_path()

        self.console.print()
        self.console.print(
            f"[bold {COLORS['results']}]Running batch prediction..."
            f"[/bold {COLORS['results']}]"
        )
        self.console.print()

        try:
            _df, stats = predict_batch(csv_path, model_path, return_stats=True)

            result_text = f"""
[bold]Batch Prediction Complete![/bold]

[bold]Summary:[/bold]
  Molecules predicted:  {stats['n_predicted']}
  Mean delta correction: {stats['mean_delta']:.4f} kcal/mol
  Std delta correction:  {stats['std_delta']:.4f} kcal/mol
  Min predicted H298:    {stats['min_predicted']:.4f} kcal/mol
  Max predicted H298:    {stats['max_predicted']:.4f} kcal/mol

Results written to: {csv_path}
"""

            self.console.print(
                Panel(
                    result_text,
                    title=f"[bold {COLORS['success']}]Prediction Results"
                    f"[/bold {COLORS['success']}]",
                    border_style=COLORS["success"],
                    padding=(1, 2),
                )
            )
            self.show_success("Batch prediction complete!")
        except FileNotFoundError as e:
            self.show_error(str(e))
        except Exception as e:
            self.show_error(f"Prediction failed: {e}")

        self.wait_for_enter()

    def _handle_run_analysis(self) -> None:
        """Run complete analysis: retrain model, predict batch, regenerate charts."""
        from grimperium.ml.charts import generate_charts
        from grimperium.ml.predictor import predict_batch
        from grimperium.ml.trainer import train

        model_path = self._get_model_path()
        csv_path = self._get_csv_path()
        charts_dir = self._get_charts_dir()

        # Capture "before" metrics if model exists
        old_metrics: dict[str, Any] | None = None
        if model_path.exists():
            try:
                old_meta = load_model_metadata(model_path)
                old_metrics = old_meta.get("metrics", {}).get("test", {})
            except Exception as exc:
                logger.debug("Could not load previous model metrics: %s", exc)
                old_metrics = None

        # Count eligible molecules
        if not csv_path.exists():
            self.show_error(f"Data file not found: {csv_path}")
            self.wait_for_enter()
            return

        try:
            df_check = pd.read_csv(csv_path)
            n_eligible = int(
                (
                    (df_check["status"] == "OK")
                    & (df_check["cbs_quality_flag"] == "OK")
                ).sum()
            )
        except Exception as exc:
            logger.debug("Could not count eligible molecules: %s", exc)
            n_eligible = 0

        # Confirmation prompt
        self.console.print()
        self.console.print(
            f"[bold {COLORS['results']}]Run New Analysis"
            f"[/bold {COLORS['results']}]\n"
            f"  Data: {csv_path} ({n_eligible:,} eligible molecules)\n"
            f"  Model: {model_path}\n\n"
            f"  This will:\n"
            f"    1. Retrain the model with current data\n"
            f"    2. Run batch predictions for all eligible molecules\n"
            f"    3. Regenerate visualization charts\n"
        )
        confirm = self.console.input("Type 'yes' to confirm: ").strip().lower()
        if confirm not in ("yes", "y"):
            self.console.print(
                f"[{COLORS['muted']}]Analysis cancelled.[/{COLORS['muted']}]"
            )
            return

        # --- Step 1: Retrain ---
        self.console.print()
        self.console.print(
            f"[bold {COLORS['results']}]Step 1/3: Retraining model..."
            f"[/bold {COLORS['results']}] (this may take a minute)"
        )

        try:
            result = train(csv_path, return_pipeline=True, random_state=42)
            if not isinstance(result, tuple) or len(result) != 4:
                msg = (
                    f"train(return_pipeline=True) returned "
                    f"{type(result).__name__} with "
                    f"{len(result) if isinstance(result, tuple) else 'N/A'} "
                    f"elements; expected a 4-tuple"
                )
                raise TypeError(msg)
            learner, train_m, test_m, pipeline = result
            test_m["gate_pass"] = evaluate_gate(test_m)

            bundle: dict[str, Any] = {
                "learner": learner,
                "pipeline": pipeline,
                "metrics": {"train": train_m, "test": test_m},
            }
            save_model(bundle, model_path)
            self.console.print(
                f"  [{COLORS['success']}]{ICONS['success']} Model retrained"
                f" and saved[/{COLORS['success']}]"
            )
        except Exception as e:
            self.show_error(f"Training failed: {e}")
            self.wait_for_enter()
            return

        # --- Step 2: Predict Batch ---
        self.console.print(
            f"[bold {COLORS['results']}]Step 2/3: Running batch predictions..."
            f"[/bold {COLORS['results']}]"
        )

        try:
            _df, pred_stats = predict_batch(csv_path, model_path, return_stats=True)
            self.console.print(
                f"  [{COLORS['success']}]{ICONS['success']} "
                f"{pred_stats['n_predicted']:,} molecules predicted"
                f"[/{COLORS['success']}]"
            )
        except Exception as e:
            self.show_error(f"Batch prediction failed: {e}")
            self.wait_for_enter()
            return

        # --- Step 3: Generate Charts ---
        self.console.print(
            f"[bold {COLORS['results']}]Step 3/3: Generating charts..."
            f"[/bold {COLORS['results']}]"
        )

        chart_result = None
        try:
            chart_result = generate_charts(csv_path, charts_dir)
            self.console.print(
                f"  [{COLORS['success']}]{ICONS['success']} Charts saved"
                f" to {charts_dir}[/{COLORS['success']}]"
            )
        except Exception as e:
            self.console.print(
                f"  [{COLORS['warning']}]{ICONS['warning']} "
                f"Chart generation failed: {e}[/{COLORS['warning']}]"
            )

        # --- Display Results ---
        self.console.print()

        # Before/after comparison table
        if old_metrics is not None:
            comp = Table(
                title="Before vs After Analysis",
                show_header=True,
                header_style=f"bold {COLORS['results']}",
                border_style=COLORS["border"],
            )
            comp.add_column("Metric")
            comp.add_column("Before", justify="right")
            comp.add_column("After", justify="right")
            comp.add_column("Delta", justify="right")
            comp.add_column("")

            def _delta_row(
                label: str,
                before: Any,
                after: Any,
                *,
                lower_is_better: bool = True,
            ) -> None:
                b = float(before) if before is not None else 0.0
                a = float(after) if after is not None else 0.0
                delta = a - b
                improved = (delta < 0) if lower_is_better else (delta > 0)
                icon = ICONS["success"] if improved else ICONS["error"]
                comp.add_row(label, f"{b:.4f}", f"{a:.4f}", f"{delta:+.4f}", icon)

            _delta_row("MAE (kcal/mol)", old_metrics.get("mae"), test_m.get("mae"))
            _delta_row("RMSE (kcal/mol)", old_metrics.get("rmse"), test_m.get("rmse"))
            _delta_row(
                "R\u00b2",
                old_metrics.get("r2"),
                test_m.get("r2"),
                lower_is_better=False,
            )
            _delta_row(
                "Max Error (kcal/mol)",
                old_metrics.get("max_error"),
                test_m.get("max_error"),
            )

            self.console.print(comp)
            self.console.print()

        # Summary panel
        gate_pass = test_m.get("gate_pass", False)
        gate_icon = ICONS["success"] if gate_pass else ICONS["error"]

        mae_val = float(test_m.get("mae", 0))
        r2_val = float(test_m.get("r2", 0))

        summary = (
            f"[bold]Analysis Complete![/bold]\n\n"
            f"[bold]Model:[/bold]\n"
            f"  MAE:         {mae_val:.3f} kcal/mol\n"
            f"  R\u00b2:          {r2_val:.4f}\n"
            f"  Gate Pass:   {gate_icon}"
            f" {'Yes' if gate_pass else 'No'}\n\n"
            f"[bold]Predictions:[/bold]\n"
            f"  Molecules:   {pred_stats['n_predicted']:,}\n"
            f"  Mean delta:  {pred_stats['mean_delta']:.4f} kcal/mol"
        )

        if chart_result is not None:
            summary += (
                f"\n\n[bold]Charts:[/bold]\n"
                f"  Parity plot:      {chart_result.parity_plot}\n"
                f"  Delta histogram:  {chart_result.delta_histogram}\n"
                f"  Residuals plot:   {chart_result.residuals_plot}"
            )

        self.console.print(
            Panel(
                summary,
                title=f"[bold {COLORS['success']}]New Analysis Results"
                f"[/bold {COLORS['success']}]",
                border_style=COLORS["success"],
                padding=(1, 2),
            )
        )
        self.show_success("Full analysis pipeline complete!")
        self.wait_for_enter()

    def run(self) -> str | None:
        """Run the results view interaction loop."""
        while True:
            self.render()
            result = show_back_menu(
                options=self.get_menu_options(),
                title="Actions",
            )

            if result is None or result == "back":
                return "main"
            else:
                next_view = self.handle_action(result)
                if next_view:
                    return next_view
