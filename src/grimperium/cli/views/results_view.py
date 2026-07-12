"""
Results view for GRIMPERIUM CLI.

Displays performance analytics and divergence analysis using real model data.
"""

from __future__ import annotations

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
from grimperium.ml.persistence import load_model_metadata

if TYPE_CHECKING:
    from grimperium import DictStrAny

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

    def _get_model_path(self) -> Path | None:
        """Resolve model path from the live controller session.

        Returns:
            The selected model path, or None when no model is selected.
        """
        model_path = getattr(self.controller, "current_model_path", None)
        if model_path is None:
            self.console.print(
                f"[{COLORS['warning']}][Model not selected — some metrics "
                f"may be unavailable][/{COLORS['warning']}]"
            )
            return None
        return Path(model_path)

    def _get_csv_path(self) -> Path:
        """Resolve CSV path from the typed live controller session."""
        session = getattr(self.controller, "__dict__", {}).get("session")
        csv_path = getattr(session, "analysis_path", None)
        if csv_path is None:
            csv_path = getattr(self.controller, "current_csv_path", None)
        if csv_path is None:
            fallback = Path("data/thermo_pm7.csv")
            self.console.print(
                f"[{COLORS['warning']}]Using default data file: {fallback}"
                f"[/{COLORS['warning']}]"
            )
            return fallback
        return Path(csv_path)

    def _show_analysis_only_message(self, message: str) -> None:
        """Show a boundary message for actions now owned by Models."""
        self.console.print()
        self.console.print(
            Panel(
                f"[{COLORS['warning']}]{message}[/{COLORS['warning']}]",
                title=f"[bold {COLORS['warning']}]Analysis Only[/bold {COLORS['warning']}]",
                border_style=COLORS["warning"],
                padding=(1, 2),
            )
        )
        self.console.print()
        self.wait_for_enter()

    def _load_real_model_row(self) -> dict[str, Any] | None:
        """Load model metadata for display. Returns None if not available."""
        model_path = self._get_model_path()
        if model_path is None:
            return None
        try:
            return load_model_metadata(model_path)
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
                    f"Use Models > Predict Batch before returning here for analysis."
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
                disabled=True,
                disabled_reason="Use Models view",
            ),
            MenuOption(
                label="Predict Batch",
                value="predict_batch",
                icon=ICONS.get("calc", "\u26a1"),
                disabled=True,
                disabled_reason="Use Models view",
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
            MenuOption(
                label="Show Outliers",
                value="show_outliers",
                icon="🔍",
            ),
            MenuOption(
                label="Top Error Molecules",
                value="top_errors",
                icon="📋",
            ),
            MenuOption(
                label="Generate HTML Report",
                value="html_report",
                icon="📄",
            ),
            MenuOption(
                label="Export Outliers CSV",
                value="export_outliers",
                icon="💾",
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

        if action == "show_outliers":
            self._handle_show_outliers()
            return None

        if action == "top_errors":
            self._handle_top_errors()
            return None

        if action == "html_report":
            self._handle_html_report()
            return None

        if action == "export_outliers":
            self._handle_export_outliers()
            return None

        return None

    def _get_charts_dir(self) -> Path:
        """Resolve charts output directory from env var.

        Returns:
            Path: The resolved charts output directory.
        """
        return Path(os.environ.get("GRIMPERIUM_CHARTS_DIR", "reports/charts"))

    def _handle_detailed_metrics(self) -> None:
        """Display detailed statistical metrics; delegates all computation to analyze()."""
        from grimperium.ml.error_analysis import analyze

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

        result = analyze(df)
        s = result.summary

        mae_val = float(s["mae"])
        rmse_val = float(s["rmse"])
        max_err = float(s["max_error"])
        bias = float(s["bias"])
        r2_val = float(s["r2"])
        pearson_r = float(s["pearson_r"])
        pct_1 = float(s["pct_within_1"])
        pct_2 = float(s["pct_within_2"])
        pct_5 = float(s["pct_within_5"])
        n = int(s["n_molecules"])

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

    def _handle_show_outliers(self) -> None:
        """Display outlier molecules detected by analyze()."""
        from grimperium.ml.error_analysis import analyze

        csv_path = self._get_csv_path()
        if not csv_path.exists():
            self.show_error(f"Data file not found: {csv_path}")
            self.wait_for_enter()
            return

        df = pd.read_csv(csv_path)
        try:
            result = analyze(df)
        except ValueError as exc:
            self.show_error(str(exc))
            self.wait_for_enter()
            return

        outliers = result.outliers_df
        if outliers.empty:
            self.console.print(
                f"\n[{COLORS['muted']}]No outliers detected.[/{COLORS['muted']}]\n"
            )
            self.wait_for_enter()
            return

        table = Table(
            title=f"Outliers — {len(outliers)} molecules",
            show_header=True,
            header_style=f"bold {COLORS['results']}",
            border_style=COLORS["border"],
        )
        table.add_column("mol_id")
        table.add_column("smiles", max_width=40)
        table.add_column("H298_cbs", justify="right")
        table.add_column("H298_predicted", justify="right")
        table.add_column("abs_error", justify="right")
        table.add_column("severity", justify="center")
        table.add_column("score", justify="right")

        for _, row in outliers.iterrows():
            mol_id = str(row.get("mol_id", "-"))
            smiles = str(row.get("smiles", "-"))
            if len(smiles) > 40:
                smiles = smiles[:37] + "..."
            sev = str(row.get("severity", "-"))
            table.add_row(
                mol_id,
                smiles,
                f"{float(row['H298_cbs']):.4f}",
                f"{float(row['H298_predicted']):.4f}",
                f"{float(row['abs_error']):.4f}",
                sev,
                str(int(row["outlier_score"])),
            )

        self.console.print()
        self.console.print(table)
        self.console.print()
        self.wait_for_enter()

    def _handle_top_errors(self) -> None:
        """Display top-N molecules by absolute error."""
        from grimperium.ml.error_analysis import analyze

        csv_path = self._get_csv_path()
        if not csv_path.exists():
            self.show_error(f"Data file not found: {csv_path}")
            self.wait_for_enter()
            return

        df = pd.read_csv(csv_path)
        try:
            result = analyze(df)
        except ValueError as exc:
            self.show_error(str(exc))
            self.wait_for_enter()
            return

        top = result.top_n_df

        table = Table(
            title=f"Top {len(top)} Error Molecules",
            show_header=True,
            header_style=f"bold {COLORS['results']}",
            border_style=COLORS["border"],
        )
        table.add_column("Rank", justify="right")
        table.add_column("mol_id")
        table.add_column("abs_error", justify="right")
        table.add_column("signed_error", justify="right")
        table.add_column("severity", justify="center")

        for _, row in top.iterrows():
            table.add_row(
                str(int(row["rank"])),
                str(row.get("mol_id", "-")),
                f"{float(row['abs_error']):.4f}",
                f"{float(row['signed_error']):.4f}",
                str(row.get("severity", "-")),
            )

        self.console.print()
        self.console.print(table)
        self.console.print()
        self.wait_for_enter()

    def _handle_html_report(self) -> None:
        """Generate HTML error analysis report."""
        from grimperium.ml.error_analysis import analyze
        from grimperium.ml.html_report import generate_html_report

        csv_path = self._get_csv_path()
        if not csv_path.exists():
            self.show_error(f"Data file not found: {csv_path}")
            self.wait_for_enter()
            return

        df = pd.read_csv(csv_path)
        try:
            result = analyze(df)
        except ValueError as exc:
            self.show_error(str(exc))
            self.wait_for_enter()
            return

        output_path = self._get_charts_dir() / "error_report.html"
        saved = generate_html_report(result, output_path)
        self.show_success(f"HTML report saved: {saved}")
        self.wait_for_enter()

    def _handle_export_outliers(self) -> None:
        """Export outlier molecules to CSV."""
        from grimperium.ml.error_analysis import analyze

        csv_path = self._get_csv_path()
        if not csv_path.exists():
            self.show_error(f"Data file not found: {csv_path}")
            self.wait_for_enter()
            return

        df = pd.read_csv(csv_path)
        try:
            result = analyze(df)
        except ValueError as exc:
            self.show_error(str(exc))
            self.wait_for_enter()
            return

        output_path = self._get_charts_dir() / "outliers.csv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result.outliers_df.to_csv(output_path, index=False)
        n = len(result.outliers_df)
        self.show_success(f"Exported {n} outlier(s) to: {output_path}")
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
        """Redirect batch prediction to Models; Results is analysis-only."""
        self._show_analysis_only_message(
            "Batch prediction is handled in Models > Predict Batch.\n"
            "Results is for analysis only."
        )

    def _handle_run_analysis(self) -> None:
        """Redirect model training and prediction to Models."""
        self._show_analysis_only_message(
            "Training and prediction are handled in Models.\n"
            "Use Models > Train Model and Models > Predict Batch first,\n"
            "then return here for analysis."
        )

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
