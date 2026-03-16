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

    def _get_model_path(self) -> Path:
        """Resolve model path from env var (read fresh on each call)."""
        return Path(
            os.environ.get("GRIMPERIUM_MODEL_PATH", "models/delta_learner_v1.joblib")
        )

    def _get_csv_path(self) -> Path:
        """Resolve CSV path from env var (read fresh on each call)."""
        return Path(os.environ.get("GRIMPERIUM_DATA_PATH", "data/thermo_pm7.csv"))

    def _load_real_model_row(self) -> dict[str, Any] | None:
        """Load model metadata for display. Returns None if not available."""
        try:
            return load_model_metadata(self._get_model_path())
        except FileNotFoundError:
            return None

    def _compute_divergence_stats(self) -> list[DictStrAny] | None:
        """Compute divergence stats from CSV with predictions.

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
                label="Predict Batch",
                value="predict_batch",
                icon=ICONS.get("calc", "\u26a1"),
            ),
            MenuOption(
                label="Detailed Metrics",
                value="detailed",
                disabled=True,
                disabled_reason="In Development",
            ),
            MenuOption(
                label="Visualization Charts",
                value="charts",
                disabled=True,
                disabled_reason="In Development",
            ),
        ]

    def handle_action(self, action: str) -> str | None:
        """Handle menu actions."""
        if action == "back":
            return "main"

        if action == "predict_batch":
            self._handle_predict_batch()
            return None

        # Handle in-development features
        if action in ["detailed", "charts"]:
            self.show_in_development(action.title())
            return None

        return None

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
