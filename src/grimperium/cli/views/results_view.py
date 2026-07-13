"""
Results view for GRIMPERIUM CLI.

Displays performance analytics and divergence analysis using real model data.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from rich.panel import Panel
from rich.table import Table

from grimperium.cli.components import EmptyStatePanel, StatusBadge
from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.session import RunRef
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView
from grimperium.ml import charts as charts_module
from grimperium.ml import html_report as html_report_module
from grimperium.results.models import ResultsAnalysisMode, ResultsAnalysisReport
from grimperium.results.service import ResultsService


class ResultsView(BaseView):
    """View for performance analytics and results."""

    name = "results"
    title = "Results"
    icon = ICONS["results"]
    color = COLORS["results"]

    def __init__(self, controller: Any) -> None:
        super().__init__(controller)
        self.results_service = ResultsService()

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

    def _get_csv_path(self) -> Path | None:
        """Resolve CSV path from the typed live controller session."""
        session = getattr(self.controller, "__dict__", {}).get("session")
        csv_path = getattr(session, "analysis_path", None)
        if csv_path is None:
            csv_path = getattr(self.controller, "current_csv_path", None)
        if csv_path is None:
            return None
        return Path(csv_path)

    def _session_model_name(self) -> str | None:
        session = getattr(self.controller, "__dict__", {}).get("session")
        model = getattr(session, "model", None)
        model_name = getattr(model, "name", None)
        if model_name:
            return str(model_name)
        current_model = getattr(self.controller, "current_model", None)
        return str(current_model) if current_model else None

    def _active_run_id(self) -> str | None:
        session = getattr(self.controller, "__dict__", {}).get("session")
        run_ref = getattr(session, "run", None)
        run_id = getattr(run_ref, "run_id", None)
        return str(run_id) if run_id else None

    def _load_analysis_report(
        self, *, show_errors: bool
    ) -> ResultsAnalysisReport | None:
        """Load the active run or dataset analysis through ResultsService."""
        run_id = self._active_run_id()
        if run_id is not None:
            try:
                return self.results_service.analyze_run(run_id)
            except FileNotFoundError as exc:
                if show_errors:
                    self.show_error(f"Output ausente: {exc}")
                return None
            except ValueError as exc:
                if show_errors:
                    message = str(exc)
                    if "incompatible" in message.lower():
                        self.show_error(f"Run incompatível: {message}")
                    elif "inválido" in message.lower() or "invalid" in message.lower():
                        self.show_error(f"Arquivo inválido: {message}")
                    elif "refer" in message.lower():
                        self.show_error(f"Sem referência disponível: {message}")
                    else:
                        self.show_error(message)
                return None

        csv_path = self._get_csv_path()
        if csv_path is None:
            if show_errors:
                self.show_error(
                    "No analysis source selected. Choose a dataset, select a run, "
                    "or run a calculation first."
                )
            return None
        if not csv_path.exists():
            if show_errors:
                self.show_error(f"Output ausente: {csv_path}")
            return None
        try:
            return self.results_service.analyze_dataset(csv_path)
        except ValueError as exc:
            if show_errors:
                self.show_error(f"Arquivo inválido: {exc}")
            return None

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

        report = self._load_analysis_report(show_errors=False)
        if report is not None and report.model_label:
            table.add_row(
                report.model_label,
                "run manifest",
                "-",
                "-",
                StatusBadge(report.run_status or "Ready").text,
            )
            self.console.print(table)
            self.console.print()
            return
        if report is not None:
            method_id = report.scientific_summary.method_id
            if method_id in {"crest_pm7", "semiempirical_am1_pm3_pm7"}:
                table.add_row(
                    "Not required",
                    str(method_id),
                    "-",
                    "-",
                    StatusBadge(report.run_status or "Ready").text,
                )
                self.console.print(table)
                self.console.print()
                return

        metadata = self.results_service.model_metadata(
            self._get_model_path(),
            model_name=self._session_model_name(),
        )
        if metadata is not None:
            mae_val = metadata.get("mae")
            r2_val = metadata.get("r2")
            mae_str = f"{mae_val:.3f}" if mae_val is not None else "-"
            r2_str = f"{r2_val:.4f}" if r2_val is not None else "-"
            status_text = str(metadata.get("status", "Ready"))
            status = StatusBadge(status_text).text
            table.add_row(
                str(metadata.get("model_label", "Selected model")),
                str(metadata.get("algorithm", "model bundle")),
                mae_str,
                r2_str,
                status,
            )
        else:
            table.add_row(
                self._session_model_name() or "No model selected",
                "-",
                "-",
                "-",
                StatusBadge("Not selected").text,
            )

        self.console.print(table)
        self.console.print()

    def _render_divergence_analysis(self) -> None:
        """Render comparative divergence or scientific summary by analysis mode."""
        report = self._load_analysis_report(show_errors=False)

        if report is None:
            self.console.print(
                EmptyStatePanel(
                    title=f"[bold {COLORS['results']}]Divergence Analysis"
                    f"[/bold {COLORS['results']}]",
                    message=(
                        "No predictions available. Run Calculate, select a dataset, "
                        "or select a saved run before returning here for analysis."
                    ),
                    border_style=COLORS["border"],
                ).render()
            )
            self.console.print()
            return

        if report.analysis_mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY:
            self._render_scientific_summary(report)
            return

        title = (
            "PM7 Baseline vs Reference Divergence"
            if report.analysis_mode is ResultsAnalysisMode.BASELINE_WITH_REFERENCE
            else "Predicted vs CBS Divergence Distribution"
        )
        if report.scientific_summary.comparison_label:
            self.console.print(
                f"[{COLORS['muted']}]Comparison: "
                f"{report.scientific_summary.comparison_label}"
                f"[/{COLORS['muted']}]"
            )
        stats = report.divergence_distribution

        # Divergence distribution table
        table = Table(
            title=title,
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

        total_molecules = sum(stat.count for stat in stats)

        for stat in stats:
            color = severity_colors.get(stat.severity, COLORS["muted"])
            pct = max(0.0, min(stat.percentage, 100.0))
            bar_length = max(0, min(int(pct / 5), 20))
            bar = (
                f"[{color}]"
                f"{'█' * bar_length}{'░' * (20 - bar_length)}"
                f"[/{color}]"
            )

            range_max = stat.range_max
            range_str = (
                f"{stat.range_min:.0f}% - {range_max:.0f}%"
                if range_max is not None
                else f"{stat.range_min:.0f}%+"
            )

            table.add_row(
                f"[{color}]{stat.severity}[/{color}]",
                range_str,
                f"{stat.count:,}",
                f"{stat.percentage:.1f}%",
                bar,
            )

        self.console.print(table)
        self.console.print()

        # Summary panel
        low_medium = sum(
            stat.count for stat in stats if stat.severity in ["LOW", "MEDIUM"]
        )
        low_medium_pct = (
            (low_medium / total_molecules) * 100 if total_molecules > 0 else 0
        )

        if report.analysis_mode is ResultsAnalysisMode.BASELINE_WITH_REFERENCE:
            low_label = "PM7 baseline is close to the reference"
            medium_label = "Moderate PM7 baseline deviation"
        else:
            low_label = "Predictions are accurate, small corrections needed"
            medium_label = "Moderate corrections, ML performs well"

        summary = f"""
[bold]Key Findings:[/bold]

\u2022 Total molecules analyzed: {total_molecules:,}
\u2022 Molecules with LOW/MEDIUM divergence: {low_medium:,} ({low_medium_pct:.1f}%)
\u2022 Source analyzed: {report.source_label}

[bold]Interpretation:[/bold]

\u2022 [{COLORS['success']}]LOW (0-10%)[/{COLORS['success']}]: {low_label}
\u2022 [{COLORS['warning']}]MEDIUM (10-25%)[/{COLORS['warning']}]: {medium_label}
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

    def _render_scientific_summary(self, report: ResultsAnalysisReport) -> None:
        """Render a non-comparative scientific overview for baseline-only runs."""
        sci = report.scientific_summary
        table = Table(
            title="Scientific Run Summary",
            show_header=True,
            header_style=f"bold {COLORS['results']}",
            border_style=COLORS["border"],
        )
        table.add_column("Field", style="bold")
        table.add_column("Value")
        table.add_row("Molecules", str(sci.n_molecules))
        table.add_row("Estimates", str(sci.n_estimates))
        table.add_row("Roles", ", ".join(sci.roles) or "-")
        table.add_row("Hamiltonians", ", ".join(sci.hamiltonians) or "-")
        table.add_row(
            "Value range (kcal/mol)",
            (
                f"{sci.value_min:.3f} … {sci.value_max:.3f}"
                if sci.value_min is not None and sci.value_max is not None
                else "-"
            ),
        )
        table.add_row(
            "Mean / median",
            (
                f"{sci.value_mean:.3f} / {sci.value_median:.3f}"
                if sci.value_mean is not None and sci.value_median is not None
                else "-"
            ),
        )
        table.add_row("Method", sci.method_id or report.method_label or "-")
        table.add_row("Run status", sci.run_status or report.run_status or "-")
        table.add_row("Started", sci.started_at or "-")
        table.add_row("Completed", sci.completed_at or "-")
        table.add_row(
            "Duration (s)",
            f"{sci.duration_seconds:.1f}" if sci.duration_seconds is not None else "-",
        )
        table.add_row("Comparative metrics", "not available (no reference)")
        self.console.print(table)
        if sci.warnings:
            for warning in sci.warnings:
                self.console.print(
                    f"[{COLORS['warning']}]{warning}[/{COLORS['warning']}]"
                )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the results view."""
        return [
            MenuOption(
                label="Analyze Active Source",
                value="analyze_source",
                icon=ICONS.get("results", "\U0001f4ca"),
            ),
            MenuOption(
                label="Select Saved Run",
                value="select_run",
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

        if action == "analyze_source":
            self._handle_analyze_source()
            return None

        if action == "select_run":
            self._handle_select_run()
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
        report = self._load_analysis_report(show_errors=True)
        if report is None:
            self.wait_for_enter()
            return

        if not report.has_comparative_metrics:
            self._render_scientific_summary(report)
            self.wait_for_enter()
            return

        s = report.summary

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
            f"{pearson_r:.5f}" if pearson_r == pearson_r else "N/A",
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
        report = self._load_analysis_report(show_errors=True)
        if report is None:
            self.wait_for_enter()
            return

        outliers = report.outliers_df
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
        report = self._load_analysis_report(show_errors=True)
        if report is None:
            self.wait_for_enter()
            return

        top = report.top_n_df

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
        report = self._load_analysis_report(show_errors=True)
        if report is None:
            self.wait_for_enter()
            return
        if report.analysis is None:
            self.show_error(
                "HTML report requires comparative metrics "
                f"(mode={report.analysis_mode.value})."
            )
            self.wait_for_enter()
            return

        output_path = self._get_charts_dir() / "error_report.html"
        saved = html_report_module.generate_html_report(report.analysis, output_path)
        self.show_success(f"HTML report saved: {saved}")
        self.wait_for_enter()

    def _handle_export_outliers(self) -> None:
        """Export outlier molecules to CSV."""
        report = self._load_analysis_report(show_errors=True)
        if report is None:
            self.wait_for_enter()
            return

        output_path = self._get_charts_dir() / "outliers.csv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        report.outliers_df.to_csv(output_path, index=False)
        n = len(report.outliers_df)
        self.show_success(f"Exported {n} outlier(s) to: {output_path}")
        self.wait_for_enter()

    def _handle_charts(self) -> None:
        """Generate visualization charts and display results."""
        report = self._load_analysis_report(show_errors=True)
        if report is None:
            self.wait_for_enter()
            return
        if report.analysis_mode is ResultsAnalysisMode.SCIENTIFIC_SUMMARY_ONLY:
            self.console.print(
                f"[{COLORS['warning']}]Comparative charts are not available for "
                "scientific_summary_only runs (no reference/prediction pair)."
                f"[/{COLORS['warning']}]"
            )
            self.wait_for_enter()
            return
        scored_df = report.scored_df
        if scored_df.empty:
            self.show_error("No comparative rows are available for chart generation.")
            self.wait_for_enter()
            return
        charts_dir = self._get_charts_dir()
        charts_dir.mkdir(parents=True, exist_ok=True)
        chart_input = charts_dir / "chart_input.csv"
        scored_df.to_csv(chart_input, index=False)

        self.console.print()
        self.console.print(
            f"[bold {COLORS['results']}]Generating charts..."
            f"[/bold {COLORS['results']}]"
        )
        self.console.print()

        try:
            result = charts_module.generate_charts(chart_input, charts_dir)

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

    def _handle_analyze_source(self) -> None:
        """Re-render analysis for the active run or dataset."""
        report = self._load_analysis_report(show_errors=True)
        if report is None:
            self.wait_for_enter()
            return
        self.clear_screen()
        self.show_header()
        self._render_model_comparison()
        self._render_divergence_analysis()
        self.wait_for_enter()

    def _handle_select_run(self) -> None:
        """Let the user pick a persisted run as the active analysis source."""
        run_service = getattr(self.controller, "run_service", None)
        if run_service is None:
            self.show_error("Run service unavailable.")
            self.wait_for_enter()
            return

        manifests = run_service.list_runs()
        if not manifests:
            self.console.print(
                EmptyStatePanel(
                    title="Saved Runs",
                    message=(
                        "No saved runs found. Execute Calculate or a batch job "
                        "first, then return here to analyze the run outputs."
                    ),
                    border_style=COLORS["border"],
                ).render()
            )
            self.wait_for_enter()
            return

        options = [
            MenuOption(
                label=f"{manifest.run_id} [{manifest.status.value}]",
                value=manifest.run_id,
                icon=ICONS.get("results", "\U0001f4ca"),
            )
            for manifest in manifests
        ]
        selected = show_back_menu(options=options, title="Select Saved Run")
        if selected is None or selected == "back":
            return

        manifest = run_service.get_run(selected)
        run_ref = RunRef(run_id=manifest.run_id, status=manifest.status.value)
        set_run = getattr(self.controller, "set_run", None)
        if callable(set_run):
            set_run(run_ref)
        else:
            self.controller.session.run = run_ref
        self.show_success(f"Active run set to {manifest.run_id}")
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
