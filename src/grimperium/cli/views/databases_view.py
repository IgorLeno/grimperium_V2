"""
Databases view for GRIMPERIUM CLI.

Displays and manages molecular databases.
"""

from __future__ import annotations

import logging
import multiprocessing
import threading
import time
from datetime import date, datetime
from pathlib import Path
from queue import Queue
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from rich.live import Live
from rich.panel import Panel
from rich.table import Table

from grimperium.cli.constants import DATA_DIR, PHASE_A_RESULTS_FILE
from grimperium.cli.database_registry import DatabaseInfo, DatabaseRegistry
from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.menu import confirm as menu_confirm
from grimperium.cli.progress_tracker import (
    CSVMonitor,
    ProgressEvent,
    ProgressTracker,
    consume_events,
)
from grimperium.cli.session_store import (
    SessionState,
    WorkerSessionInfo,
    delete_session,
    load_session,
    save_session,
)
from grimperium.cli.settings_manager import DistributedDefaults, SettingsManager
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView
from grimperium.core.metrics import mae, r2_score
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus, WorkerStatus
from grimperium.crest_pm7.database_analyzer import AnalysisReport
from grimperium.worker.client import WorkerClient, WorkerClientConfig
from grimperium.worker.runner import WorkerConfig, WorkerRunner

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController

logger = logging.getLogger(__name__)


def _compute_pm7_stats(valid: pd.DataFrame) -> dict[str, float]:
    """Compute PM7 vs CBS absolute-error statistics from a validated DataFrame.

    Args:
        valid: DataFrame with non-null ``H298_cbs`` and ``H298_pm7`` columns.

    Returns:
        Dict with keys: mare, bias, r2, n, p50, p90, p95,
        pct_lt_1, pct_lt_2, pct_lt_5.

    Note:
        Relative-error metrics (MRE%, MdRE%) are intentionally absent.
        H298 heats of formation cross zero, making percentage errors
        numerically undefined for molecules with small |H298_cbs|.
    """
    h298_cbs = valid["H298_cbs"].to_numpy(dtype=float)
    h298_pm7 = valid["H298_pm7"].to_numpy(dtype=float)

    abs_err = np.abs(h298_pm7 - h298_cbs)
    n = len(abs_err)

    return {
        "mare": float(mae(h298_cbs, h298_pm7)),
        "bias": float(np.mean(h298_pm7 - h298_cbs)),
        "r2": float(r2_score(h298_cbs, h298_pm7)),
        "n": float(n),
        "p50": float(np.median(abs_err)),
        "p90": float(np.percentile(abs_err, 90)),
        "p95": float(np.percentile(abs_err, 95)),
        "pct_lt_1": float(np.sum(abs_err < 1.0) / n * 100),
        "pct_lt_2": float(np.sum(abs_err < 2.0) / n * 100),
        "pct_lt_5": float(np.sum(abs_err < 5.0) / n * 100),
    }


def _filter_suspect_rows(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Remove rows flagged as SUSPECT CBS reference values.

    Args:
        df: DataFrame that may contain a ``cbs_quality_flag`` column.

    Returns:
        Tuple of (filtered_df, suspect_count). If the column is absent,
        returns (df, 0) unchanged — backward-compatible with datasets that
        predate the flag.
    """
    if "cbs_quality_flag" not in df.columns:
        return df, 0
    suspect_mask = df["cbs_quality_flag"] == "SUSPECT"
    suspect_count: int = int(suspect_mask.sum())
    return df[~suspect_mask].copy(), suspect_count


class DatabasesView(BaseView):
    """View for managing molecular databases."""

    name = "databases"
    title = "Databases"
    icon = ICONS["databases"]
    color = COLORS["databases"]

    def __init__(self, controller: CliController) -> None:
        """Initialize the databases view."""
        super().__init__(controller)
        self.registry = DatabaseRegistry(DATA_DIR)
        self.selected_db: DatabaseInfo | None = None
        self._server_proc: multiprocessing.Process | None = None

    def get_databases(self) -> list[DatabaseInfo]:
        """Get list of databases from the registry.

        Returns:
            List of DatabaseInfo objects enriched with live CSV data.
        """
        return self.registry.load()

    def _format_last_updated(self, value: datetime | date | None) -> str:
        """Format database last_updated timestamp consistently.

        Args:
            value: A datetime, date, or None.

        Returns:
            Formatted string or "-" if value is None.

        Note:
            If the Database.last_updated is timezone-aware, the timezone info
            is preserved in the formatted output.
        """
        if value is None:
            return "-"
        if isinstance(value, datetime):
            return value.strftime("%Y-%m-%d %H:%M:%S")
        if isinstance(value, date):
            return value.strftime("%Y-%m-%d")
        return str(value)

    def render(self) -> None:
        """Render the databases overview."""
        self.clear_screen()
        self.show_header()

        # Databases table
        table = Table(
            title="Available Databases",
            show_header=True,
            header_style=f"bold {COLORS['databases']}",
            border_style=COLORS["border"],
        )
        table.add_column("Name", style="bold")
        table.add_column("Molecules", justify="right")
        table.add_column("Last Updated")
        table.add_column("Status")

        for db in self.get_databases():
            if db.status == "ready":
                status = f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
            else:
                status = (
                    f"[{COLORS['in_dev']}]{ICONS['in_dev']} In Dev[/{COLORS['in_dev']}]"
                )

            last_updated_str = self._format_last_updated(db.last_updated)

            table.add_row(
                f"{db.alias} ({db.name})",
                f"{db.molecules:,}" if db.molecules > 0 else "-",
                last_updated_str,
                status,
            )

        self.console.print(table)
        self.console.print()

    def render_database_detail(self, db: DatabaseInfo) -> None:
        """Render detailed view for a specific database."""
        self.clear_screen()
        self.show_header()

        # Database info panel
        status_text = (
            f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
            if db.status == "ready"
            else f"[{COLORS['in_dev']}]{ICONS['in_dev']} In Development[/{COLORS['in_dev']}]"
        )

        last_updated_str = self._format_last_updated(db.last_updated)

        info = f"""
[bold]Name:[/bold]         {db.name}
[bold]Alias:[/bold]        {db.alias}
[bold]Description:[/bold]  {db.description}
[bold]Pipeline:[/bold]     {db.pipeline}
[bold]Molecules:[/bold]    {db.molecules:,}
[bold]Last Updated:[/bold] {last_updated_str}
[bold]Status:[/bold]       {status_text}

[bold]Properties:[/bold]
"""
        for prop in db.properties:
            info += f"  • {prop}\n"

        self.console.print(
            Panel(
                info,
                title=f"[bold {COLORS['databases']}]{db.alias} — {db.name}[/bold {COLORS['databases']}]",
                border_style=COLORS["databases"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the databases view."""
        options = []
        for db in self.get_databases():
            options.append(
                MenuOption(
                    label=db.alias,
                    value=f"view_{db.alias}",
                    icon=ICONS["databases"],
                )
            )

        options.extend(
            [
                MenuOption(
                    label="Refresh Databases",
                    value="refresh",
                    icon="",
                    description="Re-scan CSV files and reload registry",
                ),
                MenuOption(
                    label="Add New Database",
                    value="add",
                    disabled=True,
                    disabled_reason="In Development",
                ),
            ]
        )

        return options

    def get_detail_menu_options(self) -> list[MenuOption]:
        """Return menu options for database detail view."""
        options = [
            MenuOption(
                label="Run Calculation",
                value="calculate_run",
            ),
            MenuOption(
                label="Configuration",
                value="calculate_config",
            ),
            MenuOption(
                label="Analyze Database",
                value="analyze",
            ),
        ]

        # Conditionally add PM7 Baseline Analysis if CSV has required columns
        if self.selected_db and self.selected_db.csv_path:
            csv_path = DATA_DIR / self.selected_db.csv_path
            if csv_path.exists():
                try:
                    csv_cols = pd.read_csv(csv_path, nrows=0).columns
                    if "H298_pm7" in csv_cols and "H298_cbs" in csv_cols:
                        options.append(
                            MenuOption(
                                label="PM7 Baseline Analysis",
                                value="pm7_baseline",
                            )
                        )
                except Exception as exc:
                    logger.debug("Could not check CSV columns for PM7 menu: %s", exc)

        options.extend(
            [
                MenuOption(
                    label="Edit Database",
                    value="edit",
                    disabled=True,
                    disabled_reason="In Development",
                ),
                MenuOption(
                    label="Delete Database",
                    value="delete",
                    disabled=True,
                    disabled_reason="In Development",
                ),
            ]
        )
        return options

    def handle_action(self, action: str) -> str | None:
        """Handle menu actions."""
        if action == "back":
            if self.selected_db:
                self.selected_db = None
                return None  # Stay in databases view
            return "main"

        if action and action.startswith("view_"):
            alias = action.removeprefix("view_")
            self.selected_db = self.registry.get_by_alias(alias)
            return None

        if action == "calculate_run":
            self._handle_run_calculation_menu()
            return None

        if action == "calculate_config":
            self._handle_calculate_config()
            return None

        if action == "analyze":
            self._handle_analyze()
            return None

        if action == "pm7_baseline":
            self._handle_pm7_baseline()
            return None

        if action == "refresh":
            self.registry.reload()
            self.refresh_databases_from_filesystem()
            self.console.input("[dim]Press Enter to continue...[/dim]")
            return None

        if action in ["add", "edit", "delete"]:
            self.show_in_development(action.title())
            return None

        return None

    def _handle_calculate_config(self) -> None:
        """Show CREST/MOPAC configuration submenus."""
        settings_manager: SettingsManager | None = getattr(
            self.controller, "settings_manager", None
        )
        if settings_manager is None:
            self.console.print(
                f"[{COLORS['warning']}]Settings manager not available.[/{COLORS['warning']}]"
            )
            self.wait_for_enter()
            return

        config_options = [
            MenuOption(label="CREST Settings", value="crest_config"),
            MenuOption(label="MOPAC Settings", value="mopac_config"),
        ]
        result = show_back_menu(options=config_options, title="Configuration")

        if result == "crest_config":
            settings_manager.display_crest_menu()
        elif result == "mopac_config":
            settings_manager.display_mopac_menu()

    def _handle_analyze(self) -> None:
        """Run database analysis and display the report."""
        if self.selected_db is None:
            return

        csv_path = (
            DATA_DIR / self.selected_db.csv_path if self.selected_db.csv_path else None
        )
        if csv_path is None or not csv_path.exists():
            self.console.print(
                f"[{COLORS['warning']}]No CSV file available for {self.selected_db.alias}.[/{COLORS['warning']}]"
            )
            self.wait_for_enter()
            return

        from grimperium.crest_pm7.database_analyzer import DatabaseAnalyzer

        detail_path = DATA_DIR / "molecules_pm7" / "conformer_details"
        detail_dir: Path | None = detail_path if detail_path.exists() else None

        self.console.print(
            f"\n[bold {COLORS['databases']}]Analyzing {self.selected_db.alias}...[/bold {COLORS['databases']}]"
        )

        analyzer = DatabaseAnalyzer()
        report = analyzer.analyze(csv_path, detail_dir=detail_dir)
        self.render_analysis_report(report)

    def _handle_pm7_baseline(self) -> None:
        """Compute and display PM7 baseline error analysis vs CBS reference."""
        if self.selected_db is None:
            return

        csv_path = (
            DATA_DIR / self.selected_db.csv_path if self.selected_db.csv_path else None
        )
        if csv_path is None or not csv_path.exists():
            self.show_error(f"CSV file not found for {self.selected_db.alias}.")
            self.wait_for_enter()
            return

        df = pd.read_csv(csv_path)
        if "H298_pm7" not in df.columns or "H298_cbs" not in df.columns:
            self.show_error("CSV must contain both 'H298_pm7' and 'H298_cbs' columns.")
            self.wait_for_enter()
            return

        valid = df.dropna(subset=["H298_pm7", "H298_cbs"])
        valid, suspect_count = _filter_suspect_rows(valid)
        if len(valid) == 0:
            self.show_error("No valid rows with both PM7 and CBS values.")
            self.wait_for_enter()
            return

        stats = _compute_pm7_stats(valid)

        # Display table
        db_color = COLORS["databases"]
        table = Table(
            title="PM7 Baseline Analysis",
            show_header=True,
            header_style=f"bold {db_color}",
            border_style=COLORS["border"],
            show_edge=False,
        )
        table.add_column("Metric", style="bold")
        table.add_column("Value", justify="right")

        table.add_row("MARE (kcal/mol)", f"{stats['mare']:.4f}")
        table.add_row("Bias (kcal/mol)", f"{stats['bias']:.4f}")
        table.add_row("R\u00b2", f"{stats['r2']:.4f}")
        table.add_row("", "")
        table.add_row("Absolute Error Distribution", "")
        table.add_row("  Median |error| (P50)", f"{stats['p50']:.3f} kcal/mol")
        table.add_row("  P90", f"{stats['p90']:.3f} kcal/mol")
        table.add_row("  P95", f"{stats['p95']:.3f} kcal/mol")
        table.add_row("", "")
        table.add_row("  |error| < 1 kcal/mol", f"{stats['pct_lt_1']:.1f}%")
        table.add_row("  |error| < 2 kcal/mol", f"{stats['pct_lt_2']:.1f}%")
        table.add_row("  |error| < 5 kcal/mol", f"{stats['pct_lt_5']:.1f}%")
        table.add_row("", "")
        table.add_row("Molecules analyzed", f"{stats['n']:,}")

        self.console.print()
        self.console.print(table)
        if suspect_count > 0:
            self.console.print(
                f"[yellow]\u26a0 {suspect_count} molecule(s) with cbs_quality_flag=SUSPECT excluded from this analysis. "
                "These values originate from the CBS source dataset and are "
                "flagged as unreliable. Retained in CSV for traceability.[/yellow]"
            )
        self.console.print()

        # Interpretation panel
        interpretation = (
            f"[bold]PM7 Baseline Error Context:[/bold]\n\n"
            f"\u2022 PM7 mean absolute error: {stats['mare']:.4f} kcal/mol vs CBS reference\n"
            f"\u2022 Systematic bias: {stats['bias']:+.4f} kcal/mol "
            f"({'PM7 overestimates' if stats['bias'] > 0 else 'PM7 underestimates' if stats['bias'] < 0 else 'unbiased'})\n"
            f"\u2022 Median absolute error (P50): {stats['p50']:.3f} kcal/mol\n"
            f"\u2022 This is the error that delta-learning needs to correct"
        )
        self.console.print(
            Panel(
                interpretation,
                title=f"[bold {db_color}]Interpretation[/bold {db_color}]",
                border_style=COLORS["border"],
                padding=(1, 2),
            )
        )
        self.console.print()
        self.wait_for_enter()

    def render_analysis_report(self, report: AnalysisReport) -> None:
        """Render a full analysis report with Rich panels.

        Args:
            report: AnalysisReport from DatabaseAnalyzer.
        """
        self.clear_screen()
        self.show_header()

        db_color = COLORS["databases"]

        # ── Panel 1: STATUS OVERVIEW ──
        status_lines = [f"[bold]Total rows:[/bold] {report.total_rows:,}"]
        status_lines.append(f"[bold]OK:[/bold] {report.ok_count:,}")
        for status_name, count in sorted(report.status_counts.items()):
            if status_name != "OK":
                pct = count / report.total_rows * 100 if report.total_rows else 0
                status_lines.append(f"  {status_name}: {count:,} ({pct:.1f}%)")
        if report.orphan_running > 0:
            status_lines.append(
                f"[yellow]Orphaned RUNNING: {report.orphan_running}[/yellow]"
            )
        if report.cbs_ok_count > 0 or report.cbs_suspect_count > 0:
            total_flagged = report.cbs_ok_count + report.cbs_suspect_count
            ok_pct = report.cbs_ok_count / total_flagged * 100 if total_flagged else 0
            status_lines.append("")
            status_lines.append("[bold]CBS Quality:[/bold]")
            status_lines.append(
                f"  Available for training: {report.cbs_ok_count:,} ({ok_pct:.1f}%)"
            )
            if report.cbs_suspect_count > 0:
                status_lines.append(
                    f"  [yellow]SUSPECT (excluded): {report.cbs_suspect_count}[/yellow]"
                )
        self.console.print(
            Panel(
                "\n".join(status_lines),
                title=f"[bold {db_color}]Status Overview[/bold {db_color}]",
                border_style=db_color,
            )
        )

        # ── Panel 2: QUALITY GRADES ──
        if report.grade_distribution:
            grade_lines = []
            total_graded = sum(report.grade_distribution.values())
            bar_width = 30
            for grade in ["A", "B", "C", "FAILED"]:
                count = report.grade_distribution.get(grade, 0)
                pct = count / total_graded if total_graded else 0
                filled = int(pct * bar_width)
                bar = "\u2588" * filled + "\u2591" * (bar_width - filled)
                grade_lines.append(f"  {grade:>6}: {bar} {count:,} ({pct:.1%})")
            threshold_color = (
                COLORS["success"] if report.grade_ab_rate >= 0.67 else COLORS["warning"]
            )
            grade_lines.append(
                f"\n[{threshold_color}]A+B rate: {report.grade_ab_rate:.1%} "
                f"(threshold: 67%)[/{threshold_color}]"
            )
            self.console.print(
                Panel(
                    "\n".join(grade_lines),
                    title=f"[bold {db_color}]Quality Grades[/bold {db_color}]",
                    border_style=db_color,
                )
            )

        # ── Panel 3: H298_PM7 ──
        h298_lines = []
        if report.h298_pm7_mean is not None:
            h298_lines.append(
                f"Mean: {report.h298_pm7_mean:.2f}  Std: {report.h298_pm7_std:.2f}"
            )
            h298_lines.append(
                f"Min:  {report.h298_pm7_min:.2f}  Max: {report.h298_pm7_max:.2f}"
            )
            if report.h298_pm7_positive_count > 0:
                h298_lines.append(
                    f"[yellow]Positive H298: {report.h298_pm7_positive_count}[/yellow]"
                )
            if report.h298_pm7_extreme_count > 0:
                h298_lines.append(
                    f"[yellow]Extreme (per nheavy): {report.h298_pm7_extreme_count}[/yellow]"
                )
            if report.null_h298_pm7 > 0:
                h298_lines.append(
                    f"[yellow]Null in OK rows: {report.null_h298_pm7}[/yellow]"
                )
        else:
            h298_lines.append("[dim]No H298_pm7 data available[/dim]")
        self.console.print(
            Panel(
                "\n".join(h298_lines),
                title=f"[bold {db_color}]H298_PM7[/bold {db_color}]",
                border_style=db_color,
            )
        )

        # ── Panel 4: TARGET DELTA ──
        delta_lines = []
        if report.target_delta_mean is not None:
            delta_lines.append(
                f"Mean: {report.target_delta_mean:.2f}  Std: {report.target_delta_std:.2f}"
            )
            delta_lines.append(
                f"Min:  {report.target_delta_min:.2f}  Max: {report.target_delta_max:.2f}"
            )
            delta_lines.append(
                f"Outliers >200: {report.target_delta_outliers_200}  "
                f">500: {report.target_delta_outliers_500}"
            )
            if report.top_delta_outliers:
                delta_lines.append("\n[bold]Top outliers:[/bold]")
                for entry in report.top_delta_outliers[:10]:
                    mol_id = entry.get("mol_id", "?")
                    delta_val = entry.get("delta", 0)
                    smiles = entry.get("smiles", "")
                    nheavy = entry.get("nheavy", "?")
                    delta_lines.append(
                        f"  {mol_id}: {delta_val:+.1f} kcal/mol "
                        f"(nheavy={nheavy}) {smiles}"
                    )
        else:
            delta_lines.append("[dim]No target_delta data available[/dim]")
        self.console.print(
            Panel(
                "\n".join(delta_lines),
                title=f"[bold {db_color}]Target Delta[/bold {db_color}]",
                border_style=db_color,
            )
        )

        # ── Panel 5: ELECTRONIC DESCRIPTORS ──
        elec_lines = []
        if report.homo_mean is not None:
            elec_lines.append(
                f"HOMO mean: {report.homo_mean:.2f} eV  std: {report.homo_std:.2f}"
            )
        if report.lumo_mean is not None:
            elec_lines.append(f"LUMO mean: {report.lumo_mean:.2f} eV")
        if report.gap_mean is not None:
            elec_lines.append(f"Gap mean:  {report.gap_mean:.2f} eV")
        if report.homo_suspicious > 0:
            elec_lines.append(
                f"[yellow]|HOMO| > 100 eV: {report.homo_suspicious}[/yellow]"
            )
        if report.gap_negative_count > 0:
            elec_lines.append(f"[red]Negative gap: {report.gap_negative_count}[/red]")
        if report.gap_suspicious > 0:
            elec_lines.append(f"[yellow]Gap > 100 eV: {report.gap_suspicious}[/yellow]")
        if not elec_lines:
            elec_lines.append("[dim]No electronic descriptor data available[/dim]")
        self.console.print(
            Panel(
                "\n".join(elec_lines),
                title=f"[bold {db_color}]Electronic Descriptors[/bold {db_color}]",
                border_style=db_color,
            )
        )

        # ── Panel 6: PIPELINE HEALTH ──
        pipeline_lines = []
        if report.reruns_distribution:
            pipeline_lines.append("[bold]Reruns distribution:[/bold]")
            for reruns, count in sorted(report.reruns_distribution.items()):
                pipeline_lines.append(f"  reruns={reruns}: {count:,}")
        if report.max_reruns_count > 0:
            pipeline_lines.append(
                f"[yellow]Max reruns (=3): {report.max_reruns_count}[/yellow]"
            )
        if report.crest_failed_count > 0:
            pipeline_lines.append(f"CREST failed: {report.crest_failed_count}")
        if report.mopac_failed_count > 0:
            pipeline_lines.append(f"MOPAC failed: {report.mopac_failed_count}")
        if report.crest_timeout_count > 0:
            pipeline_lines.append(
                f"[yellow]CREST timeouts: {report.crest_timeout_count}[/yellow]"
            )
        if report.mopac_timeout_count > 0:
            pipeline_lines.append(
                f"[yellow]MOPAC timeouts: {report.mopac_timeout_count}[/yellow]"
            )
        if report.conformers_mean is not None:
            pipeline_lines.append(
                f"Conformers mean: {report.conformers_mean:.1f}  "
                f"zero: {report.conformers_zero_count}  "
                f"one: {report.conformers_one_count}"
            )
        if report.crest_time_mean is not None:
            pipeline_lines.append(
                f"CREST time mean: {report.crest_time_mean:.0f}s  "
                f"max: {report.crest_time_max:.0f}s"
            )
        if not pipeline_lines:
            pipeline_lines.append("[dim]No pipeline health data available[/dim]")
        self.console.print(
            Panel(
                "\n".join(pipeline_lines),
                title=f"[bold {db_color}]Pipeline Health[/bold {db_color}]",
                border_style=db_color,
            )
        )

        # ── Panel 7: ALERTS ──
        if report.alerts:
            alert_lines = []
            for alert in report.alerts:
                alert_lines.append(alert)
            self.console.print(
                Panel(
                    "\n".join(alert_lines),
                    title="[bold yellow]Alerts[/bold yellow]",
                    border_style="yellow",
                )
            )
        else:
            self.console.print(
                Panel(
                    "[green]No alerts — database looks healthy![/green]",
                    title="[bold green]Alerts[/bold green]",
                    border_style="green",
                )
            )

        self.console.print()
        self.wait_for_enter()

    def _refresh_database(self) -> None:
        """Trigger database refresh from source.

        Resyncs the working CSV from the source-of-truth CSV,
        resetting all status fields to PENDING.
        """
        from grimperium.cli.dataset_manager import DatasetManager

        manager = DatasetManager(
            source_csv=DATA_DIR / "thermo_cbs_chon.csv",
            working_csv=DATA_DIR / "thermo_pm7.csv",
        )
        try:
            manager.refresh_database()
            self.console.print("[green]Database refreshed[/green]")
        except FileNotFoundError as e:
            self.console.print(f"[red]Refresh failed: {e}[/red]")
        except Exception as e:
            self.console.print(f"[red]Refresh failed: {e}[/red]")
        self.console.input("[dim]Press Enter to continue...[/dim]")

    def refresh_databases_from_filesystem(self) -> int:
        """Scan data/ directory for CSV files and display database info.

        Discovers available database CSV files in the data/ directory,
        counts molecules in each file, and displays summary.

        Returns:
            Number of database CSV files discovered.
        """
        try:
            if not DATA_DIR.exists():
                self.console.print(
                    f"[{COLORS['warning']}]⚠️  Data directory not found: {DATA_DIR}[/{COLORS['warning']}]"
                )
                return 0

            # Find all CSV files in data/ directory
            csv_files = list(DATA_DIR.glob("*.csv"))

            if not csv_files:
                self.console.print(
                    f"[{COLORS['muted']}]No CSV files found in {DATA_DIR}[/{COLORS['muted']}]"
                )
                return 0

            self.console.print()
            self.console.print(
                f"[bold {COLORS['databases']}]Discovered {len(csv_files)} database(s):[/bold {COLORS['databases']}]"
            )
            self.console.print()

            # Display each database with molecule count
            for csv_file in sorted(csv_files):
                try:
                    df = pd.read_csv(csv_file)

                    # ✅ FIX: Show breakdown of calculated (OK) vs pending
                    if "status" in df.columns:
                        ok_count = len(df[df["status"] == "OK"])
                        total_count = len(df)
                        pending_count = total_count - ok_count

                        if ok_count > 0:
                            msg = f"{ok_count} calculated"
                            if pending_count > 0:
                                msg += f" ({pending_count} pending)"
                        else:
                            msg = f"{total_count} pending (none calculated)"
                    else:
                        # Legacy CSV without status - just show total
                        msg = f"{len(df)} molecules"

                    self.console.print(
                        f"  [{COLORS['success']}]✓[/{COLORS['success']}] "
                        f"[bold]{csv_file.name}[/bold]: "
                        f"[{COLORS['databases']}]{msg}[/{COLORS['databases']}]"
                    )
                except Exception as e:
                    self.console.print(
                        f"  [{COLORS['error']}]✗[/{COLORS['error']}] "
                        f"[bold]{csv_file.name}[/bold]: "
                        f"[{COLORS['muted']}]Error reading file: {e}[/{COLORS['muted']}]"
                    )

            self.console.print()
            self.console.print(
                f"[{COLORS['success']}]✓ Refresh complete: {len(csv_files)} database(s) found[/{COLORS['success']}]"
            )

            return len(csv_files)

        except Exception as e:
            self.console.print(
                f"[{COLORS['error']}]❌ Refresh failed: {e}[/{COLORS['error']}]"
            )
            return 0

    def handle_calculate_pm7(self) -> None:
        """Interactive handler for CREST PM7 batch calculations.

        Collects user inputs, runs batch processor with progress bar,
        updates results and returns to database view.
        """
        self.console.print()
        self.console.print(
            Panel(
                "[bold cyan]CREST PM7 Batch Calculation Configuration[/bold cyan]",
                border_style="cyan",
            )
        )
        self.console.print()

        try:
            num_molecules = self._prompt_positive_int(
                "How many molecules to calculate?",
                default=10,
                max_value=30026,
            )
            if num_molecules is None:
                return

            crest_timeout = self._prompt_positive_int(
                "CREST timeout per molecule (minutes)?",
                default=30,
            )
            if crest_timeout is None:
                return

            mopac_timeout = self._prompt_positive_int(
                "MOPAC/PM7 timeout per molecule (minutes)?",
                default=60,
            )
            if mopac_timeout is None:
                return

            self.console.print()
            self.console.print("[bold]Configuration Summary:[/bold]")
            self.console.print(f"  • Molecules to calculate: {num_molecules}")
            self.console.print(f"  • CREST timeout: {crest_timeout} min")
            self.console.print(f"  • MOPAC timeout: {mopac_timeout} min")
            self.console.print()

            confirm = (
                self.console.input(
                    "[yellow]Proceed with calculation? (yes/no) [/yellow]"
                )
                .strip()
                .lower()
            )
            if confirm not in ("yes", "y"):
                self.console.print("[yellow]Calculation cancelled.[/yellow]")
                self.console.input("[dim]Press Enter to continue...[/dim]")
                return

            self._run_pm7_batch(
                num_molecules=num_molecules,
                crest_timeout_minutes=crest_timeout,
                mopac_timeout_minutes=mopac_timeout,
            )

        except KeyboardInterrupt:
            self.console.print("\n[yellow]Calculation interrupted.[/yellow]")
            self.console.input("[dim]Press Enter to continue...[/dim]")

    def _prompt_positive_int(
        self,
        prompt: str,
        default: int,
        max_value: int | None = None,
    ) -> int | None:
        """Prompt user for a positive integer with validation.

        Args:
            prompt: Question to display
            default: Default value if user presses Enter
            max_value: Maximum allowed value (optional)

        Returns:
            Validated integer or None if user cancelled
        """
        while True:
            try:
                value_str = self.console.input(
                    f"[yellow]{prompt} [default: {default}][/yellow] "
                ).strip()

                if not value_str:
                    return default

                value = int(value_str)
                if value <= 0:
                    self.console.print("[red]✗ Must be a positive number[/red]")
                    continue
                if max_value and value > max_value:
                    self.console.print(f"[red]✗ Cannot exceed {max_value}[/red]")
                    continue
                return value

            except ValueError:
                self.console.print("[red]✗ Invalid number. Please try again.[/red]")
            except KeyboardInterrupt:
                return None

    def _run_pm7_batch(
        self,
        num_molecules: int,
        crest_timeout_minutes: int,
        mopac_timeout_minutes: int,
    ) -> None:
        """Execute PM7 batch processing with progress display.

        Args:
            num_molecules: Total molecules to process (also used as batch_size)
            crest_timeout_minutes: CREST timeout per molecule
            mopac_timeout_minutes: MOPAC timeout per molecule
        """
        from grimperium.crest_pm7.batch import (
            BatchCSVManager,
            BatchExecutionManager,
            BatchSortingStrategy,
            BatchStateManager,
            ConformerDetailManager,
            FixedTimeoutProcessor,
        )
        from grimperium.crest_pm7.config import PM7Config

        csv_path = DATA_DIR / "thermo_pm7.csv"
        detail_dir = DATA_DIR / "molecules_pm7" / "conformer_details"

        if not csv_path.exists():
            self.console.print(
                f"[red]✗ Batch CSV not found: {csv_path}[/red]\n"
                "[dim]Create a batch CSV with columns: mol_id, smiles, nheavy, status[/dim]"
            )
            self.console.input("[dim]Press Enter to continue...[/dim]")
            return

        try:
            self.console.print()
            self.console.print("[bold cyan]Initializing batch processor...[/bold cyan]")

            pm7_config = PM7Config()
            settings_manager: SettingsManager | None = getattr(
                self.controller, "settings_manager", None
            )
            if settings_manager is not None:
                settings_manager.apply_to_pm7_config(pm7_config)
            csv_manager = BatchCSVManager(csv_path)
            state_manager = BatchStateManager(
                csv_path.parent / "batch_state.csv", pm7_config
            )
            csv_manager.load_csv()

            detail_dir.mkdir(parents=True, exist_ok=True)
            detail_manager = ConformerDetailManager(detail_dir)

            processor = FixedTimeoutProcessor(
                config=pm7_config,
                crest_timeout_minutes=float(crest_timeout_minutes),
                mopac_timeout_minutes=float(mopac_timeout_minutes),
                enable_xtb_preopt=pm7_config.xtb_preopt,
                xtb_timeout_seconds=(
                    settings_manager.xtb.timeout_seconds
                    if settings_manager is not None
                    else None
                ),
            )

            exec_manager = BatchExecutionManager(
                csv_manager=csv_manager,
                state_manager=state_manager,
                detail_manager=detail_manager,
                pm7_config=pm7_config,
                processor_adapter=processor,
            )

            # Recover any molecules stuck in RUNNING from interrupted batches
            stuck_count = csv_manager.reset_stuck_running()
            if stuck_count > 0:
                self.console.print(
                    f"[yellow]Recovered {stuck_count} stuck RUNNING molecule(s) → PENDING[/yellow]"
                )

            batch_id = csv_manager.generate_batch_id()
            batch = csv_manager.select_batch(
                batch_id=batch_id,
                batch_size=num_molecules,
                crest_timeout_minutes=crest_timeout_minutes,
                mopac_timeout_minutes=mopac_timeout_minutes,
                strategy=BatchSortingStrategy.RERUN_FIRST_THEN_EASY,
            )

            if batch.is_empty:
                self.console.print(
                    "[yellow]No molecules available for processing.[/yellow]\n"
                    "[dim]All molecules may already be processed or none are PENDING.[/dim]"
                )
                self.console.input("[dim]Press Enter to continue...[/dim]")
                return

            self.console.print(
                f"[bold cyan]Starting batch {batch_id}: "
                f"{batch.size} molecules[/bold cyan]"
            )
            self.console.print()

            event_queue: Queue[ProgressEvent] = Queue()
            tracker = ProgressTracker(
                console=self.console,
                batch_size=batch.size,
            )
            csv_monitor = CSVMonitor(
                csv_path=csv_path,
                event_queue=event_queue,
                poll_interval_ms=500,
            )

            for mol in batch.molecules:
                tracker.register_molecule(mol.mol_id)
                csv_monitor.register_molecule(mol.mol_id)

            frame_idx = 0
            result = None
            csv_monitor.start()

            # Suppress noisy logs during Rich.Live rendering to avoid display corruption.
            previous_disable = logging.root.manager.disable
            logging.disable(logging.INFO)

            # Define batch_complete and batch_thread outside Live context for cleanup access
            batch_complete = threading.Event()
            batch_error: Exception | None = None
            batch_thread: threading.Thread | None = None

            def run_batch() -> None:
                nonlocal result, batch_error
                try:
                    result = exec_manager.execute_batch(batch)
                except Exception as e:
                    batch_error = e
                finally:
                    batch_complete.set()

            try:
                with Live(console=self.console, refresh_per_second=10) as live:
                    batch_thread = threading.Thread(target=run_batch, daemon=True)
                    batch_thread.start()

                    try:
                        while not batch_complete.is_set():
                            consume_events(event_queue, tracker)
                            display = self._render_pm7_batch_display(tracker, frame_idx)
                            live.update(display)
                            frame_idx += 1

                        if result is not None:
                            tracker.successful = result.success_count
                            tracker.failed = result.failed_count
                            tracker.skipped = result.skip_count
                            tracker.total_processed = result.total_count
                            display = self._render_pm7_batch_display(tracker, frame_idx)
                            live.update(display)

                        if batch_error is not None:
                            raise batch_error
                    except (KeyboardInterrupt, Exception):
                        batch_complete.set()
                        raise
            finally:
                logging.disable(previous_disable)
                # Ensure batch_thread completes and flushes CSV writes before cleanup
                if batch_thread is not None:
                    batch_complete.set()  # Ensure event is set even if exception occurred
                    batch_thread.join()
                csv_monitor.stop(timeout=1.0)

            self.console.print()
            self.console.print("[bold green]✓ Batch completed![/bold green]")
            self.console.print()
            self.console.print("[bold]Results Summary:[/bold]")
            if result is not None:
                self.console.print(f"  • Total processed: {result.total_count}")
                self.console.print(f"  • Successful: {result.success_count}")
                self.console.print(f"  • Failed: {result.failed_count}")
                self.console.print(f"  • Skipped: {result.skip_count}")

                if result.total_count > 0:
                    rate = result.success_count / result.total_count * 100
                    self.console.print(f"  • Success rate: {rate:.1f}%")
            else:
                self.console.print("  • No results available")

            self.console.print()
            self.console.print(
                f"[cyan]Results saved to:[/cyan] {PHASE_A_RESULTS_FILE.parent}"
            )

        except FileNotFoundError as e:
            self.console.print(f"[red]✗ File not found: {e}[/red]")
        except RuntimeError as e:
            self.console.print(f"[red]✗ Runtime error: {e}[/red]")
        except Exception as e:
            self.console.print(f"[red]✗ Unexpected error: {e}[/red]")

        self.console.print()
        self.console.input("[dim]Press Enter to continue...[/dim]")

    def _render_pm7_batch_display(
        self, tracker: ProgressTracker, frame_idx: int
    ) -> Panel:
        """Render PM7 batch display with single active progress and history.

        Args:
            tracker: ProgressTracker with current state
            frame_idx: Animation frame index for spinner

        Returns:
            Rich Panel with formatted display
        """
        header = tracker.render_batch_header()
        lines = [header, ""]

        current_mol = tracker.get_current_molecule_id()
        if current_mol:
            lines.append(tracker.render_molecule_line(current_mol, frame_idx))
            lines.append("")

        completed = tracker.get_completed_molecules()
        if completed:
            lines.append("[bold]Recent Completions:[/bold]")
            for mol_id, success, _elapsed_minutes in completed[-5:]:
                icon = ICONS["success"] if success else ICONS["error"]
                color = COLORS["success"] if success else COLORS["error"]
                label = "Completed successfully" if success else "Failed"
                lines.append(f"  [{color}]{icon} {mol_id} {label}[/{color}]")

        content = "\n".join(lines)

        return Panel(
            content,
            title=f"[bold {COLORS['batch']}]Batch Processing[/bold {COLORS['batch']}]",
            border_style=COLORS["batch"],
        )

    # ── Distributed mode (Phase 4) ────────────────────────────────────────────

    def _handle_run_calculation_menu(self) -> None:
        """Route between Local Mode and Distributed Mode."""
        options = [
            MenuOption(label="Local Mode", value="local"),
            MenuOption(label="Distributed Mode", value="distributed"),
        ]
        choice = show_back_menu(options=options, title="Run Calculation")
        if choice == "local":
            self.handle_calculate_pm7()
        elif choice == "distributed":
            self._handle_distributed_mode()

    def _handle_distributed_mode(self) -> None:
        """Orchestrate distributed mode as a state machine."""
        state: str | None = "check_port"
        while state is not None:
            if state == "check_port":
                state = self._state_check_port()
            elif state == "check_session":
                state = self._state_check_session()
            elif state == "config_menu":
                state = self._state_config_menu()
            elif state == "monitoring":
                state = self._state_monitoring()
            else:
                break

    @staticmethod
    def _run_server_process(csv_path_str: str, port: int) -> None:
        """Entry point for the server subprocess."""
        import uvicorn

        from grimperium.server.app import create_app
        from grimperium.server.config import ServerConfig

        config = ServerConfig(csv_path=csv_path_str)
        app = create_app(config)
        uvicorn.run(app, host="0.0.0.0", port=port, log_level="warning")  # noqa: S104

    def _start_server_in_background(
        self, csv_path: Path, port: int = 8000
    ) -> multiprocessing.Process:
        proc = multiprocessing.Process(
            target=DatabasesView._run_server_process,
            args=(str(csv_path), port),
            daemon=True,
        )
        proc.start()
        return proc

    def _fetch_worker_status(
        self, server_url: str, max_wait_s: float = 5.0
    ) -> list[dict[str, Any]]:
        """Fetch worker list from the server, retrying until ready."""
        import httpx

        deadline = time.monotonic() + max_wait_s
        while time.monotonic() < deadline:
            try:
                r = httpx.get(f"{server_url}/status", timeout=2.0)
                if r.status_code == 200:
                    return list(r.json().get("workers", []))
            except Exception as exc:  # noqa: BLE001
                logger.debug("Server not yet ready: %s", exc)
            time.sleep(0.5)
        return []

    def _render_workers_table(self, workers: list[dict[str, Any]]) -> Table:
        table = Table(
            title="Connected Workers",
            show_header=True,
            header_style=f"bold {COLORS['databases']}",
            border_style=COLORS["border"],
        )
        table.add_column("Worker ID")
        table.add_column("Hostname")
        table.add_column("Last Seen")
        table.add_column("Processing")

        if not workers:
            table.add_row("[dim]No workers connected[/dim]", "", "", "")
        else:
            for w in workers:
                table.add_row(
                    str(w.get("worker_id", "")),
                    str(w.get("hostname", "")),
                    str(w.get("last_seen") or "-"),
                    str(w.get("mol_id_current") or "-"),
                )
        return table

    def _screen_check_offline_workers(self, csv_path: Path) -> None:
        """Screen 3: if offline molecules exist, offer reassignment."""
        csv_manager = BatchCSVManager(csv_path)
        csv_manager.load_csv()
        df = csv_manager.df

        if df is None or "worker_status" not in df.columns:
            return

        offline_mask = (df["status"] == MoleculeStatus.ASSIGNED.value) & (
            df["worker_status"] == WorkerStatus.OFFLINE.value
        )
        offline_count = int(offline_mask.sum())

        if offline_count == 0:
            return

        self.console.print()
        if menu_confirm(
            f"⚠ {offline_count} molecules assigned to offline workers. Reassign?"
        ):
            n = csv_manager.reassign_offline_molecules()
            self.console.print(
                f"[green]✓ {n} molecules returned to pending pool[/green]"
            )
        self.wait_for_enter()

    # ── State machine helpers ─────────────────────────────────────────────────

    @staticmethod
    def _is_port_free(port: int = 8000) -> bool:
        import socket as _socket

        with _socket.socket(_socket.AF_INET, _socket.SOCK_STREAM) as s:
            return s.connect_ex(("localhost", port)) != 0

    @staticmethod
    def _server_is_responding(server_url: str) -> bool:
        import httpx

        try:
            r = httpx.get(f"{server_url}/status", timeout=2.0)
            return bool(r.status_code == 200)
        except Exception:  # noqa: BLE001
            return False

    def _state_check_port(self, port: int = 8000) -> str | None:
        server_url = f"http://localhost:{port}"
        while True:
            if self._is_port_free(port) or self._server_is_responding(server_url):
                return "check_session"
            options = [
                MenuOption("Verificar novamente", "retry"),
                MenuOption("Encerrar", "exit"),
            ]
            choice = show_back_menu(
                options=options,
                title=f"Porta {port} ocupada por outro processo",
            )
            if choice != "retry":
                return None

    def _state_check_session(self) -> str | None:
        session = load_session()
        if session is None:
            return "config_menu"

        if not self._server_is_responding(session.server_url):
            self.console.print(
                "[yellow]⚠ Sessão anterior encontrada mas servidor não responde. "
                "Iniciando nova sessão.[/yellow]"
            )
            delete_session()
            return "config_menu"

        options = [
            MenuOption("Entrar na sessão (processar localmente)", "join"),
            MenuOption("Acompanhar sessão (apenas monitorar)", "monitor"),
            MenuOption("Iniciar nova sessão", "new"),
            MenuOption("Encerrar sessão", "end"),
        ]
        choice = show_back_menu(options=options, title="Sessão ativa encontrada")

        if choice == "join":
            self._start_local_worker(session.server_url)
            return "monitoring"
        if choice == "monitor":
            return "monitoring"
        if choice == "new":
            self._shutdown_all_workers(session.server_url)
            delete_session()
            return "config_menu"
        if choice == "end":
            self._shutdown_all_workers(session.server_url)
            delete_session()
            return None
        return None

    def _shutdown_all_workers(self, server_url: str) -> None:
        import httpx

        try:
            httpx.post(f"{server_url}/shutdown/all", timeout=5.0)
        except Exception as exc:  # noqa: BLE001
            logger.debug("Shutdown request failed (server may be gone): %s", exc)

    def _start_local_worker(self, server_url: str) -> threading.Thread:
        client_cfg = WorkerClientConfig(server_url=server_url, worker_id="local")
        client = WorkerClient(client_cfg)

        server_cfg: dict[str, Any] = {}
        try:
            server_cfg = client.register()
        except Exception as exc:  # noqa: BLE001
            logger.debug("Local worker register failed, using defaults: %s", exc)

        config = WorkerConfig(
            server_url=server_url,
            worker_id="local",
            crest_timeout_minutes=int(server_cfg.get("crest_timeout_minutes", 60)),
            mopac_timeout_minutes=int(server_cfg.get("mopac_timeout_minutes", 30)),
            batch_size=int(server_cfg.get("batch_size", 10)),
        )
        runner = WorkerRunner(config=config, client=client)
        t = threading.Thread(
            target=runner.run,
            daemon=True,
            name="grimperium-local-worker",
        )
        t.start()
        return t

    def _configure_worker(
        self, server_url: str, worker_id: str, defaults: DistributedDefaults
    ) -> None:
        import httpx

        try:
            httpx.post(
                f"{server_url}/configure/{worker_id}",
                json={
                    "batch_size": defaults.batch_size,
                    "crest_timeout_minutes": defaults.crest_timeout_minutes,
                    "mopac_timeout_minutes": defaults.mopac_timeout_minutes,
                    "profile_name": defaults.profile_name,
                },
                timeout=5.0,
            )
        except Exception as exc:  # noqa: BLE001
            logger.debug("Configure worker %r failed: %s", worker_id, exc)

    def _fetch_workers_extended(self, server_url: str) -> list[dict[str, Any]]:
        import httpx

        try:
            r = httpx.get(f"{server_url}/workers/status", timeout=5.0)
            if r.status_code == 200:
                return list(r.json())
        except Exception as exc:  # noqa: BLE001
            logger.debug("GET /workers/status failed: %s", exc)
        return []

    def _render_monitoring_table(self, workers: list[dict[str, Any]]) -> Table:
        table = Table(
            title="Modo Distribuído — Monitoramento",
            show_header=True,
            header_style=f"bold {COLORS['databases']}",
            border_style=COLORS["border"],
        )
        table.add_column("Worker ID")
        table.add_column("Hostname")
        table.add_column("Proc", justify="right")
        table.add_column("OK", justify="right")
        table.add_column("Fail", justify="right")
        table.add_column("Skip", justify="right")
        table.add_column("Mol Ativo")

        if not workers:
            table.add_row("[dim]Nenhum worker[/dim]", "", "", "", "", "", "")
        else:
            for w in workers:
                table.add_row(
                    str(w.get("worker_id", "")),
                    str(w.get("hostname", "")),
                    str(w.get("processed", 0)),
                    str(w.get("successful", 0)),
                    str(w.get("failed", 0)),
                    str(w.get("skipped", 0)),
                    str(w.get("current_mol_id") or "—"),
                )
        return table

    def _state_config_menu(self) -> str | None:
        csv_path = DATA_DIR / "thermo_pm7.csv"
        server_port = 8000
        server_url = f"http://localhost:{server_port}"

        if not self._server_is_responding(server_url):
            if self._server_proc is None or not self._server_proc.is_alive():
                self._server_proc = self._start_server_in_background(
                    csv_path, server_port
                )
            deadline = time.monotonic() + 5.0
            while time.monotonic() < deadline:
                if self._server_is_responding(server_url):
                    break
                time.sleep(0.3)

        while True:
            workers = self._fetch_worker_status(server_url, max_wait_s=2.0)
            self.console.print()
            self.console.print(self._render_workers_table(workers))

            options = [
                MenuOption("Executar (iniciar processamento)", "run"),
                MenuOption("Atualizar lista de workers", "refresh"),
            ]
            choice = show_back_menu(
                options=options, title="Modo Distribuído — Configuração"
            )

            if choice == "refresh":
                continue
            if choice != "run":
                return None

            if not workers:
                self.console.print(
                    "[yellow]⚠ Nenhum worker conectado. Aguarde conexões e atualize.[/yellow]"
                )
                continue

            defaults = SettingsManager.load_distributed_defaults()
            for w in workers:
                worker_id = str(w.get("worker_id", ""))
                if worker_id:
                    self._configure_worker(server_url, worker_id, defaults)

            import httpx

            try:
                httpx.post(f"{server_url}/dispatch/start", timeout=5.0)
            except Exception as exc:  # noqa: BLE001
                logger.debug("POST /dispatch/start failed: %s", exc)

            worker_infos = [
                WorkerSessionInfo(
                    worker_id=str(w.get("worker_id", "")),
                    hostname=str(w.get("hostname", "")),
                    batch_size=defaults.batch_size,
                    profile_name=defaults.profile_name,
                    crest_timeout_minutes=defaults.crest_timeout_minutes,
                    mopac_timeout_minutes=defaults.mopac_timeout_minutes,
                )
                for w in workers
                if w.get("worker_id")
            ]
            from datetime import timezone

            session = SessionState(
                started_at=datetime.now(timezone.utc).isoformat(),
                server_url=server_url,
                workers=worker_infos,
            )
            save_session(session)

            self.console.print(
                f"[green]✓ Processamento iniciado com {len(workers)} worker(s)[/green]"
            )
            return "monitoring"

    def _state_monitoring(self) -> str | None:
        session = load_session()
        server_url = session.server_url if session else "http://localhost:8000"

        while True:
            workers = self._fetch_workers_extended(server_url)
            self.console.print()
            self.console.print(self._render_monitoring_table(workers))

            options = [
                MenuOption("Atualizar", "refresh"),
                MenuOption("Sair (servidor continua rodando)", "quit"),
            ]
            choice = show_back_menu(options=options, title="Monitoramento")
            if choice != "refresh":
                return None

    def run(self) -> str | None:
        """Run the databases view interaction loop."""
        while True:
            if self.selected_db:
                self.render_database_detail(self.selected_db)
                result = show_back_menu(
                    options=self.get_detail_menu_options(),
                    title="Actions",
                )
            else:
                self.render()
                result = show_back_menu(
                    options=self.get_menu_options(),
                    title="Select Database",
                )

            if result is None or result == "back":  # Ctrl+C
                if self.selected_db:
                    self.selected_db = None
                else:
                    return "main"
            else:
                next_view = self.handle_action(result)
                if next_view:
                    return next_view
