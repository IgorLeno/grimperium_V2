"""
Main application for GRIMPERIUM CLI.

This module provides the main entry point and application loop.
"""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

from rich.panel import Panel

from grimperium.cli.controller import CliController
from grimperium.cli.menu import show_main_menu
from grimperium.cli.styles import (
    BANNER,
    COLORS,
    FOOTER_NAVIGATION,
    SUBTITLE,
    VERSION_LINE,
)
from grimperium.cli.views import (
    AboutView,
    BatchView,
    CalcView,
    DatabasesView,
    MethodsView,
    ModelsView,
    ResultsView,
    SettingsView,
)


class GrimperiumCLI:
    """
    Main CLI application class.

    Manages the application lifecycle, view registration, and main loop.
    """

    def __init__(self) -> None:
        """Initialize the CLI application."""
        self.controller = CliController()
        self.console = self.controller.console
        self._register_views()

    def _register_views(self) -> None:
        """Register all views with the controller."""
        self.controller.register_view("batch", BatchView(self.controller))
        self.controller.register_view("calc", CalcView(self.controller))
        self.controller.register_view("methods", MethodsView(self.controller))
        self.controller.register_view("databases", DatabasesView(self.controller))
        self.controller.register_view("models", ModelsView(self.controller))
        self.controller.register_view("results", ResultsView(self.controller))
        self.controller.register_view("settings", SettingsView(self.controller))
        self.controller.register_view("about", AboutView(self.controller))

    def show_welcome(self) -> None:
        """Display the welcome screen with ASCII banner."""
        self.console.clear()
        self.console.print()

        # Banner panel
        banner_content = f"{BANNER}\n{SUBTITLE}\n{VERSION_LINE}"
        self.console.print(
            Panel(
                banner_content,
                border_style=COLORS["calc"],
                padding=(1, 4),
            )
        )
        self.console.print()

    def display_main_menu(self) -> str | None:
        """
        Display and handle the main menu.

        Returns:
            Selected menu option or None if cancelled
        """
        self.console.print(
            Panel(
                f"[bold]MAIN MENU[/bold]\n\n{FOOTER_NAVIGATION}",
                border_style=COLORS["border"],
                padding=(0, 2),
            )
        )
        self.console.print()

        summary = self.controller.session_summary()
        return show_main_menu(
            property_label=summary["property"],
            method_label=summary["method"],
            dataset_label=summary["dataset"],
            model_label=summary["model"],
            status=summary["status"],
        )

    def run_view(self, view_name: str) -> str | None:
        """
        Run a specific view.

        Args:
            view_name: Name of the view to run

        Returns:
            Next view to navigate to, or None to stay at main menu
        """
        view = self.controller.get_view(view_name)
        if view is None:
            self.console.print(f"[error]View '{view_name}' not found[/error]")
            return None

        try:
            return view.run()
        except KeyboardInterrupt:
            return "main"

    def run(self, skip_preflight: bool = False) -> int:
        """
        Run the main application loop.

        Args:
            skip_preflight: If True, skip the environment check entirely.

        Returns:
            Exit code (0 for success, non-zero for errors)
        """
        self.controller.start()

        # -- Preflight environment check ------------------------------------
        if not skip_preflight:
            from grimperium.cli.preflight import run_preflight

            self.show_welcome()
            force = "--force-preflight" in sys.argv
            if not run_preflight(self.console, force=force):
                self.console.print(
                    "[muted]Exiting due to failed preflight checks.[/muted]"
                )
                self.controller.stop()
                return 2
        # ---------------------------------------------------------------------

        try:
            while self.controller.is_running():
                # Show welcome screen and main menu
                self.show_welcome()
                selection = self.display_main_menu()

                if selection is None:
                    # User pressed Ctrl+C at main menu - just break, stop() in finally
                    break

                # Navigate to selected view
                if selection == "main":
                    continue

                # Run the selected view
                next_view = self.run_view(selection)

                # Handle navigation from view
                while next_view and next_view != "main":
                    next_view = self.run_view(next_view)

        except KeyboardInterrupt:
            self.console.print()
            self.console.print("[muted]Interrupted by user[/muted]")

        except Exception as e:
            self.console.print()
            self.console.print(f"[error]Unexpected error: {e}[/error]")
            return 1

        finally:
            self.controller.stop()

        return 0


def main() -> int:
    """
    Main entry point for the GRIMPERIUM CLI.

    Returns:
        Exit code (0 for success, non-zero for errors)
    """
    from grimperium.cli.constants import get_project_root
    from grimperium.utils.logging import setup_logging

    log_file = get_project_root() / "logs" / "grimperium.log"
    setup_logging(
        level="DEBUG",
        log_file=log_file,
        console_level="WARNING",
    )

    # Non-interactive subcommand — runs and exits without loading the TUI
    if len(sys.argv) > 1 and sys.argv[1] == "analyze-errors":
        return _run_analyze_errors(sys.argv[2:])

    skip_preflight = "--skip-preflight" in sys.argv
    app = GrimperiumCLI()
    return app.run(skip_preflight=skip_preflight)


def _run_analyze_errors(args: list[str]) -> int:
    """Execute the analyze-errors subcommand and return an exit code."""
    import argparse
    from pathlib import Path

    import pandas as pd

    parser = argparse.ArgumentParser(
        prog="grimperium analyze-errors",
        description="Analyze prediction errors and detect outliers.",
    )
    parser.add_argument("--input", required=True, type=Path, metavar="PATH")
    parser.add_argument("--show-outliers", action="store_true")
    parser.add_argument("--top", type=int, default=20, metavar="N")
    parser.add_argument("--export-outliers", type=Path, metavar="PATH")
    parser.add_argument("--html-report", type=Path, metavar="PATH")
    parser.add_argument("--threshold", type=float, default=10.0)
    parser.add_argument("--zscore", type=float, default=3.0)
    parser.add_argument("--percentile", type=float, default=95.0)

    parsed = parser.parse_args(args)

    from grimperium.ml.error_analysis import ErrorAnalysisConfig, analyze

    input_path: Path = parsed.input
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        return 1

    try:
        df = pd.read_csv(input_path)
    except Exception as exc:
        print(f"Error reading CSV: {exc}", file=sys.stderr)
        return 1

    config = ErrorAnalysisConfig(
        threshold_kcalmol=parsed.threshold,
        zscore_threshold=parsed.zscore,
        percentile_outlier=parsed.percentile,
    )

    try:
        result = analyze(df, config)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    s = result.summary

    # Mandatory summary (always printed)
    print("Error analysis complete")
    print(f"Input: {input_path}")
    print(f"Valid molecules analyzed: {int(s['n_molecules']):,}")
    print(
        f"MAE: {float(s['mae']):.2f} kcal/mol | "
        f"RMSE: {float(s['rmse']):.2f} kcal/mol | "
        f"R²: {float(s['r2']):.4f}"
    )
    print(f"Max absolute error: {float(s['max_error']):.2f} kcal/mol")
    print()
    print(f"Outliers detected: {int(s['n_outliers'])}")
    print(
        f"Critical molecules: {int(s.get('n_critical', 0))} | "
        f"Extreme molecules: {int(s.get('n_extreme', 0))}"
    )
    print()
    print("Top anomalous molecule IDs:")
    for _, row in result.top_n_df.head(parsed.top).iterrows():
        mol_id = str(row.get("mol_id", "unknown"))
        abs_err = float(row["abs_error"])
        severity = str(row.get("severity", "-"))
        zscore = float(row.get("abs_error_zscore", 0.0))
        print(
            f"  {int(row['rank'])}. {mol_id} | "
            f"abs_error={abs_err:.2f} | "
            f"severity={severity} | "
            f"zscore={zscore:.2f}"
        )

    # Optional: outlier detail table
    if parsed.show_outliers:
        _cli_print_outlier_table(result.outliers_df)

    # Optional: export outliers CSV
    if parsed.export_outliers:
        out: Path = parsed.export_outliers
        out.parent.mkdir(parents=True, exist_ok=True)
        result.outliers_df.to_csv(out, index=False)
        print(f"\nOutliers exported to: {out}")

    # Optional: HTML report
    if parsed.html_report:
        from grimperium.ml.html_report import generate_html_report

        html_path: Path = parsed.html_report
        generate_html_report(result, html_path)
        print(f"\nHTML report saved to: {html_path}")

    return 0


def _cli_print_outlier_table(outliers_df: pd.DataFrame) -> None:
    if outliers_df.empty:
        print("\nNo outliers detected.")
        return

    cols = ["mol_id", "abs_error", "severity", "outlier_score"]
    present = [c for c in cols if c in outliers_df.columns]
    print(f"\nOutlier Details ({len(outliers_df)} molecules):")
    print("-" * 70)
    print("  ".join(c.ljust(20) for c in present))
    print("-" * 70)
    for _, row in outliers_df.iterrows():
        parts = []
        for c in present:
            val = row[c]
            if isinstance(val, float):
                parts.append(f"{val:.4f}".ljust(20))
            else:
                parts.append(str(val).ljust(20))
        print("  ".join(parts))


if __name__ == "__main__":
    sys.exit(main())
