"""
Results view for GRIMPERIUM CLI.

Displays performance analytics and divergence analysis.
"""

from rich.markup import escape
from rich.panel import Panel
from rich.table import Table

from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.mock_data import DIVERGENCE_STATS, MODELS
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView


class ResultsView(BaseView):
    """View for performance analytics and results."""

    name = "results"
    title = "Results"
    icon = ICONS["results"]
    color = COLORS["results"]

    def render(self) -> None:
        """Render the results overview."""
        self.clear_screen()
        self.show_header()

        self._render_model_comparison()
        self._render_divergence_analysis()

    def _render_model_comparison(self) -> None:
        """Render model performance comparison table."""
        table = Table(
            title="Model Performance Comparison",
            show_header=True,
            header_style=f"bold {COLORS['results']}",
            border_style=COLORS["border"],
        )
        table.add_column("Model", style="bold")
        table.add_column("Algorithm")
        table.add_column("MAE (kcal/mol)", justify="right")
        table.add_column("R²", justify="right")
        table.add_column("Rank", justify="center")

        # Sort models by MAE (best first), filter for ready models with mae and r2
        ready_models = [
            m
            for m in MODELS
            if m.status == "ready" and m.mae is not None and m.r2 is not None
        ]
        # Type narrowing: mae is guaranteed non-None by filter above
        sorted_models = sorted(ready_models, key=lambda x: x.mae or 0.0)

        for rank, model in enumerate(sorted_models, 1):
            # Escape user-provided strings to prevent Rich markup injection
            name = escape(model.name)
            algo = escape(model.algorithm)
            mae_str = f"{model.mae:.2f}"
            r2_str = f"{model.r2:.4f}" if model.r2 is not None else "-"

            # Choose rank string based on position
            if rank == 1:
                rank_str = "🥇"
            elif rank == 2:
                rank_str = "🥈"
            elif rank == 3:
                rank_str = "🥉"
            else:
                rank_str = str(rank)

            # Apply highlighting only for rank 1
            if rank == 1:
                name = f"[bold {COLORS['success']}]{name}[/bold {COLORS['success']}]"
                mae_str = (
                    f"[bold {COLORS['success']}]{mae_str}[/bold {COLORS['success']}]"
                )
                r2_str = (
                    f"[bold {COLORS['success']}]{r2_str}[/bold {COLORS['success']}]"
                )
                rank_str = (
                    f"[bold {COLORS['success']}]{rank_str} 1[/bold {COLORS['success']}]"
                )
            elif rank == 2:
                rank_str = f"{rank_str} 2"
            elif rank == 3:
                rank_str = f"{rank_str} 3"

            table.add_row(name, algo, mae_str, r2_str, rank_str)

        self.console.print(table)
        self.console.print()

    def _render_divergence_analysis(self) -> None:
        """Render CBS vs PM7 divergence analysis."""
        # Divergence distribution table
        table = Table(
            title="CBS vs PM7 Divergence Distribution",
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

        total_molecules = sum(d.count for d in DIVERGENCE_STATS)

        for stat in DIVERGENCE_STATS:
            color = severity_colors.get(stat.severity, COLORS["muted"])

            # Clamp percentage to [0, 100] before creating bar
            pct = max(0, min(stat.percentage, 100))
            bar_length = int(pct / 5)  # Scale to max ~20 chars
            bar_length = max(0, min(bar_length, 20))
            fill_count = bar_length
            empty_count = 20 - fill_count
            bar = f"[{color}]{'█' * fill_count}{'░' * empty_count}[/{color}]"

            table.add_row(
                f"[{color}]{stat.severity}[/{color}]",
                f"{stat.range_min:.0f}% - {stat.range_max:.0f}%",
                f"{stat.count:,}",
                f"{stat.percentage:.1f}%",
                bar,
            )

        self.console.print(table)
        self.console.print()

        # Summary panel
        low_medium = sum(
            d.count for d in DIVERGENCE_STATS if d.severity in ["LOW", "MEDIUM"]
        )
        low_medium_pct = (
            (low_medium / total_molecules) * 100 if total_molecules > 0 else 0
        )

        summary = f"""
[bold]Key Findings:[/bold]

• Total molecules analyzed: {total_molecules:,}
• Molecules with LOW/MEDIUM divergence: {low_medium:,} ({low_medium_pct:.1f}%)
• The Delta-Learning approach effectively corrects PM7 predictions

[bold]Interpretation:[/bold]

• [{COLORS['success']}]LOW (0-10%)[/{COLORS['success']}]: PM7 predictions are accurate, small corrections needed
• [{COLORS['warning']}]MEDIUM (10-25%)[/{COLORS['warning']}]: Moderate corrections, ML performs well
• [{COLORS['high']}]HIGH (25-50%)[/{COLORS['high']}]: Significant corrections needed, challenging cases
• [{COLORS['error']}]CRITICAL (>50%)[/{COLORS['error']}]: Large deviations, may require special handling
"""
        self.console.print(
            Panel(
                summary,
                title=f"[bold {COLORS['results']}]Analysis Summary[/bold {COLORS['results']}]",
                border_style=COLORS["border"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the results view."""
        return [
            MenuOption(
                label="Export Report",
                value="export",
                disabled=True,
                disabled_reason="In Development",
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

        # Handle in-development features
        if action in ["export", "detailed", "charts"]:
            self.show_in_development(action.title())
            return None

        return None

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
