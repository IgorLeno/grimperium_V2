"""
Results view for GRIMPERIUM CLI.

Displays performance analytics and divergence analysis.
"""

from typing import Optional

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

        # Sort models by MAE (best first)
        ready_models = [m for m in MODELS if m.status == "ready" and m.mae is not None]
        sorted_models = sorted(ready_models, key=lambda x: x.mae or float("inf"))

        for rank, model in enumerate(sorted_models, 1):
            # Highlight best model
            if rank == 1:
                name = (
                    f"[bold {COLORS['success']}]{model.name}[/bold {COLORS['success']}]"
                )
                mae = f"[bold {COLORS['success']}]{model.mae:.2f}[/bold {COLORS['success']}]"
                r2 = f"[bold {COLORS['success']}]{model.r2:.4f}[/bold {COLORS['success']}]"
                rank_str = f"[bold {COLORS['success']}]🥇 1[/bold {COLORS['success']}]"
            elif rank == 2:
                name = model.name
                mae = f"{model.mae:.2f}"
                r2 = f"{model.r2:.4f}"
                rank_str = "🥈 2"
            elif rank == 3:
                name = model.name
                mae = f"{model.mae:.2f}"
                r2 = f"{model.r2:.4f}"
                rank_str = "🥉 3"
            else:
                name = model.name
                mae = f"{model.mae:.2f}"
                r2 = f"{model.r2:.4f}"
                rank_str = str(rank)

            table.add_row(name, model.algorithm, mae, r2, rank_str)

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
            "HIGH": "#FF8800",
            "CRITICAL": COLORS["error"],
        }

        total_molecules = sum(d.count for d in DIVERGENCE_STATS)

        for stat in DIVERGENCE_STATS:
            color = severity_colors.get(stat.severity, COLORS["muted"])

            # Create visual bar
            bar_length = int(stat.percentage / 5)  # Scale to max ~20 chars
            bar = f"[{color}]{'█' * bar_length}{'░' * (20 - bar_length)}[/{color}]"

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
• [#FF8800]HIGH (25-50%)[/#FF8800]: Significant corrections needed, challenging cases
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

    def handle_action(self, action: str) -> Optional[str]:
        """Handle menu actions."""
        if action == "back":
            return "main"

        # Handle in-development features
        if action in ["export", "detailed", "charts"]:
            self.show_in_development(action.title())
            return None

        return None

    def run(self) -> Optional[str]:
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
