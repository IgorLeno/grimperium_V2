"""
Models view for GRIMPERIUM CLI.

Displays and manages trained ML models.
"""

from typing import TYPE_CHECKING, Optional

from rich.panel import Panel
from rich.table import Table

from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.mock_data import MODELS, Model
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController


class ModelsView(BaseView):
    """View for managing ML models."""

    name = "models"
    title = "Models"
    icon = ICONS["models"]
    color = COLORS["models"]

    def __init__(self, controller: "CliController") -> None:
        """Initialize the models view."""
        super().__init__(controller)
        self.selected_model: Optional[Model] = None

    def render(self) -> None:
        """Render the models overview."""
        self.clear_screen()
        self.show_header()

        # Models table
        table = Table(
            title="Trained Models",
            show_header=True,
            header_style=f"bold {COLORS['models']}",
            border_style=COLORS["border"],
        )
        table.add_column("Name", style="bold")
        table.add_column("Algorithm")
        table.add_column("MAE", justify="right")
        table.add_column("R²", justify="right")
        table.add_column("Status")

        for model in MODELS:
            if model.status == "ready":
                status = f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
                mae = f"{model.mae:.2f}" if model.mae else "-"
                r2 = f"{model.r2:.3f}" if model.r2 else "-"
            else:
                status = (
                    f"[{COLORS['in_dev']}]{ICONS['in_dev']} In Dev[/{COLORS['in_dev']}]"
                )
                mae = "-"
                r2 = "-"

            # Highlight default model
            name = model.name
            if model.name == self.controller.current_model:
                name = f"[bold]{model.name}[/bold] ★"

            table.add_row(
                name,
                model.algorithm,
                mae,
                r2,
                status,
            )

        self.console.print(table)
        self.console.print()
        self.console.print(
            f"[{COLORS['muted']}]★ = Default model for predictions[/{COLORS['muted']}]"
        )
        self.console.print()

    def render_model_detail(self, model: Model) -> None:
        """Render detailed view for a specific model."""
        self.clear_screen()
        self.show_header()

        # Model info panel
        status_text = (
            f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
            if model.status == "ready"
            else f"[{COLORS['in_dev']}]{ICONS['in_dev']} In Development[/{COLORS['in_dev']}]"
        )

        info = f"""
[bold]Name:[/bold]          {model.name}
[bold]Algorithm:[/bold]     {model.algorithm}
[bold]Status:[/bold]        {status_text}
"""
        if model.status == "ready":
            info += f"""
[bold]Performance Metrics:[/bold]
  MAE:             {model.mae:.2f} kcal/mol
  R² Score:        {model.r2:.4f}
  Training Date:   {model.training_date}
  File Size:       {model.file_size}

[bold]Hyperparameters:[/bold]
"""
            for key, value in model.hyperparameters.items():
                info += f"  {key}: {value}\n"
        else:
            info += """
[bold]Hyperparameters:[/bold] (planned)
"""
            for key, value in model.hyperparameters.items():
                info += f"  {key}: {value}\n"

        self.console.print(
            Panel(
                info,
                title=f"[bold {COLORS['models']}]{model.name}[/bold {COLORS['models']}]",
                border_style=COLORS["models"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the models view."""
        options = []
        for model in MODELS:
            disabled = model.status != "ready"
            options.append(
                MenuOption(
                    label=f"View {model.name}",
                    value=f"view_{model.name}",
                    icon=ICONS["models"],
                    disabled=disabled,
                    disabled_reason="In Development" if disabled else "",
                )
            )

        # Additional options (in development)
        options.extend(
            [
                MenuOption(
                    label="Train New Model",
                    value="train",
                    disabled=True,
                    disabled_reason="In Development",
                ),
                MenuOption(
                    label="Compare Models",
                    value="compare",
                    disabled=True,
                    disabled_reason="In Development",
                ),
            ]
        )

        return options

    def get_detail_menu_options(self, model: Model) -> list[MenuOption]:
        """Return menu options for model detail view."""
        options = []

        if model.status == "ready":
            is_default = model.name == self.controller.current_model
            options.append(
                MenuOption(
                    label="Set as Default" if not is_default else "Already Default",
                    value="set_default",
                    disabled=is_default,
                    disabled_reason="Current default" if is_default else "",
                )
            )

        options.extend(
            [
                MenuOption(
                    label="Test Model",
                    value="test",
                    disabled=True,
                    disabled_reason="In Development",
                ),
                MenuOption(
                    label="Retrain Model",
                    value="retrain",
                    disabled=True,
                    disabled_reason="In Development",
                ),
                MenuOption(
                    label="Export Model",
                    value="export",
                    disabled=True,
                    disabled_reason="In Development",
                ),
            ]
        )

        return options

    def handle_action(self, action: str) -> Optional[str]:
        """Handle menu actions."""
        if action == "back":
            if self.selected_model:
                self.selected_model = None
                return None
            return "main"

        if action.startswith("view_"):
            model_name = action.replace("view_", "")
            for model in MODELS:
                if model.name == model_name:
                    self.selected_model = model
                    return None

        if action == "set_default" and self.selected_model:
            self.controller.set_model(self.selected_model.name)
            self.show_success(f"Default model set to {self.selected_model.name}")
            return None

        # Handle in-development features
        if action in ["train", "compare", "test", "retrain", "export"]:
            self.show_in_development(action.title())
            return None

        return None

    def run(self) -> Optional[str]:
        """Run the models view interaction loop."""
        while True:
            if self.selected_model:
                self.render_model_detail(self.selected_model)
                result = show_back_menu(
                    options=self.get_detail_menu_options(self.selected_model),
                    title="Actions",
                )
            else:
                self.render()
                result = show_back_menu(
                    options=self.get_menu_options(),
                    title="Select Model",
                )

            if result is None or result == "back":  # Ctrl+C
                if self.selected_model:
                    self.selected_model = None
                else:
                    return "main"
            else:
                next_view = self.handle_action(result)
                if next_view:
                    return next_view
