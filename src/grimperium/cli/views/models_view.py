"""
Models view for GRIMPERIUM CLI.

Displays and manages trained ML models using real persistence data.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING, Any

from rich.panel import Panel
from rich.table import Table

from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView
from grimperium.ml.persistence import load_model_metadata, save_model

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController

_MODEL_NAME = "DeltaLearner v1"
_MODEL_ALGORITHM = "KRR + XGBoost Ensemble"
_MODEL_PATH = Path(
    os.environ.get("GRIMPERIUM_MODEL_PATH", "models/delta_learner_v1.joblib")
)
_DATA_PATH = Path(os.environ.get("GRIMPERIUM_DATA_PATH", "data/thermo_pm7.csv"))


def _safe_metric(value: Any, default: float = 0.0) -> float:
    """Return *value* if numeric, otherwise *default*."""
    return default if value is None else float(value)


class ModelsView(BaseView):
    """View for managing ML models."""

    name = "models"
    title = "Models"
    icon = ICONS["models"]
    color = COLORS["models"]

    def __init__(self, controller: CliController) -> None:
        """Initialize the models view."""
        super().__init__(controller)
        self.selected_model: bool = False

    def _load_model_info(self) -> dict[str, Any] | None:
        """Try to load model metadata. Returns None if model not trained yet."""
        try:
            return load_model_metadata(_MODEL_PATH)
        except FileNotFoundError:
            return None

    def render(self) -> None:
        """Render the models overview."""
        self.clear_screen()
        self.show_header()

        metadata = self._load_model_info()

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

        if metadata is not None:
            test_metrics = metadata.get("metrics", {}).get("test", {})
            mae_val = test_metrics.get("mae")
            r2_val = test_metrics.get("r2")
            mae = f"{mae_val:.3f}" if mae_val is not None else "-"
            r2 = f"{r2_val:.4f}" if r2_val is not None else "-"
            status = (
                f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
            )
        else:
            mae = "-"
            r2 = "-"
            status = f"[{COLORS['in_dev']}]{ICONS['in_dev']} Not Trained[/{COLORS['in_dev']}]"

        # Highlight default model
        name = _MODEL_NAME
        if self.controller.current_model == _MODEL_NAME:
            name = f"[bold]{_MODEL_NAME}[/bold] ★"

        table.add_row(name, _MODEL_ALGORITHM, mae, r2, status)

        self.console.print(table)
        self.console.print()
        self.console.print(
            f"[{COLORS['muted']}]★ = Default model for predictions[/{COLORS['muted']}]"
        )
        self.console.print()

    def render_model_detail(self, metadata: dict[str, Any] | None) -> None:
        """Render detailed view for a specific model."""
        self.clear_screen()
        self.show_header()

        if metadata is None:
            info = f"""
[bold]Name:[/bold]          {_MODEL_NAME}
[bold]Algorithm:[/bold]     {_MODEL_ALGORITHM}
[bold]Status:[/bold]        [{COLORS['in_dev']}]{ICONS['in_dev']} Not Trained[/{COLORS['in_dev']}]

Model has not been trained yet. Use "Train New Model" to train.
"""
        else:
            test_metrics = metadata.get("metrics", {}).get("test", {})
            train_metrics = metadata.get("metrics", {}).get("train", {})

            mae_str = f"{_safe_metric(test_metrics.get('mae')):.3f}"
            r2_str = f"{_safe_metric(test_metrics.get('r2')):.4f}"
            rmse_str = f"{_safe_metric(test_metrics.get('rmse')):.3f}"
            mape_str = f"{_safe_metric(test_metrics.get('mape')):.2f}"
            max_err_str = f"{_safe_metric(test_metrics.get('max_error')):.3f}"

            train_rmse = f"{_safe_metric(train_metrics.get('rmse')):.3f}"
            train_mae = f"{_safe_metric(train_metrics.get('mae')):.3f}"
            train_r2 = f"{_safe_metric(train_metrics.get('r2')):.4f}"

            info = f"""
[bold]Name:[/bold]          {_MODEL_NAME}
[bold]Algorithm:[/bold]     {_MODEL_ALGORITHM}
[bold]Version:[/bold]       {metadata.get('version', 'unknown')}
[bold]Trained at:[/bold]    {metadata.get('trained_at', 'unknown')}
[bold]Status:[/bold]        [{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]

[bold]Test Metrics:[/bold]
  RMSE:            {rmse_str} kcal/mol
  MAE:             {mae_str} kcal/mol
  R² Score:        {r2_str}
  MAPE:            {mape_str}%
  Max Error:       {max_err_str} kcal/mol

[bold]Train Metrics:[/bold]
  RMSE:            {train_rmse} kcal/mol
  MAE:             {train_mae} kcal/mol
  R² Score:        {train_r2}
"""

        self.console.print(
            Panel(
                info,
                title=f"[bold {COLORS['models']}]{_MODEL_NAME}[/bold {COLORS['models']}]",
                border_style=COLORS["models"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the models view."""
        options = [
            MenuOption(
                label=f"View {_MODEL_NAME}",
                value="view_model",
                icon=ICONS["models"],
            ),
            MenuOption(
                label="Train New Model",
                value="train",
            ),
            MenuOption(
                label="Compare Models",
                value="compare",
                disabled=True,
                disabled_reason="In Development",
            ),
        ]
        return options

    def get_detail_menu_options(self) -> list[MenuOption]:
        """Return menu options for model detail view."""
        options = []

        is_default = self.controller.current_model == _MODEL_NAME
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

    def handle_action(self, action: str | None) -> str | None:
        """Handle menu actions."""
        # Handle None or "back" action
        if action is None or action == "back":
            if self.selected_model:
                self.selected_model = False
                return "models"  # Stay in models view, just return to list
            return "main"

        if action == "view_model":
            self.selected_model = True
            return None

        if action == "set_default":
            self.controller.set_model(_MODEL_NAME)
            self.show_success(f"Default model set to {_MODEL_NAME}")
            return None

        if action == "train":
            self._handle_train()
            return None

        if action == "test":
            self._handle_test()
            return None

        # Handle in-development features
        if action in ["compare", "retrain", "export"]:
            self.show_in_development(action.title())
            return None

        return None

    def _handle_train(self) -> None:
        """Train a new model and save it."""
        from grimperium.ml.trainer import train as train_model

        self.console.print()
        self.console.print(
            f"[bold {COLORS['models']}]Training {_MODEL_NAME}...[/bold {COLORS['models']}] "
            "(this may take a minute)"
        )
        self.console.print()

        try:
            result = train_model(
                _DATA_PATH,
                return_pipeline=True,
                random_state=42,
            )
            if not isinstance(result, tuple) or len(result) != 4:
                msg = (
                    f"train_model(return_pipeline=True) returned "
                    f"{type(result).__name__} with "
                    f"{len(result) if isinstance(result, tuple) else 'N/A'} "
                    f"elements; expected a 4-tuple"
                )
                raise TypeError(msg)
            learner, train_m, test_m, pipeline = result

            bundle: dict[str, Any] = {
                "learner": learner,
                "pipeline": pipeline,
                "metrics": {
                    "train": train_m,
                    "test": test_m,
                },
            }
            save_model(bundle, _MODEL_PATH)

            gate_pass = test_m.get("gate_pass", False)
            gate_icon = ICONS["success"] if gate_pass else ICONS["error"]

            result_text = f"""
[bold]Training Complete![/bold]

[bold]Test Results:[/bold]
  RMSE:        {_safe_metric(test_m.get('rmse')):.3f} kcal/mol
  MAE:         {_safe_metric(test_m.get('mae')):.3f} kcal/mol
  R² Score:    {_safe_metric(test_m.get('r2')):.4f}
  Gate Pass:   {gate_icon} {'Yes' if gate_pass else 'No'}
"""

            self.console.print(
                Panel(
                    result_text,
                    title=f"[bold {COLORS['success']}]Training Results[/bold {COLORS['success']}]",
                    border_style=COLORS["success"],
                    padding=(1, 2),
                )
            )
            self.show_success(f"Model trained and saved to {_MODEL_PATH}")
        except Exception as e:
            self.show_error(f"Training failed: {e}")

        self.wait_for_enter()

    def _handle_test(self) -> None:
        """Display test metrics for the trained model."""
        metadata = self._load_model_info()
        if metadata is None:
            self.show_error("No trained model found. Train first.")
            self.wait_for_enter()
            return

        test_metrics = metadata.get("metrics", {}).get("test", {})
        gate_pass = test_metrics.get("gate_pass", False)
        gate_icon = ICONS["success"] if gate_pass else ICONS["error"]

        result_text = f"""
[bold]Model:[/bold]       {_MODEL_NAME}
[bold]Version:[/bold]     {metadata.get('version', 'unknown')}
[bold]Trained at:[/bold]  {metadata.get('trained_at', 'unknown')}

[bold]Test Metrics:[/bold]
  RMSE:        {_safe_metric(test_metrics.get('rmse')):.3f} kcal/mol
  MAE:         {_safe_metric(test_metrics.get('mae')):.3f} kcal/mol
  R² Score:    {_safe_metric(test_metrics.get('r2')):.4f}
  MAPE:        {_safe_metric(test_metrics.get('mape')):.2f}%
  Max Error:   {_safe_metric(test_metrics.get('max_error')):.3f} kcal/mol
  Gate Pass:   {gate_icon} {'Yes' if gate_pass else 'No'}
"""

        self.console.print(
            Panel(
                result_text,
                title=f"[bold {COLORS['models']}]Model Test Results[/bold {COLORS['models']}]",
                border_style=COLORS["models"],
                padding=(1, 2),
            )
        )
        self.wait_for_enter()

    def run(self) -> str | None:
        """Run the models view interaction loop."""
        while True:
            if self.selected_model:
                metadata = self._load_model_info()
                self.render_model_detail(metadata)
                result = show_back_menu(
                    options=self.get_detail_menu_options(),
                    title="Actions",
                )
            else:
                self.render()
                result = show_back_menu(
                    options=self.get_menu_options(),
                    title="Select Model",
                )

            # Always delegate to handle_action
            next_view = self.handle_action(result)
            if next_view:
                return next_view
