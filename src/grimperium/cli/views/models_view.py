"""
Models view for GRIMPERIUM CLI.

Displays and manages trained ML models using real persistence data.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING, Any

from rich.panel import Panel
from rich.table import Table

from grimperium.cli.menu import MenuOption, show_back_menu
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView
from grimperium.ml.gate import evaluate_gate
from grimperium.ml.persistence import load_model_metadata, save_model
from grimperium.ml.trainer import train

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController

logger = logging.getLogger(__name__)

_MODEL_ALGORITHM = "KRR + XGBoost Ensemble"


def discover_available_models(
    models_dir: Path | None = None,
    active_model_path: Path | None = None,
) -> list[dict[str, Any]]:
    """Return metadata dicts for all valid ``.joblib`` models in *models_dir*.

    Args:
        models_dir: Directory to scan; defaults to ``models/`` relative to CWD.
        active_model_path: Resolved path of the currently active model,
            used to populate the ``is_active`` flag.

    Returns:
        List of dicts with keys: ``path``, ``display_name``, ``algorithm``,
        ``mae``, ``r2``, ``trained_at``, ``version``, ``is_active``.
        Returns an empty list when the directory does not exist.
    """
    search_dir = models_dir or Path("models/")
    if not search_dir.is_dir():
        return []
    results: list[dict[str, Any]] = []
    for p in sorted(search_dir.glob("*.joblib")):
        p = p.resolve()
        try:
            meta = load_model_metadata(p)
        except Exception:
            logger.warning("Skipping corrupted or unreadable model: %s", p)
            continue
        test_metrics = meta.get("metrics", {}).get("test", {})
        results.append(
            {
                "path": p,
                "display_name": p.stem.replace("_", " ").title(),
                "algorithm": _MODEL_ALGORITHM,
                "mae": test_metrics.get("mae"),
                "r2": test_metrics.get("r2"),
                "trained_at": meta.get("trained_at", "unknown"),
                "version": meta.get("version", "unknown"),
                "is_active": p == active_model_path,
            }
        )
    return results


def _safe_metric(value: Any, default: float = 0.0) -> float:
    """Return *value* if numeric, otherwise *default*."""
    return default if value is None else float(value)


def _format_split_summary(meta: dict[str, Any]) -> str:
    """Render training-volume metadata when present in the bundle."""
    n_total = meta.get("n_total")
    n_train = meta.get("n_train")
    n_test = meta.get("n_test")
    test_size = meta.get("test_size")

    if n_total is None or n_train is None or n_test is None or test_size is None:
        return ""

    train_ratio = 100 * (1 - float(test_size))
    test_ratio = 100 * float(test_size)
    return (
        f"\n  Total de moléculas: {int(n_total):,}"
        f"\n  Divisão treino / teste: {int(n_train):,} / {int(n_test):,}"
        f"  ({train_ratio:.0f}% / {test_ratio:.0f}%)"
    )


class ModelsView(BaseView):
    """View for managing ML models."""

    name = "models"
    title = "Models"
    icon = ICONS["models"]
    color = COLORS["models"]

    def __init__(self, controller: CliController) -> None:
        """Initialize the models view."""
        super().__init__(controller)
        self.selected_model_path: Path | None = None

    def _get_model_path(self) -> Path:
        """Resolve default model path from env var (read fresh on each call)."""
        return Path(
            os.environ.get("GRIMPERIUM_MODEL_PATH", "models/delta_learner_v1.joblib")
        )

    def _get_data_path(self) -> Path:
        """Resolve data path from env var (read fresh on each call)."""
        return Path(os.environ.get("GRIMPERIUM_DATA_PATH", "data/thermo_pm7.csv"))

    def _discover_models(self) -> list[dict[str, Any]]:
        """Discover all .joblib model files in the models/ directory."""
        return discover_available_models(
            active_model_path=self.controller.current_model_path,
        )

    def render(self) -> None:
        """Render the models overview."""
        self.clear_screen()
        self.show_header()

        models = self._discover_models()

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
        table.add_column("Active")
        table.add_column("Status")

        if not models:
            self.console.print(
                f"[{COLORS['muted']}]No trained models found. "
                f"Use 'Train New Model' to create one.[/{COLORS['muted']}]"
            )
        else:
            for m in models:
                mae = f"{m['mae']:.3f}" if m["mae"] is not None else "-"
                r2 = f"{m['r2']:.4f}" if m["r2"] is not None else "-"
                active = "★" if m["is_active"] else ""
                status = f"[{COLORS['success']}]{ICONS['success']} Ready[/{COLORS['success']}]"
                table.add_row(
                    m["display_name"], m["algorithm"], mae, r2, active, status
                )
            self.console.print(table)

        self.console.print()
        self.console.print(
            f"[{COLORS['muted']}]★ = Active model for predictions[/{COLORS['muted']}]"
        )
        self.console.print()

    def render_model_detail(self, metadata: dict[str, Any] | None) -> None:
        """Render detailed view for the selected model."""
        self.clear_screen()
        self.show_header()

        display_name = (
            self.selected_model_path.stem.replace("_", " ").title()
            if self.selected_model_path is not None
            else "Unknown"
        )

        if metadata is None:
            info = f"""
[bold]Name:[/bold]          {display_name}
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
            split_summary = _format_split_summary(metadata)

            info = f"""
[bold]Name:[/bold]          {display_name}
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
{split_summary}
"""

        self.console.print(
            Panel(
                info,
                title=f"[bold {COLORS['models']}]{display_name}[/bold {COLORS['models']}]",
                border_style=COLORS["models"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the models view — one per discovered model."""
        models = self._discover_models()
        options: list[MenuOption] = []
        for m in models:
            options.append(
                MenuOption(
                    label=f"View {m['display_name']}",
                    value=f"view_{m['path'].stem}",
                    icon=ICONS["models"],
                )
            )
        options.append(
            MenuOption(
                label="Train New Model",
                value="train",
            )
        )
        return options

    def get_detail_menu_options(self) -> list[MenuOption]:
        """Return menu options for model detail view."""
        options = []

        is_default = (
            self.selected_model_path is not None
            and self.selected_model_path == self.controller.current_model_path
        )
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
        if action is None or action == "back":
            if self.selected_model_path is not None:
                self.selected_model_path = None
                return "models"
            return "main"

        if action.startswith("view_"):
            stem = action[len("view_") :]
            for m in self._discover_models():
                if m["path"].stem == stem:
                    self.selected_model_path = m["path"]
                    return None
            return None

        if action == "set_default":
            if self.selected_model_path is not None:
                name = self.selected_model_path.stem.replace("_", " ").title()
                self.controller.set_model(name, model_path=self.selected_model_path)
                self.show_success(f"Active model set to: {name}")
            return None

        if action == "train":
            self._handle_train()
            return None

        if action == "test":
            self._handle_test()
            return None

        if action == "retrain":
            self._handle_retrain()
            return None

        if action == "export":
            self.show_in_development("Export")
            return None

        return None

    def _handle_train(self) -> None:
        """Train a new model and save it."""
        self.console.print()
        self.console.print(
            f"[bold {COLORS['models']}]Training new model...[/bold {COLORS['models']}] "
            "(this may take a minute)"
        )
        self.console.print()

        try:
            result = train(
                self._get_data_path(),
                return_pipeline=True,
                random_state=42,
            )
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
                "metrics": {
                    "train": train_m,
                    "test": test_m,
                },
            }
            save_model(bundle, self._get_model_path())

            gate_pass = test_m.get("gate_pass", False)
            gate_icon = ICONS["success"] if gate_pass else ICONS["error"]
            split_summary = _format_split_summary(
                {
                    "n_total": test_m.get("n_total"),
                    "n_train": train_m.get("n_samples"),
                    "n_test": test_m.get("n_samples"),
                    "test_size": test_m.get("test_size"),
                }
            )

            result_text = f"""
[bold]Training Complete![/bold]

[bold]Test Results:[/bold]
  RMSE:        {_safe_metric(test_m.get('rmse')):.3f} kcal/mol
  MAE:         {_safe_metric(test_m.get('mae')):.3f} kcal/mol
  R² Score:    {_safe_metric(test_m.get('r2')):.4f}
  Gate Pass:   {gate_icon} {'Yes' if gate_pass else 'No'}
{split_summary}
"""

            self.console.print(
                Panel(
                    result_text,
                    title=f"[bold {COLORS['success']}]Training Results[/bold {COLORS['success']}]",
                    border_style=COLORS["success"],
                    padding=(1, 2),
                )
            )
            self.show_success(f"Model trained and saved to {self._get_model_path()}")
        except Exception as e:
            self.show_error(f"Training failed: {e}")

        self.wait_for_enter()

    def _handle_retrain(self) -> None:
        """Retrain the selected model in-place with before/after comparison."""
        if self.selected_model_path is None:
            self.show_error("No model selected.")
            return

        try:
            old_meta = load_model_metadata(self.selected_model_path)
        except Exception as e:
            self.show_error(f"Could not load model metadata: {e}")
            return

        old_test = old_meta.get("metrics", {}).get("test", {})

        try:
            n_molecules = sum(1 for _ in self._get_data_path().open()) - 1
            if n_molecules < 0:
                n_molecules = 0
        except Exception as exc:
            logger.debug("Could not count molecules in data file: %s", exc)
            n_molecules = 0

        display_name = self.selected_model_path.stem.replace("_", " ").title()

        self.console.print()
        self.console.print(
            f"[bold]Retrain {display_name}?[/bold]\n"
            f"  Data: {self._get_data_path()} ({n_molecules} molecules)\n"
            f"  This will overwrite the existing model file.\n"
        )
        confirm = self.console.input("Type 'yes' to confirm: ").strip().lower()
        if confirm not in ("yes", "y"):
            self.console.print(
                f"[{COLORS['muted']}]Retrain cancelled.[/{COLORS['muted']}]"
            )
            return

        self.console.print()
        self.console.print(
            f"[bold {COLORS['models']}]Retraining {display_name}...[/bold {COLORS['models']}] "
            "(this may take a minute)"
        )

        try:
            result = train(
                self._get_data_path(),
                return_pipeline=True,
                random_state=42,
            )
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
                "metrics": {
                    "train": train_m,
                    "test": test_m,
                },
            }
            save_model(bundle, self.selected_model_path)

            # Before/after comparison table
            comp = Table(
                title="Before vs After Retrain",
                show_header=True,
                header_style=f"bold {COLORS['models']}",
                border_style=COLORS["border"],
            )
            comp.add_column("Metric")
            comp.add_column("Before", justify="right")
            comp.add_column("After", justify="right")
            comp.add_column("Delta", justify="right")
            comp.add_column("")

            def _delta_row(
                label: str,
                before: float | None,
                after: float | None,
                lower_is_better: bool = True,
            ) -> None:
                b = _safe_metric(before)
                a = _safe_metric(after)
                delta = a - b
                improved = (delta < 0) if lower_is_better else (delta > 0)
                icon = ICONS["success"] if improved else ICONS["error"]
                comp.add_row(label, f"{b:.4f}", f"{a:.4f}", f"{delta:+.4f}", icon)

            _delta_row("MAE (kcal/mol)", old_test.get("mae"), test_m.get("mae"))
            _delta_row("RMSE (kcal/mol)", old_test.get("rmse"), test_m.get("rmse"))
            _delta_row(
                "R²",
                old_test.get("r2"),
                test_m.get("r2"),
                lower_is_better=False,
            )
            _delta_row(
                "Max Error (kcal/mol)",
                old_test.get("max_error"),
                test_m.get("max_error"),
            )

            self.console.print(comp)
            gate_pass = test_m.get("gate_pass", False)
            gate_icon = ICONS["success"] if gate_pass else ICONS["error"]
            self.console.print(
                f"  Gate Pass: {gate_icon} {'Yes' if gate_pass else 'No'}"
            )
            split_summary = _format_split_summary(
                {
                    "n_total": test_m.get("n_total"),
                    "n_train": train_m.get("n_samples"),
                    "n_test": test_m.get("n_samples"),
                    "test_size": test_m.get("test_size"),
                }
            )
            if split_summary:
                self.console.print(split_summary)
            self.show_success(
                f"Model retrained and saved to {self.selected_model_path}"
            )
        except Exception as e:
            self.show_error(f"Retraining failed: {e}")

        self.wait_for_enter()

    def _handle_test(self) -> None:
        """Display test metrics for the selected model."""
        model_path = self.selected_model_path or self._get_model_path()
        try:
            metadata = load_model_metadata(model_path)
        except FileNotFoundError:
            self.show_error("No trained model found. Train first.")
            self.wait_for_enter()
            return

        display_name = model_path.stem.replace("_", " ").title()
        test_metrics = metadata.get("metrics", {}).get("test", {})
        gate_pass = test_metrics.get("gate_pass", False)
        gate_icon = ICONS["success"] if gate_pass else ICONS["error"]
        split_summary = _format_split_summary(metadata)

        result_text = f"""
[bold]Model:[/bold]       {display_name}
[bold]Version:[/bold]     {metadata.get('version', 'unknown')}
[bold]Trained at:[/bold]  {metadata.get('trained_at', 'unknown')}

[bold]Test Metrics:[/bold]
  RMSE:        {_safe_metric(test_metrics.get('rmse')):.3f} kcal/mol
  MAE:         {_safe_metric(test_metrics.get('mae')):.3f} kcal/mol
  R² Score:    {_safe_metric(test_metrics.get('r2')):.4f}
  MAPE:        {_safe_metric(test_metrics.get('mape')):.2f}%
  Max Error:   {_safe_metric(test_metrics.get('max_error')):.3f} kcal/mol
{split_summary}
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
            if self.selected_model_path is not None:
                try:
                    metadata: dict[str, Any] | None = load_model_metadata(
                        self.selected_model_path
                    )
                except Exception as exc:
                    logger.debug("Could not load model metadata: %s", exc)
                    metadata = None
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

            next_view = self.handle_action(result)
            if next_view:
                return next_view
