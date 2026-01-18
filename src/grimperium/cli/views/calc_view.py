"""
Calc view for GRIMPERIUM CLI.

Handles molecular property predictions.
"""

from typing import TYPE_CHECKING

from rdkit import Chem
from rich.panel import Panel
from rich.table import Table

from grimperium.cli.menu import MenuOption, show_back_menu, text_input
from grimperium.cli.mock_data import PredictionResult, mock_predict
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController


class CalcView(BaseView):
    """View for molecular property predictions."""

    name = "calc"
    title = "Prediction Engine"
    icon = ICONS["calc"]
    color = COLORS["calc"]

    def __init__(self, controller: "CliController") -> None:
        """Initialize the calc view."""
        super().__init__(controller)
        self.last_result: PredictionResult | None = None
        self.history: list[PredictionResult] = []

    def render(self) -> None:
        """Render the calc view."""
        self.clear_screen()
        self.show_header()

        intro = f"""
[bold]Molecular Property Prediction[/bold]

Enter a SMILES string to predict the Heat of Formation (HOF)
using the Delta-Learning model.

[{COLORS['muted']}]Current Model:[/{COLORS['muted']}] [{COLORS['calc']}]{self.controller.current_model}[/{COLORS['calc']}]
"""
        self.console.print(
            Panel(
                intro,
                border_style=COLORS["calc"],
                padding=(1, 2),
            )
        )
        self.console.print()

    def render_result(self, result: PredictionResult) -> None:
        """Render a prediction result."""
        self.console.print()

        # Result table
        table = Table(
            show_header=False,
            border_style=COLORS["calc"],
            padding=(0, 2),
        )
        table.add_column("Property", style="bold")
        table.add_column("Value")

        table.add_row("SMILES", f"[bold]{result.smiles}[/bold]")
        table.add_row("Molecule", result.molecule_name)
        table.add_row("Property", result.property_name)
        table.add_row(
            "Predicted Value",
            f"[bold {COLORS['success']}]{result.predicted_value:.2f} {result.unit}[/bold {COLORS['success']}]",
        )
        table.add_row(
            "Confidence",
            f"[{COLORS['success']}]{result.confidence * 100:.1f}%[/{COLORS['success']}]",
        )
        table.add_row("Model Used", result.model_used)

        self.console.print(
            Panel(
                table,
                title=f"[bold {COLORS['calc']}]{ICONS['success']} Prediction Result[/bold {COLORS['calc']}]",
                border_style=COLORS["success"],
                padding=(1, 1),
            )
        )
        self.console.print()

    def render_history(self) -> None:
        """Render prediction history."""
        if not self.history:
            self.console.print(
                f"[{COLORS['muted']}]No predictions yet.[/{COLORS['muted']}]"
            )
            return

        table = Table(
            title="Prediction History",
            show_header=True,
            header_style=f"bold {COLORS['calc']}",
            border_style=COLORS["border"],
        )
        table.add_column("#", justify="right")
        table.add_column("SMILES")
        table.add_column("Molecule")
        table.add_column("HOF (kcal/mol)", justify="right")
        table.add_column("Confidence", justify="right")

        for i, result in enumerate(reversed(self.history[-10:]), 1):
            table.add_row(
                str(i),
                (
                    result.smiles[:20] + "..."
                    if len(result.smiles) > 20
                    else result.smiles
                ),
                result.molecule_name,
                f"{result.predicted_value:.2f}",
                f"{result.confidence * 100:.1f}%",
            )

        self.console.print(table)
        self.console.print()

    def validate_smiles(self, smiles: str) -> bool | str:
        """
        Validate SMILES string using RDKit.

        Args:
            smiles: SMILES string to validate.

        Returns:
            True if valid, error message string if invalid.
        """
        # Strip whitespace
        smiles_clean = smiles.strip() if smiles else ""

        # Check for empty string
        if not smiles_clean:
            return "Please enter a valid SMILES string"

        # Use RDKit to validate SMILES
        try:
            mol = Chem.MolFromSmiles(smiles_clean)
            if mol is None:
                return f"Invalid SMILES: '{smiles_clean}' cannot be parsed by RDKit"
            return True
        except Exception as e:
            return f"Error validating SMILES: {str(e)}"

    def validate_molecules(
        self, molecules: list[dict[str, str]]
    ) -> tuple[list[dict[str, str]], str]:
        """
        Validate a list of molecules and remove duplicates/invalid SMILES.

        Args:
            molecules: List of dicts with 'smiles' and optional 'name' keys.

        Returns:
            Tuple of (valid_molecules, summary_message).
        """
        if not molecules:
            return [], "No molecules provided"

        valid: list[dict[str, str]] = []
        seen_smiles: set[str] = set()
        rejected_count = 0
        duplicate_count = 0

        for mol_dict in molecules:
            smiles = mol_dict.get("smiles", "").strip()

            # Skip empty SMILES
            if not smiles:
                rejected_count += 1
                continue

            # Check for duplicates
            if smiles in seen_smiles:
                duplicate_count += 1
                continue

            # Validate with RDKit
            validation_result = self.validate_smiles(smiles)
            if validation_result is True:
                valid.append(mol_dict)
                seen_smiles.add(smiles)
            else:
                rejected_count += 1

        # Build summary
        summary_parts = [f"{len(valid)} valid molecule(s)"]
        if rejected_count > 0:
            summary_parts.append(f"{rejected_count} rejected (invalid SMILES)")
        if duplicate_count > 0:
            summary_parts.append(f"{duplicate_count} duplicate(s) removed")

        summary = ", ".join(summary_parts)
        return valid, summary

    def do_prediction(self) -> bool:
        """
        Perform a prediction interaction.

        Returns:
            True to continue, False to go back
        """
        self.render()

        # Get SMILES input
        smiles = text_input(
            message="Enter SMILES",
            validate=self.validate_smiles,
        )

        if smiles is None:  # Ctrl+C
            return True  # Stay in calc view

        smiles = smiles.strip()
        if not smiles:
            return True

        # Perform prediction
        self.console.print()
        self.console.print(
            f"[{COLORS['muted']}]Calculating prediction...[/{COLORS['muted']}]"
        )

        result = mock_predict(smiles, self.controller.current_model)
        self.last_result = result
        self.history.append(result)

        self.console.print(
            f"[green]✓ Calculation complete for SMILES: {smiles}[/green]"
        )
        self.render_result(result)

        return True

    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for the calc view."""
        options = [
            MenuOption(
                label="Predict New Molecule",
                value="predict",
                icon=ICONS["calc"],
            ),
        ]

        if self.history:
            options.append(
                MenuOption(
                    label="View History",
                    value="history",
                    icon="📜",
                )
            )

        options.extend(
            [
                MenuOption(
                    label="Batch Processing",
                    value="batch",
                    icon="📁",
                    disabled=True,
                    disabled_reason="In Development",
                ),
                MenuOption(
                    label="Export Results",
                    value="export",
                    icon="💾",
                    disabled=True,
                    disabled_reason="In Development",
                ),
            ]
        )

        return options

    def handle_action(self, action: str) -> str | None:
        """Handle menu actions."""
        if action == "back":
            return "main"

        if action == "predict":
            self.do_prediction()
            return None

        if action == "history":
            self.clear_screen()
            self.show_header()
            self.render_history()
            self.wait_for_enter()
            return None

        # Handle in-development features
        if action in ["batch", "export"]:
            self.show_in_development(action.title())
            return None

        return None

    def run(self) -> str | None:
        """Run the calc view interaction loop."""
        while True:
            self.render()

            # Show last result if available
            if self.last_result:
                self.console.print(
                    f"[{COLORS['muted']}]Last prediction: {self.last_result.smiles} → "
                    f"{self.last_result.predicted_value:.2f} {self.last_result.unit}[/{COLORS['muted']}]"
                )
                self.console.print()

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
