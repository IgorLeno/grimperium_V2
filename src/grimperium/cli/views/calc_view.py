"""
Calc view for GRIMPERIUM CLI.

Handles molecular property predictions.
"""

import hashlib
import os
from pathlib import Path
from typing import TYPE_CHECKING

from rdkit import Chem
from rich.panel import Panel
from rich.table import Table

from grimperium.calculation.contracts.enums import PropertyRole
from grimperium.calculation.contracts.models import (
    MoleculeCalculationResult,
    PropertyEstimate,
)
from grimperium.calculation.methods import (
    CalculationMethodDefinition,
    get_calculation_method,
    list_calculation_methods,
)
from grimperium.calculation.runners import SemiempiricalFormationEnthalpyRunner
from grimperium.cli.calc_pipeline import (
    CalcPipelineError,
    CalcPipelineResult,
    run_single_molecule_prediction,
)
from grimperium.cli.menu import MenuOption, show_back_menu, show_menu, text_input
from grimperium.cli.mock_data import PredictionResult
from grimperium.cli.model_compatibility import (
    ModelCompatibilityError,
    validate_model_for_method,
)
from grimperium.cli.styles import COLORS, ICONS
from grimperium.cli.views.base_view import BaseView
from grimperium.crest_pm7.config import PM7Config

if TYPE_CHECKING:
    from grimperium.cli.controller import CliController

KCAL_TO_KJ: float = 4.184  # exact per IUPAC definition


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
        self.last_calculation_result: MoleculeCalculationResult | None = None
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
        table.add_row(
            "H298 (PM7)",
            f"[bold {COLORS['success']}]{result.h298_pm7:.2f} kcal/mol[/bold {COLORS['success']}]",
        )
        table.add_row(
            "Delta Correction",
            f"{result.delta_correction:.2f} kcal/mol",
        )
        h298_kj = result.h298_corrected * KCAL_TO_KJ
        table.add_row(
            "H298 (Corrected)",
            f"[bold {COLORS['success']}]{result.h298_corrected:.2f} kcal/mol  "
            f"({h298_kj:.2f} kJ/mol)[/bold {COLORS['success']}]",
        )
        table.add_row("Model Used", result.model_name)
        table.add_row("Model Version", result.model_version)
        table.add_row("Conformers", str(result.n_conformers))
        exec_min = result.execution_time / 60
        table.add_row(
            "Execution Time",
            f"{exec_min:.2f} min",
        )

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
        table.add_column("H298 PM7", justify="right")
        table.add_column("Delta", justify="right")
        table.add_column("H298 Corrected", justify="right")
        table.add_column("Model", justify="center")
        table.add_column("Conformers", justify="right")
        table.add_column("Time (min)", justify="right")

        for i, result in enumerate(reversed(self.history[-10:]), 1):
            table.add_row(
                str(i),
                (
                    result.smiles[:15] + "..."
                    if len(result.smiles) > 15
                    else result.smiles
                ),
                f"{result.h298_pm7:.2f}",
                f"{result.delta_correction:.2f}",
                f"{result.h298_corrected:.2f}",
                result.model_name.split("_")[0] if result.model_name else "-",
                str(result.n_conformers),
                f"{result.execution_time / 60:.2f}",
            )

        self.console.print(table)
        self.console.print()

    def render_available_methods(self) -> None:
        """Render available standard enthalpy calculation methods."""
        methods = list_calculation_methods("standard_enthalpy_of_formation")
        table = Table(
            title="Calculation Methods",
            show_header=True,
            header_style=f"bold {COLORS['calc']}",
            border_style=COLORS["border"],
        )
        table.add_column("Method ID", style="bold", min_width=30, no_wrap=True)
        table.add_column("Name")
        table.add_column("Model")
        table.add_column("Conformer Strategy")
        table.add_column("xTB")

        for method in methods:
            model = "Required" if method.model_requirement.model_required else "No"
            xtb = (
                "Default"
                if method.xtb.optional and method.xtb.enabled_by_default
                else "Disabled"
            )
            table.add_row(
                method.method_id,
                method.display_name,
                model,
                method.conformer_selection.strategy,
                xtb,
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

    def _resolve_model_path(self) -> Path | None:
        """
        Resolve the path to the current model file.

        Checks in order:
        1. Controller's current_model_path (if set and exists)
        2. GRIMPERIUM_MODEL_PATH environment variable
        3. Returns None if no valid path found

        Returns:
            Path to model file, or None if not found
        """
        # Check controller path
        if (
            self.controller.current_model_path
            and self.controller.current_model_path.exists()
        ):
            return self.controller.current_model_path

        # Check environment variable
        env_path = os.environ.get("GRIMPERIUM_MODEL_PATH")
        if env_path:
            path = Path(env_path)
            if path.exists():
                return path

        return None

    def _select_method(self) -> CalculationMethodDefinition | None:
        """Prompt for a standard enthalpy calculation method."""
        methods = list_calculation_methods("standard_enthalpy_of_formation")
        selected = show_menu(
            [
                MenuOption(
                    label=method.display_name,
                    value=method.method_id,
                    icon=ICONS["calc"],
                )
                for method in methods
            ],
            title="Calculation Method",
        )
        if selected is None:
            return None
        return get_calculation_method(
            selected,
            property_id="standard_enthalpy_of_formation",
        )

    def _resolve_required_model(
        self,
        method: CalculationMethodDefinition,
    ) -> Path | None:
        """Resolve and validate the live-session model required by a method."""
        if not method.model_requirement.model_required:
            return None

        model_path = self.controller.current_model_path
        if model_path is None:
            self.show_error(
                "Select a compatible model in Models before running this method."
            )
            return None
        if not model_path.exists():
            self.show_error(f"Selected model file no longer exists: {model_path}")
            return None

        try:
            validate_model_for_method(model_path, method)
        except ModelCompatibilityError as exc:
            self.show_error(f"Selected model is not compatible: {exc}")
            return None

        return model_path

    def _pm7_config_from_settings(self) -> PM7Config:
        config = PM7Config()
        self.controller.settings_manager.apply_to_pm7_config(config)
        return config

    def render_method_a_result(
        self,
        result: MoleculeCalculationResult,
        *,
        units: str = "both",
    ) -> None:
        """Render a compact Method A canonical result summary."""
        table = Table(
            show_header=True,
            header_style=f"bold {COLORS['calc']}",
            border_style=COLORS["calc"],
        )
        table.add_column("Hamiltonian")
        if units in ("kcal/mol", "both"):
            table.add_column("H298", justify="right")
        if units in ("kJ/mol", "both"):
            table.add_column("H298 (kJ/mol)", justify="right")

        for estimate in result.estimates:
            if estimate.value_kcal_mol is None:
                continue
            kj_value = (
                estimate.value_kj_mol
                if estimate.value_kj_mol is not None
                else estimate.value_kcal_mol * KCAL_TO_KJ
            )
            row = [estimate.hamiltonian or "-"]
            if units in ("kcal/mol", "both"):
                row.append(f"{estimate.value_kcal_mol:.2f} kcal/mol")
            if units in ("kJ/mol", "both"):
                row.append(f"{kj_value:.2f} kJ/mol")
            table.add_row(*row)

        self.console.print(
            Panel(
                table,
                title=f"[bold {COLORS['calc']}]{ICONS['success']} Calculation Result[/bold {COLORS['calc']}]",
                border_style=COLORS["success"],
                padding=(1, 1),
            )
        )
        self.console.print()

    def _run_method_a(
        self,
        smiles: str,
        method: CalculationMethodDefinition,
        units: str = "both",
    ) -> bool:
        mol_id = "calc_" + hashlib.sha1(smiles.encode()).hexdigest()[:6]
        runner = SemiempiricalFormationEnthalpyRunner(
            config=self._pm7_config_from_settings(),
            xtb_enabled=method.xtb.enabled_by_default,
        )
        try:
            result = runner.calculate_single_smiles(
                smiles,
                molecule_id=mol_id,
                name=mol_id,
            )
        except Exception as exc:
            self.show_error(f"Calculation failed: {exc}")
            return True

        self.last_calculation_result = result
        self.last_result = self._method_a_prediction_summary(result, method)
        self.console.print(
            f"[green]✓ Calculation complete for SMILES: {smiles}[/green]"
        )
        self.render_method_a_result(result, units=units)
        return True

    def _method_a_prediction_summary(
        self,
        result: MoleculeCalculationResult,
        method: CalculationMethodDefinition,
    ) -> PredictionResult:
        estimate = self._select_method_a_display_estimate(result)
        value_kcal = (
            estimate.value_kcal_mol
            if estimate.value_kcal_mol is not None
            else estimate.value.to_kcal_mol()
        )
        execution_time = sum(
            stage.execution_time_s or 0.0 for stage in result.stage_executions
        )
        return PredictionResult(
            smiles=result.molecule.smiles,
            h298_pm7=value_kcal,
            delta_correction=0.0,
            h298_corrected=value_kcal,
            model_name="Method A",
            model_version=method.version,
            execution_time=execution_time,
            n_conformers=len(result.conformers),
        )

    @staticmethod
    def _select_method_a_display_estimate(
        result: MoleculeCalculationResult,
    ) -> PropertyEstimate:
        final_estimates = [
            estimate
            for estimate in result.estimates
            if estimate.role is PropertyRole.FINAL
        ]
        for hamiltonian in ("PM7", "PM3", "AM1"):
            for estimate in final_estimates:
                if estimate.hamiltonian == hamiltonian:
                    return estimate
        if final_estimates:
            return final_estimates[0]
        raise ValueError("Method A result has no final property estimate")

    def _prediction_result_from_pipeline(
        self,
        smiles: str,
        pipeline_result: CalcPipelineResult,
    ) -> PredictionResult:
        return PredictionResult(
            smiles=smiles,
            h298_pm7=pipeline_result.h298_pm7,
            delta_correction=pipeline_result.delta_correction,
            h298_corrected=pipeline_result.h298_corrected,
            model_name=self.controller.current_model,
            model_version=pipeline_result.model_version,
            execution_time=pipeline_result.execution_time,
            n_conformers=pipeline_result.n_conformers,
        )

    def _run_method_b(
        self,
        smiles: str,
        method: CalculationMethodDefinition,
    ) -> bool:
        model_path = self._resolve_required_model(method)
        if model_path is None:
            return True

        def progress_update(msg: str) -> None:
            self.console.print(msg)

        mol_id = "calc_" + hashlib.sha1(smiles.encode()).hexdigest()[:6]
        try:
            pipeline_result = run_single_molecule_prediction(
                smiles,
                mol_id,
                model_path,
                self._pm7_config_from_settings(),
                progress_update,
            )
        except CalcPipelineError as exc:
            self.show_error(str(exc))
            return True
        except Exception as exc:
            self.show_error(f"Unexpected error: {exc}")
            return True

        result = self._prediction_result_from_pipeline(smiles, pipeline_result)
        self.last_result = result
        self.history.append(result)
        self.console.print(
            f"[green]✓ Calculation complete for SMILES: {smiles}[/green]"
        )
        self.render_result(result)
        return True

    def _select_units(self) -> str:
        selected = show_menu(
            [
                MenuOption(label="Both", value="both"),
                MenuOption(label="kcal/mol", value="kcal/mol"),
                MenuOption(label="kJ/mol", value="kJ/mol"),
            ],
            title="Display Units",
        )
        return selected or "both"

    def render_review_panel(
        self,
        *,
        method: CalculationMethodDefinition,
        smiles: str,
        units: str,
    ) -> None:
        model_status = (
            str(self.controller.current_model_path)
            if method.model_requirement.model_required
            and self.controller.current_model_path is not None
            else "Not required"
        )
        if (
            method.model_requirement.model_required
            and self.controller.current_model_path is None
        ):
            model_status = "Required but not selected"
        xtb_status = (
            "enabled by default"
            if method.xtb.optional and method.xtb.enabled_by_default
            else "not enabled by default"
        )
        table = Table(show_header=False, border_style=COLORS["calc"])
        table.add_column("Field", style="bold")
        table.add_column("Value")
        table.add_row("Property", method.property_name)
        table.add_row("Method", method.display_name)
        table.add_row("SMILES", smiles)
        table.add_row("Model", model_status)
        table.add_row("xTB", xtb_status)
        table.add_row("Units", units)
        self.console.print(
            Panel(
                table,
                title=f"[bold {COLORS['calc']}]Review Calculation[/bold {COLORS['calc']}]",
                border_style=COLORS["calc"],
                padding=(1, 1),
            )
        )

    def do_prediction(self) -> bool:
        """
        Perform a prediction interaction using the real pipeline.

        Returns:
            True to continue, False to go back
        """
        self.render()

        method = self._select_method()
        if method is None:
            return True

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

        units = self._select_units()
        self.render_review_panel(method=method, smiles=smiles, units=units)
        self.console.print()
        self.console.print()

        if method.method_id == "semiempirical_am1_pm3_pm7":
            return self._run_method_a(smiles, method, units=units)
        return self._run_method_b(smiles, method)

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
                    label="Calculation Methods",
                    value="methods",
                    icon=ICONS["settings"],
                ),
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

        if action == "methods":
            self.clear_screen()
            self.show_header()
            self.render_available_methods()
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
                _h298_kj = self.last_result.h298_corrected * KCAL_TO_KJ
                self.console.print(
                    f"[{COLORS['muted']}]Last prediction: {self.last_result.smiles} → "
                    f"H298={self.last_result.h298_corrected:.2f} kcal/mol "
                    f"({_h298_kj:.2f} kJ/mol)[/{COLORS['muted']}]"
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
