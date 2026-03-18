"""
Settings Manager for GRIMPERIUM CLI.

Manages CREST, MOPAC, and xTB configuration with interactive menus.
"""

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import questionary
from prompt_toolkit.styles import Style
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from grimperium.cli.styles import COLORS

if TYPE_CHECKING:
    from grimperium.crest_pm7.config import PM7Config


@dataclass
class CRESTSettings:
    """CREST configuration options for conformer search.

    Attributes:
        v3: Use iMTD-GC v3 algorithm (recommended for better conformer sampling).
        nci: NCI mode for molecular complexes and weak interactions.
        crest_method: Choose CREST quantum method (gfn2/gfnff/gfn2//gfnff).
        quick_mode: Speed/accuracy tradeoff (off/quick/squick/mquick).
        ewin: Energy window in kcal/mol for conformer selection.
        rthr: RMSD threshold in Angstroms for geometry deduplication.
        optlev: Optimization level (loose/normal/tight/vtight/extreme).
        threads: Number of parallel threads.
    """

    # Valid options for dropdown menus
    CREST_METHOD_OPTIONS: ClassVar[list[str]] = ["gfn2", "gfnff", "gfn2//gfnff"]
    QUICK_MODE_OPTIONS: ClassVar[list[str]] = ["off", "quick", "squick", "mquick"]
    OPTLEV_CHOICES: ClassVar[list[str]] = [
        "loose",
        "normal",
        "tight",
        "vtight",
        "extreme",
    ]

    # Settings
    v3: bool = True
    nci: bool = False
    crest_method: str = "gfn2"
    quick_mode: str = "off"
    ewin: float = 5.0
    rthr: float = 0.125
    optlev: str = "normal"
    threads: int = 4


@dataclass
class MOPACSettings:
    """MOPAC/PM7 configuration options.

    Attributes:
        precise: 100x tighter SCF convergence.
        scfcrt: SCF convergence threshold.
        itry: Maximum SCF iterations.
        pulay: PULAY acceleration for SCF.
        prtall: Verbose output.
        archive: Save archive file.
    """

    precise: bool = False
    scfcrt: float = 1.0e-4
    itry: int = 1000
    pulay: bool = False
    prtall: bool = False
    archive: bool = False


@dataclass
class xTBSettings:
    """xTB pre-optimization configuration.

    Attributes:
        preopt: Pre-optimize structures with xTB before CREST.
        timeout_seconds: Timeout for xTB pre-optimization.
    """

    preopt: bool = True
    timeout_seconds: int = 300


# Questionary style for settings menus
SETTINGS_STYLE = Style.from_dict(
    {
        "qmark": COLORS["settings"],
        "question": "bold",
        "answer": COLORS["success"],
        "pointer": f"{COLORS['settings']} bold",
        "highlighted": f"{COLORS['settings']} bold",
        "selected": COLORS["success"],
        "separator": COLORS["muted"],
        "instruction": COLORS["muted"],
        "text": "",
        "disabled": COLORS["in_dev"],
    }
)


@dataclass
class SettingsManager:
    """Manage CREST, MOPAC, and xTB settings with interactive menus.

    Provides dataclasses for settings, help text, interactive menus,
    and serialization to/from dict for CSV storage.

    Example:
        >>> manager = SettingsManager()
        >>> manager.crest_settings.ewin = 6.0
        >>> settings_dict = manager.to_dict()
        >>> settings_dict['crest_ewin']
        6.0
    """

    crest: CRESTSettings = field(default_factory=CRESTSettings)
    mopac: MOPACSettings = field(default_factory=MOPACSettings)
    xtb: xTBSettings = field(default_factory=xTBSettings)
    console: Console = field(default_factory=Console)

    HELP_TEXT: ClassVar[dict[str, str]] = {
        "crest_v3": (
            "Use iMTD-GC v3 algorithm (recommended for better conformer sampling)"
        ),
        "crest_nci": "NCI mode for molecular complexes and weak interactions",
        "crest_method": (
            "Choose CREST quantum method: gfn2 (default, balanced), "
            "gfnff (faster), gfn2//gfnff (two-step refinement)"
        ),
        "crest_quick_mode": (
            "Choose speed/accuracy tradeoff: off (full), quick (fast), "
            "squick (super-fast), mquick (fastest)"
        ),
        "crest_ewin": (
            "Energy window in kcal/mol for conformer selection "
            "(higher = more conformers)"
        ),
        "crest_rthr": "RMSD threshold in Angstroms for geometry deduplication",
        "crest_optlev": (
            "Optimization level: loose/normal/tight/vtight/extreme "
            "(tighter = more thorough)"
        ),
        "crest_threads": (
            "Number of parallel threads (higher = faster, uses more CPU)"
        ),
        "mopac_precise": (
            "100x tighter SCF convergence - slower but more accurate "
            "for critical molecules"
        ),
        "mopac_scfcrt": "SCF convergence threshold - lower = tighter, more accurate",
        "mopac_itry": "Maximum SCF iterations - higher allows more difficult convergence",
        "mopac_pulay": (
            "PULAY acceleration for SCF (experimental, may help difficult cases)"
        ),
        "mopac_prtall": "Verbose output - useful for debugging",
        "mopac_archive": "Save archive file - useful for post-processing",
        "xtb_preopt": (
            "Pre-optimize structure with xTB before CREST "
            "(recommended by CREST docs)"
        ),
    }

    @staticmethod
    def _parse_bool(value: str | int | bool) -> bool:
        """Parse boolean from string, int, or bool value.

        Handles CSV round-trip correctly: value → to_dict() → CSV → from_dict() → value

        Args:
            value: Value to parse (string, int, or bool).

        Returns:
            Parsed boolean value.

        Raises:
            ValueError: If string value cannot be parsed.
            TypeError: If value type is not supported.
        """
        if isinstance(value, bool):
            return value

        if isinstance(value, str):
            lower_val = value.strip().lower()
            if lower_val in ("true", "yes", "1", "on"):
                return True
            elif lower_val in ("false", "no", "0", "off"):
                return False
            else:
                raise ValueError(f"Cannot parse boolean from: '{value}'")

        if isinstance(value, int):
            return value != 0

        raise TypeError(f"Cannot parse boolean from {type(value).__name__}")

    def to_dict(self) -> dict[str, Any]:
        """Convert all settings to a flat dictionary.

        Returns:
            Dictionary with prefixed keys for all settings.
        """
        return {
            "crest_v3": self.crest.v3,
            "crest_nci": self.crest.nci,
            "crest_method": self.crest.crest_method,
            "crest_quick_mode": self.crest.quick_mode,
            "crest_ewin": self.crest.ewin,
            "crest_rthr": self.crest.rthr,
            "crest_optlev": self.crest.optlev,
            "crest_threads": self.crest.threads,
            "mopac_precise": self.mopac.precise,
            "mopac_scfcrt": self.mopac.scfcrt,
            "mopac_itry": self.mopac.itry,
            "mopac_pulay": self.mopac.pulay,
            "mopac_prtall": self.mopac.prtall,
            "mopac_archive": self.mopac.archive,
            "crest_xtb_preopt": self.xtb.preopt,
        }

    def from_dict(self, settings_dict: dict[str, Any]) -> None:
        """Load settings from a dictionary.

        Args:
            settings_dict: Dictionary with prefixed keys for settings.
        """
        # Backward compatibility: convert old boolean toggles to new dropdown values
        if "crest_gfnff" in settings_dict and "crest_method" not in settings_dict:
            if self._parse_bool(settings_dict["crest_gfnff"]):
                settings_dict["crest_method"] = "gfnff"
            else:
                settings_dict["crest_method"] = "gfn2"

        if "crest_quick" in settings_dict and "crest_quick_mode" not in settings_dict:
            if self._parse_bool(settings_dict["crest_quick"]):
                settings_dict["crest_quick_mode"] = "quick"
            else:
                settings_dict["crest_quick_mode"] = "off"

        # CREST boolean fields
        for key, attr in [
            ("crest_v3", "v3"),
            ("crest_nci", "nci"),
        ]:
            if key in settings_dict:
                try:
                    setattr(self.crest, attr, self._parse_bool(settings_dict[key]))
                except (ValueError, TypeError):
                    pass

        # CREST numeric fields
        if "crest_ewin" in settings_dict:
            try:
                self.crest.ewin = float(settings_dict["crest_ewin"])
            except (ValueError, TypeError):
                pass

        if "crest_rthr" in settings_dict:
            try:
                self.crest.rthr = float(settings_dict["crest_rthr"])
            except (ValueError, TypeError):
                pass

        if "crest_optlev" in settings_dict:
            val = str(settings_dict["crest_optlev"]).strip().lower()
            if val in CRESTSettings.OPTLEV_CHOICES:
                self.crest.optlev = val

        if "crest_method" in settings_dict:
            val = str(settings_dict["crest_method"]).strip().lower()
            if val in CRESTSettings.CREST_METHOD_OPTIONS:
                self.crest.crest_method = val

        if "crest_quick_mode" in settings_dict:
            val = str(settings_dict["crest_quick_mode"]).strip().lower()
            if val in CRESTSettings.QUICK_MODE_OPTIONS:
                self.crest.quick_mode = val

        if "crest_threads" in settings_dict:
            try:
                self.crest.threads = int(settings_dict["crest_threads"])
            except (ValueError, TypeError):
                pass

        # MOPAC boolean fields
        for key, attr in [
            ("mopac_precise", "precise"),
            ("mopac_pulay", "pulay"),
            ("mopac_prtall", "prtall"),
            ("mopac_archive", "archive"),
        ]:
            if key in settings_dict:
                try:
                    setattr(self.mopac, attr, self._parse_bool(settings_dict[key]))
                except (ValueError, TypeError):
                    pass

        # MOPAC numeric fields
        if "mopac_scfcrt" in settings_dict:
            try:
                self.mopac.scfcrt = float(settings_dict["mopac_scfcrt"])
            except (ValueError, TypeError):
                pass

        if "mopac_itry" in settings_dict:
            try:
                self.mopac.itry = int(settings_dict["mopac_itry"])
            except (ValueError, TypeError):
                pass

        # xTB fields
        if "crest_xtb_preopt" in settings_dict:
            try:
                self.xtb.preopt = self._parse_bool(settings_dict["crest_xtb_preopt"])
            except (ValueError, TypeError):
                pass

        if "crest_xtb_timeout_seconds" in settings_dict:
            try:
                self.xtb.timeout_seconds = int(
                    settings_dict["crest_xtb_timeout_seconds"]
                )
            except (ValueError, TypeError):
                pass

    @staticmethod
    def _optlev_to_level(optlev: str) -> int:
        """Convert CREST optlev label to numeric level.

        Args:
            optlev: Optimization level label (loose/normal/tight/vtight/extreme)

        Returns:
            Numeric level (0, 1, or 2)
        """
        mapping = {
            "loose": 0,
            "normal": 1,
            "tight": 2,
            "vtight": 2,
            "extreme": 2,
        }
        return mapping.get(optlev, 1)

    def apply_to_pm7_config(self, config: "PM7Config") -> None:
        """Apply current settings to a PM7Config instance.

        Args:
            config: PM7Config instance to update
        """
        config.energy_window = self.crest.ewin
        config.crest_rmsd_threshold = self.crest.rthr
        config.crest_opt_level = self._optlev_to_level(self.crest.optlev)
        config.crest_optlev_label = self.crest.optlev
        config.crest_threads = self.crest.threads
        config.crest_method = self.crest.crest_method
        config.crest_quick_mode = self.crest.quick_mode
        config.crest_v3 = self.crest.v3
        config.crest_nci = self.crest.nci
        config.xtb_preopt = self.xtb.preopt
        config.mopac_precise_scf = self.mopac.precise
        config.mopac_scf_threshold = self.mopac.scfcrt

    def reset_crest(self) -> None:
        """Reset CREST settings to defaults."""
        self.crest = CRESTSettings()
        self.console.print("[green]✓ CREST settings reset to defaults[/green]")

    def reset_mopac(self) -> None:
        """Reset MOPAC settings to defaults."""
        self.mopac = MOPACSettings()
        self.console.print("[green]✓ MOPAC settings reset to defaults[/green]")

    def reset_xtb(self) -> None:
        """Reset xTB settings to defaults."""
        self.xtb = xTBSettings()

    def reset_all(self) -> None:
        """Reset all settings to defaults."""
        self.reset_crest()
        self.reset_mopac()
        self.reset_xtb()

    # ── Persistence ────────────────────────────────────────────────────

    @staticmethod
    def default_config_path() -> Path:
        """Return the default path for the settings config file.

        The file is stored in the current working directory (project root).

        Returns:
            Path to grimperium_settings.json.
        """
        return Path.cwd() / "grimperium_settings.json"

    def save_to_file(self, path: Path | None = None) -> bool:
        """Persist current settings to a JSON file.

        Args:
            path: Target file path. Defaults to default_config_path().

        Returns:
            True on success, False on failure.
        """
        target = path or self.default_config_path()
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(json.dumps(self.to_dict(), indent=2, ensure_ascii=False))
            return True
        except OSError:
            return False

    def load_from_file(self, path: Path | None = None) -> bool:
        """Load settings from a JSON file.

        If the file does not exist or contains invalid JSON, settings
        remain unchanged and False is returned.

        Args:
            path: Source file path. Defaults to default_config_path().

        Returns:
            True if settings were loaded, False otherwise.
        """
        target = path or self.default_config_path()
        if not target.is_file():
            return False
        try:
            data = json.loads(target.read_text())
        except (json.JSONDecodeError, OSError):
            return False
        self.from_dict(data)
        return True

    def delete_config_file(self, path: Path | None = None) -> bool:
        """Remove the persisted config file.

        Args:
            path: File to delete. Defaults to default_config_path().

        Returns:
            True if file was removed, False if it didn't exist.
        """
        target = path or self.default_config_path()
        try:
            target.unlink()
            return True
        except FileNotFoundError:
            return False

    def _status_icon(self, value: bool) -> str:
        """Return status icon for boolean value."""
        return f"[{COLORS['success']}]✓[/{COLORS['success']}]" if value else "○"

    def show_crest_summary(self) -> Table:
        """Create a summary table for CREST settings.

        Returns:
            Rich Table with current CREST settings.
        """
        table = Table(
            show_header=False,
            box=None,
            padding=(0, 2),
            expand=True,
        )
        table.add_column("Setting", style="bold")
        table.add_column("Value", style=COLORS["settings"])

        table.add_row("v3 Algorithm", self._status_icon(self.crest.v3))
        table.add_row("NCI Mode", self._status_icon(self.crest.nci))
        table.add_row("CREST Method", self.crest.crest_method.upper())
        table.add_row("Quick Mode", self.crest.quick_mode.upper())
        table.add_row("Energy Window", f"{self.crest.ewin} kcal/mol")
        table.add_row("RMSD Threshold", f"{self.crest.rthr} Å")
        table.add_row("Optimization Level", self.crest.optlev)
        table.add_row("Threads", str(self.crest.threads))
        table.add_row("", "")  # Separator
        table.add_row("⚡ xTB Pre-optimization", self._status_icon(self.xtb.preopt))

        return table

    def show_mopac_summary(self) -> Table:
        """Create a summary table for MOPAC settings.

        Returns:
            Rich Table with current MOPAC settings.
        """
        table = Table(
            show_header=False,
            box=None,
            padding=(0, 2),
            expand=True,
        )
        table.add_column("Setting", style="bold")
        table.add_column("Value", style=COLORS["settings"])

        table.add_row("Precise SCF", self._status_icon(self.mopac.precise))
        table.add_row("SCF Threshold", f"{self.mopac.scfcrt:.1e}")
        table.add_row("Max Iterations", str(self.mopac.itry))
        table.add_row("PULAY Acceleration", self._status_icon(self.mopac.pulay))
        table.add_row("Verbose Output", self._status_icon(self.mopac.prtall))
        table.add_row("Archive Output", self._status_icon(self.mopac.archive))

        return table

    def show_xtb_summary(self) -> Table:
        """Create a summary table for xTB settings.

        Returns:
            Rich Table with current xTB settings.
        """
        table = Table(
            show_header=False,
            box=None,
            padding=(0, 2),
            expand=True,
        )
        table.add_column("Setting", style="bold")
        table.add_column("Value", style=COLORS["settings"])

        table.add_row("xTB Pre-optimization", self._status_icon(self.xtb.preopt))
        table.add_row("Timeout", f"{self.xtb.timeout_seconds}s")

        return table

    def show_help(self, section: str = "ALL") -> None:
        """Display help text for settings.

        Args:
            section: Section to show help for (CREST, MOPAC, XTB, or ALL).
        """
        self.console.print()

        if section in ("CREST", "ALL"):
            self.console.print(
                Panel(
                    "\n".join(
                        f"[bold]{k}[/bold]: {v}"
                        for k, v in self.HELP_TEXT.items()
                        if k.startswith("crest_")
                    ),
                    title="[bold]CREST Help[/bold]",
                    border_style=COLORS["settings"],
                )
            )

        if section in ("MOPAC", "ALL"):
            self.console.print(
                Panel(
                    "\n".join(
                        f"[bold]{k}[/bold]: {v}"
                        for k, v in self.HELP_TEXT.items()
                        if k.startswith("mopac_")
                    ),
                    title="[bold]MOPAC Help[/bold]",
                    border_style=COLORS["settings"],
                )
            )

        if section in ("XTB", "ALL"):
            self.console.print(
                Panel(
                    "\n".join(
                        f"[bold]{k}[/bold]: {v}"
                        for k, v in self.HELP_TEXT.items()
                        if k.startswith("xtb_")
                    ),
                    title="[bold]xTB Help[/bold]",
                    border_style=COLORS["settings"],
                )
            )

        self.console.print()

    def display_crest_menu(self) -> bool:
        """Interactive CREST settings menu.

        Returns:
            True if settings were saved, False if cancelled.
        """
        while True:
            self.console.clear()
            self.console.print()
            self.console.print(
                Panel(
                    self.show_crest_summary(),
                    title="[bold]CREST Configuration[/bold]",
                    subtitle="[dim]Conformer Search Settings[/dim]",
                    border_style=COLORS["settings"],
                )
            )

            # Format toggle labels with current state
            v3_state = "✓ ON" if self.crest.v3 else "○ OFF"
            nci_state = "✓ ON" if self.crest.nci else "○ OFF"
            xtb_state = "✓ ON" if self.xtb.preopt else "○ OFF"

            choices = [
                questionary.Choice(f"Toggle v3 Algorithm [{v3_state}]", value="v3"),
                questionary.Choice(f"Toggle NCI Mode [{nci_state}]", value="nci"),
                questionary.Choice(
                    f"Set CREST Method (current: {self.crest.crest_method})",
                    value="crest_method",
                ),
                questionary.Choice(
                    f"Set Quick Mode (current: {self.crest.quick_mode})",
                    value="quick_mode",
                ),
                questionary.Choice(
                    f"Set Energy Window (current: {self.crest.ewin})", value="ewin"
                ),
                questionary.Choice(
                    f"Set RMSD Threshold (current: {self.crest.rthr})", value="rthr"
                ),
                questionary.Choice(
                    f"Set Optimization Level (current: {self.crest.optlev})",
                    value="optlev",
                ),
                questionary.Choice(
                    f"Set Threads (current: {self.crest.threads})", value="threads"
                ),
                questionary.Separator(),
                questionary.Choice(
                    f"⚡ xTB Pre-optimization [{xtb_state}]", value="xtb_preopt"
                ),
                questionary.Separator(),
                questionary.Choice("❓ Help", value="help"),
                questionary.Choice("🔄 Reset to Defaults", value="reset"),
                questionary.Separator(),
                questionary.Choice("💾 Save & Return", value="save"),
                questionary.Choice("◀ Cancel", value="cancel"),
            ]

            choice = questionary.select(
                "Select option:",
                choices=choices,
                style=SETTINGS_STYLE,
                qmark="",
                pointer="❯",
            ).ask()

            if choice is None or choice == "cancel":
                return False
            if choice == "save":
                return True
            if choice == "help":
                self.show_help("CREST")
                self.console.input("[dim]Press Enter to continue...[/dim]")
            elif choice == "reset":
                self.reset_crest()
            elif choice == "v3":
                self.crest.v3 = not self.crest.v3
                state = "ON" if self.crest.v3 else "OFF"
                self.console.print(
                    f"[green]✓ Setting updated: v3 Algorithm {state}[/green]"
                )
            elif choice == "crest_method":
                method = questionary.select(
                    "Choose CREST method:",
                    choices=CRESTSettings.CREST_METHOD_OPTIONS,
                    default=self.crest.crest_method,
                    style=SETTINGS_STYLE,
                ).ask()
                if method:
                    self.crest.crest_method = method
                    self.console.print(
                        f"[green]✓ CREST Method set to: {method}[/green]"
                    )
            elif choice == "quick_mode":
                mode = questionary.select(
                    "Choose quick mode:",
                    choices=CRESTSettings.QUICK_MODE_OPTIONS,
                    default=self.crest.quick_mode,
                    style=SETTINGS_STYLE,
                ).ask()
                if mode:
                    self.crest.quick_mode = mode
                    self.console.print(f"[green]✓ Quick Mode set to: {mode}[/green]")
            elif choice == "nci":
                self.crest.nci = not self.crest.nci
                state = "ON" if self.crest.nci else "OFF"
                self.console.print(
                    f"[green]✓ Setting updated: NCI Mode {state}[/green]"
                )
            elif choice == "ewin":
                val = questionary.text(
                    "Energy window (kcal/mol):",
                    default=str(self.crest.ewin),
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        self.crest.ewin = float(val)
                        self.console.print(
                            f"[green]✓ Energy window set to {self.crest.ewin} kcal/mol[/green]"
                        )
                    except ValueError:
                        self.console.print("[red]❌ Invalid number[/red]")
            elif choice == "rthr":
                val = questionary.text(
                    "RMSD threshold (Å):",
                    default=str(self.crest.rthr),
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        self.crest.rthr = float(val)
                        self.console.print(
                            f"[green]✓ RMSD threshold set to {self.crest.rthr} Å[/green]"
                        )
                    except ValueError:
                        self.console.print("[red]❌ Invalid number[/red]")
            elif choice == "optlev":
                opt = questionary.select(
                    "Optimization level:",
                    choices=CRESTSettings.OPTLEV_CHOICES,
                    default=self.crest.optlev,
                    style=SETTINGS_STYLE,
                ).ask()
                if opt:
                    self.crest.optlev = opt
                    self.console.print(
                        f"[green]✓ Optimization level set to {opt}[/green]"
                    )
            elif choice == "threads":
                val = questionary.text(
                    "Number of threads:",
                    default=str(self.crest.threads),
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        threads_val = int(val)
                        max_threads = os.cpu_count() or 4
                        if threads_val < 1:
                            self.console.print("[red]Threads must be >= 1[/red]")
                        elif threads_val > max_threads:
                            self.console.print(
                                f"[red]Threads exceeds system max ({max_threads})[/red]"
                            )
                        else:
                            self.crest.threads = threads_val
                            self.console.print(
                                f"[green]✓ Threads set to {threads_val}[/green]"
                            )
                    except ValueError:
                        self.console.print("[red]Invalid number[/red]")
            elif choice == "xtb_preopt":
                self.xtb.preopt = not self.xtb.preopt
                state = "ON" if self.xtb.preopt else "OFF"
                self.console.print(
                    f"[green]✓ Setting updated: xTB Pre-optimization {state}[/green]"
                )

    def display_mopac_menu(self) -> bool:
        """Interactive MOPAC settings menu.

        Returns:
            True if settings were saved, False if cancelled.
        """
        while True:
            self.console.clear()
            self.console.print()
            self.console.print(
                Panel(
                    self.show_mopac_summary(),
                    title="[bold]MOPAC Configuration[/bold]",
                    subtitle="[dim]PM7 Optimization Settings[/dim]",
                    border_style=COLORS["settings"],
                )
            )

            # Format toggle labels with current state
            precise_state = "✓ ON" if self.mopac.precise else "○ OFF"
            pulay_state = "✓ ON" if self.mopac.pulay else "○ OFF"
            prtall_state = "✓ ON" if self.mopac.prtall else "○ OFF"
            archive_state = "✓ ON" if self.mopac.archive else "○ OFF"

            choices = [
                questionary.Choice(
                    f"Toggle Precise SCF [{precise_state}]", value="precise"
                ),
                questionary.Choice(
                    f"Set SCF Threshold (current: {self.mopac.scfcrt:.1e})",
                    value="scfcrt",
                ),
                questionary.Choice(
                    f"Set Max Iterations (current: {self.mopac.itry})", value="itry"
                ),
                questionary.Choice(
                    f"Toggle PULAY Acceleration [{pulay_state}]", value="pulay"
                ),
                questionary.Choice(
                    f"Toggle Verbose Output [{prtall_state}]", value="prtall"
                ),
                questionary.Choice(
                    f"Toggle Archive Output [{archive_state}]", value="archive"
                ),
                questionary.Separator(),
                questionary.Choice("❓ Help", value="help"),
                questionary.Choice("🔄 Reset to Defaults", value="reset"),
                questionary.Separator(),
                questionary.Choice("💾 Save & Return", value="save"),
                questionary.Choice("◀ Cancel", value="cancel"),
            ]

            choice = questionary.select(
                "Select option:",
                choices=choices,
                style=SETTINGS_STYLE,
                qmark="",
                pointer="❯",
            ).ask()

            if choice is None or choice == "cancel":
                return False
            if choice == "save":
                return True
            if choice == "help":
                self.show_help("MOPAC")
                self.console.input("[dim]Press Enter to continue...[/dim]")
            elif choice == "reset":
                self.reset_mopac()
            elif choice == "precise":
                self.mopac.precise = not self.mopac.precise
                state = "ON" if self.mopac.precise else "OFF"
                self.console.print(
                    f"[green]✓ Setting updated: Precise SCF {state}[/green]"
                )
            elif choice == "pulay":
                self.mopac.pulay = not self.mopac.pulay
                state = "ON" if self.mopac.pulay else "OFF"
                self.console.print(
                    f"[green]✓ Setting updated: PULAY Acceleration {state}[/green]"
                )
            elif choice == "prtall":
                self.mopac.prtall = not self.mopac.prtall
                state = "ON" if self.mopac.prtall else "OFF"
                self.console.print(
                    f"[green]✓ Setting updated: Verbose Output {state}[/green]"
                )
            elif choice == "archive":
                self.mopac.archive = not self.mopac.archive
                state = "ON" if self.mopac.archive else "OFF"
                self.console.print(
                    f"[green]✓ Setting updated: Archive Output {state}[/green]"
                )
            elif choice == "scfcrt":
                val = questionary.text(
                    "SCF threshold (e.g., 1e-4):",
                    default=f"{self.mopac.scfcrt:.1e}",
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        self.mopac.scfcrt = float(val)
                        self.console.print(
                            f"[green]✓ SCF threshold set to {self.mopac.scfcrt:.1e}[/green]"
                        )
                    except ValueError:
                        self.console.print("[red]❌ Invalid number[/red]")
            elif choice == "itry":
                val = questionary.text(
                    "Max iterations:",
                    default=str(self.mopac.itry),
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        self.mopac.itry = int(val)
                        self.console.print(
                            f"[green]✓ Max iterations set to {self.mopac.itry}[/green]"
                        )
                    except ValueError:
                        self.console.print("[red]❌ Invalid number[/red]")

    def display_xtb_menu(self) -> bool:
        """Interactive xTB pre-optimization menu.

        Returns:
            True if settings were saved, False if cancelled.
        """
        while True:
            self.console.clear()
            self.console.print()
            self.console.print(
                Panel(
                    self.show_xtb_summary(),
                    title="[bold]xTB Configuration[/bold]",
                    subtitle="[dim]Pre-optimization Settings[/dim]",
                    border_style=COLORS["settings"],
                )
            )

            choices = [
                questionary.Choice("Toggle xTB Pre-optimization", value="toggle"),
                questionary.Separator(),
                questionary.Choice("❓ Help", value="help"),
                questionary.Separator(),
                questionary.Choice("💾 Save & Return", value="save"),
                questionary.Choice("◀ Cancel", value="cancel"),
            ]

            choice = questionary.select(
                "Select option:",
                choices=choices,
                style=SETTINGS_STYLE,
                qmark="",
                pointer="❯",
            ).ask()

            if choice is None or choice == "cancel":
                return False
            if choice == "save":
                return True
            if choice == "help":
                self.show_help("xTB")
                self.console.input("[dim]Press Enter to continue...[/dim]")
            elif choice == "toggle":
                self.xtb.preopt = not self.xtb.preopt
