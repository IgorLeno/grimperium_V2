"""
Settings Manager for GRIMPERIUM CLI.

Manages CREST, MOPAC, and xTB configuration with interactive menus.
"""

from dataclasses import dataclass, field
from typing import Any, ClassVar

import questionary
from prompt_toolkit.styles import Style
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from grimperium.cli.styles import COLORS


@dataclass
class CRESTSettings:
    """CREST configuration options for conformer search.

    Attributes:
        v3: Use iMTD-GC v3 algorithm (recommended for better conformer sampling).
        quick: Fast mode - reduces search time but may miss structures.
        nci: NCI mode for molecular complexes and weak interactions.
        gfnff: Use GFN-FF force field (default: GFN2-xTB).
        ewin: Energy window in kcal/mol for conformer selection.
        rthr: RMSD threshold in Angstroms for geometry deduplication.
        optlev: Optimization level (loose/normal/tight/vtight/extreme).
        threads: Number of parallel threads.
    """

    v3: bool = True
    quick: bool = False
    nci: bool = False
    gfnff: bool = False
    ewin: float = 5.0
    rthr: float = 0.125
    optlev: str = "normal"
    threads: int = 4

    OPTLEV_CHOICES: ClassVar[list[str]] = [
        "loose",
        "normal",
        "tight",
        "vtight",
        "extreme",
    ]


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

    preopt: bool = False
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
        "crest_quick": "Fast mode - reduces search time but may miss structures",
        "crest_nci": "NCI mode for molecular complexes and weak interactions",
        "crest_gfnff": "Use GFN-FF force field (default: GFN2-xTB)",
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

    def to_dict(self) -> dict[str, Any]:
        """Convert all settings to a flat dictionary.

        Returns:
            Dictionary with prefixed keys for all settings.
        """
        return {
            "crest_v3": self.crest.v3,
            "crest_quick": self.crest.quick,
            "crest_nci": self.crest.nci,
            "crest_gfnff": self.crest.gfnff,
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
        if "crest_v3" in settings_dict:
            self.crest.v3 = bool(settings_dict["crest_v3"])
        if "crest_quick" in settings_dict:
            self.crest.quick = bool(settings_dict["crest_quick"])
        if "crest_nci" in settings_dict:
            self.crest.nci = bool(settings_dict["crest_nci"])
        if "crest_gfnff" in settings_dict:
            self.crest.gfnff = bool(settings_dict["crest_gfnff"])
        if "crest_ewin" in settings_dict:
            self.crest.ewin = float(settings_dict["crest_ewin"])
        if "crest_rthr" in settings_dict:
            self.crest.rthr = float(settings_dict["crest_rthr"])
        if "crest_optlev" in settings_dict:
            val = str(settings_dict["crest_optlev"])
            if val in CRESTSettings.OPTLEV_CHOICES:
                self.crest.optlev = val
        if "crest_threads" in settings_dict:
            self.crest.threads = int(settings_dict["crest_threads"])
        if "mopac_precise" in settings_dict:
            self.mopac.precise = bool(settings_dict["mopac_precise"])
        if "mopac_scfcrt" in settings_dict:
            self.mopac.scfcrt = float(settings_dict["mopac_scfcrt"])
        if "mopac_itry" in settings_dict:
            self.mopac.itry = int(settings_dict["mopac_itry"])
        if "mopac_pulay" in settings_dict:
            self.mopac.pulay = bool(settings_dict["mopac_pulay"])
        if "mopac_prtall" in settings_dict:
            self.mopac.prtall = bool(settings_dict["mopac_prtall"])
        if "mopac_archive" in settings_dict:
            self.mopac.archive = bool(settings_dict["mopac_archive"])
        if "crest_xtb_preopt" in settings_dict:
            self.xtb.preopt = bool(settings_dict["crest_xtb_preopt"])

    def reset_crest(self) -> None:
        """Reset CREST settings to defaults."""
        self.crest = CRESTSettings()

    def reset_mopac(self) -> None:
        """Reset MOPAC settings to defaults."""
        self.mopac = MOPACSettings()

    def reset_xtb(self) -> None:
        """Reset xTB settings to defaults."""
        self.xtb = xTBSettings()

    def reset_all(self) -> None:
        """Reset all settings to defaults."""
        self.reset_crest()
        self.reset_mopac()
        self.reset_xtb()

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
        table.add_row("Quick Mode", self._status_icon(self.crest.quick))
        table.add_row("NCI Mode", self._status_icon(self.crest.nci))
        table.add_row("GFN-FF Force Field", self._status_icon(self.crest.gfnff))
        table.add_row("Energy Window", f"{self.crest.ewin} kcal/mol")
        table.add_row("RMSD Threshold", f"{self.crest.rthr} Å")
        table.add_row("Optimization Level", self.crest.optlev)
        table.add_row("Threads", str(self.crest.threads))

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

            choices = [
                questionary.Choice("Toggle v3 Algorithm", value="v3"),
                questionary.Choice("Toggle Quick Mode", value="quick"),
                questionary.Choice("Toggle NCI Mode", value="nci"),
                questionary.Choice("Toggle GFN-FF", value="gfnff"),
                questionary.Choice("Set Energy Window", value="ewin"),
                questionary.Choice("Set RMSD Threshold", value="rthr"),
                questionary.Choice("Set Optimization Level", value="optlev"),
                questionary.Choice("Set Threads", value="threads"),
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
            elif choice == "quick":
                self.crest.quick = not self.crest.quick
            elif choice == "nci":
                self.crest.nci = not self.crest.nci
            elif choice == "gfnff":
                self.crest.gfnff = not self.crest.gfnff
            elif choice == "ewin":
                val = questionary.text(
                    "Energy window (kcal/mol):",
                    default=str(self.crest.ewin),
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        self.crest.ewin = float(val)
                    except ValueError:
                        self.console.print("[red]Invalid number[/red]")
            elif choice == "rthr":
                val = questionary.text(
                    "RMSD threshold (Å):",
                    default=str(self.crest.rthr),
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        self.crest.rthr = float(val)
                    except ValueError:
                        self.console.print("[red]Invalid number[/red]")
            elif choice == "optlev":
                opt = questionary.select(
                    "Optimization level:",
                    choices=CRESTSettings.OPTLEV_CHOICES,
                    default=self.crest.optlev,
                    style=SETTINGS_STYLE,
                ).ask()
                if opt:
                    self.crest.optlev = opt
            elif choice == "threads":
                val = questionary.text(
                    "Number of threads:",
                    default=str(self.crest.threads),
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        self.crest.threads = int(val)
                    except ValueError:
                        self.console.print("[red]Invalid number[/red]")

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

            choices = [
                questionary.Choice("Toggle Precise SCF", value="precise"),
                questionary.Choice("Set SCF Threshold", value="scfcrt"),
                questionary.Choice("Set Max Iterations", value="itry"),
                questionary.Choice("Toggle PULAY Acceleration", value="pulay"),
                questionary.Choice("Toggle Verbose Output", value="prtall"),
                questionary.Choice("Toggle Archive Output", value="archive"),
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
            elif choice == "pulay":
                self.mopac.pulay = not self.mopac.pulay
            elif choice == "prtall":
                self.mopac.prtall = not self.mopac.prtall
            elif choice == "archive":
                self.mopac.archive = not self.mopac.archive
            elif choice == "scfcrt":
                val = questionary.text(
                    "SCF threshold (e.g., 1e-4):",
                    default=f"{self.mopac.scfcrt:.1e}",
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        self.mopac.scfcrt = float(val)
                    except ValueError:
                        self.console.print("[red]Invalid number[/red]")
            elif choice == "itry":
                val = questionary.text(
                    "Max iterations:",
                    default=str(self.mopac.itry),
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    try:
                        self.mopac.itry = int(val)
                    except ValueError:
                        self.console.print("[red]Invalid number[/red]")

    def display_xtb_menu(self) -> bool:
        """Interactive xTB pre-optimization menu.

        Returns:
            True if settings were saved, False if cancelled.
        """
        self.console.clear()
        self.console.print()
        self.console.print(
            Panel(
                self.show_xtb_summary(),
                title="[bold]xTB Pre-optimization[/bold]",
                subtitle="[dim]Pre-optimize structures before CREST[/dim]",
                border_style=COLORS["settings"],
            )
        )

        self.console.print()
        self.console.print("[dim]Recommended: Yes (per CREST documentation)[/dim]")
        self.console.print("[dim]Cost: ~15 seconds per molecule (optional)[/dim]")
        self.console.print()

        choices = [
            questionary.Choice(
                "Enable xTB Pre-optimization",
                value="enable",
            ),
            questionary.Choice(
                "Disable xTB Pre-optimization",
                value="disable",
            ),
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
            self.console.print()
            self.console.print(
                Panel(
                    self.HELP_TEXT["xtb_preopt"],
                    title="[bold]xTB Help[/bold]",
                    border_style=COLORS["settings"],
                )
            )
            self.console.input("[dim]Press Enter to continue...[/dim]")
            return self.display_xtb_menu()
        if choice == "enable":
            self.xtb.preopt = True
            return True
        if choice == "disable":
            self.xtb.preopt = False
            return True

        return False
