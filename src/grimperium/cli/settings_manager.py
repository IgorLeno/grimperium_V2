"""
Settings Manager for GRIMPERIUM CLI.

Manages CREST, MOPAC, and xTB configuration with interactive menus.
"""

import json
import os
from dataclasses import asdict, dataclass, field, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

from prompt_toolkit.styles import Style
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from grimperium.cli._questionary import questionary
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


@dataclass
class CalculationProfile:
    """Named calculation profile storing CREST and MOPAC parameters.

    Profiles are saved globally in ~/.grimperium/profiles.json and can be
    reused across projects or assigned to remote workers.
    """

    name: str = "default"
    crest_ewin: float = 6.0
    crest_rthr: float = 0.125
    crest_opt: int = 2
    crest_threads: int = 4
    crest_timeout_minutes: int = 60
    mopac_keywords: str = "PM7 EF GNORM=0.01"
    mopac_timeout_minutes: int = 30
    is_standard: bool = False
    created_at: str = ""

    def __post_init__(self) -> None:
        if not self.created_at:
            self.created_at = datetime.now(timezone.utc).isoformat()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "CalculationProfile":
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in known})


@dataclass
class DistributedDefaults:
    """Default parameters used by the server when workers register.

    Persisted in ~/.grimperium/distributed_defaults.json.
    """

    profile_name: str = "default"
    batch_size: int = 10
    crest_timeout_minutes: int = 60
    mopac_timeout_minutes: int = 30

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "DistributedDefaults":
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in data.items() if k in known})


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

    # ── Global distributed config (~/.grimperium/) ─────────────────────

    @staticmethod
    def grimperium_home_dir() -> Path:
        """Return ~/.grimperium/, creating it if necessary."""
        try:
            home = Path.home() / ".grimperium"
        except RuntimeError:
            raise RuntimeError(
                "Cannot determine home directory ($HOME not set). "
                "Set $HOME before running Grimperium."
            ) from None
        home.mkdir(parents=True, exist_ok=True)
        return home

    @staticmethod
    def profiles_path() -> Path:
        """Path to ~/.grimperium/profiles.json."""
        return SettingsManager.grimperium_home_dir() / "profiles.json"

    @staticmethod
    def distributed_defaults_path() -> Path:
        """Path to ~/.grimperium/distributed_defaults.json."""
        return SettingsManager.grimperium_home_dir() / "distributed_defaults.json"

    @staticmethod
    def load_profiles() -> list[CalculationProfile]:
        """Load profiles from disk, auto-creating the default if absent.

        Returns:
            List of CalculationProfile objects; always contains at least
            the immutable 'default' profile.
        """
        path = SettingsManager.profiles_path()
        profiles: list[CalculationProfile] = []
        if path.is_file():
            try:
                raw = json.loads(path.read_text(encoding="utf-8"))
                profiles = [
                    CalculationProfile.from_dict(p) for p in raw.get("profiles", [])
                ]
            except (json.JSONDecodeError, OSError, TypeError, KeyError):
                profiles = []

        if not any(p.name == "default" for p in profiles):
            profiles.insert(0, CalculationProfile(name="default", is_standard=True))

        return profiles

    @staticmethod
    def save_profiles(profiles: list[CalculationProfile]) -> bool:
        """Persist profiles to ~/.grimperium/profiles.json.

        Args:
            profiles: List of profiles to save.

        Returns:
            True on success, False on I/O error.
        """
        path = SettingsManager.profiles_path()
        try:
            path.write_text(
                json.dumps(
                    {"profiles": [p.to_dict() for p in profiles]},
                    indent=2,
                    ensure_ascii=False,
                ),
                encoding="utf-8",
            )
            return True
        except OSError:
            return False

    @staticmethod
    def load_distributed_defaults() -> DistributedDefaults:
        """Load distributed defaults from disk, falling back to dataclass defaults."""
        path = SettingsManager.distributed_defaults_path()
        if path.is_file():
            try:
                raw = json.loads(path.read_text(encoding="utf-8"))
                return DistributedDefaults.from_dict(raw)
            except (json.JSONDecodeError, OSError, TypeError):
                pass
        return DistributedDefaults()

    @staticmethod
    def save_distributed_defaults(defaults: DistributedDefaults) -> bool:
        """Persist distributed defaults to ~/.grimperium/distributed_defaults.json.

        Args:
            defaults: DistributedDefaults instance to save.

        Returns:
            True on success, False on I/O error.
        """
        path = SettingsManager.distributed_defaults_path()
        try:
            path.write_text(
                json.dumps(defaults.to_dict(), indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            return True
        except OSError:
            return False

    def display_readme_updater_menu(self) -> None:
        """Launch the README updater workflow from Settings."""
        from grimperium.cli.readme_updater import ReadmeUpdater

        ReadmeUpdater(console=self.console).display_menu()

    # ── Distributed Settings ──────────────────────────────────────────────────

    def display_distributed_menu(self) -> None:
        """Top-level Distributed Settings menu: Profiles and Standard Values."""
        while True:
            self.console.print()
            self.console.print(
                Panel(
                    "[dim]Manage calculation profiles and default dispatch settings "
                    "for distributed processing.[/dim]",
                    title="[bold]Distributed Settings[/bold]",
                    border_style=COLORS["settings"],
                )
            )
            choices = [
                questionary.Choice("Calculation Profiles", value="profiles"),
                questionary.Choice("Standard Values (Defaults)", value="defaults"),
                questionary.Separator(),
                questionary.Choice("◀ Back", value="back"),
            ]
            choice = questionary.select(
                "Select option:",
                choices=choices,
                style=SETTINGS_STYLE,
                qmark="",
                pointer="❯",
            ).ask()

            if choice is None or choice == "back":
                return
            if choice == "profiles":
                self._display_profiles_menu()
            elif choice == "defaults":
                self._display_standard_values_menu()

    def _display_profiles_menu(self) -> None:
        """List and manage calculation profiles."""
        while True:
            profiles = SettingsManager.load_profiles()
            # Ensure "default" profile always exists
            if not any(p.name == "default" for p in profiles):
                profiles.insert(0, CalculationProfile(name="default", is_standard=True))
                SettingsManager.save_profiles(profiles)

            self.console.print()
            table = Table(
                title="Calculation Profiles",
                show_header=True,
                header_style=f"bold {COLORS['settings']}",
            )
            table.add_column("Name")
            table.add_column("CREST Timeout", justify="right")
            table.add_column("MOPAC Timeout", justify="right")
            table.add_column("Batch Size", justify="right")
            table.add_column("Standard")
            for p in profiles:
                table.add_row(
                    p.name,
                    f"{p.crest_timeout_minutes}m",
                    f"{p.mopac_timeout_minutes}m",
                    "—",
                    "✓" if p.is_standard else "",
                )
            self.console.print(table)

            choices = [
                questionary.Choice("Create new profile", value="create"),
                questionary.Separator(),
                questionary.Choice("◀ Back", value="back"),
            ]
            choice = questionary.select(
                "Select option:",
                choices=choices,
                style=SETTINGS_STYLE,
                qmark="",
                pointer="❯",
            ).ask()

            if choice is None or choice == "back":
                return
            if choice == "create":
                self._create_profile_interactive(profiles)

    def _create_profile_interactive(self, profiles: list[CalculationProfile]) -> None:
        """Interactively create a new calculation profile."""
        name = questionary.text(
            "Profile name (leave blank to cancel):",
            style=SETTINGS_STYLE,
        ).ask()
        if not name or not name.strip():
            return
        name = name.strip()
        if any(p.name == name for p in profiles):
            self.console.print(f"[red]✗ Profile '{name}' already exists.[/red]")
            return

        crest_t_str = questionary.text(
            "CREST timeout (minutes) [default: 60]:",
            default="60",
            style=SETTINGS_STYLE,
        ).ask()
        mopac_t_str = questionary.text(
            "MOPAC timeout (minutes) [default: 30]:",
            default="30",
            style=SETTINGS_STYLE,
        ).ask()

        try:
            crest_t = max(1, int(crest_t_str or "60"))
            mopac_t = max(1, int(mopac_t_str or "30"))
        except ValueError:
            self.console.print(
                "[red]✗ Invalid timeout value — profile not created.[/red]"
            )
            return

        new_profile = CalculationProfile(
            name=name,
            crest_timeout_minutes=crest_t,
            mopac_timeout_minutes=mopac_t,
        )
        profiles.append(new_profile)
        if SettingsManager.save_profiles(profiles):
            self.console.print(f"[green]✓ Profile '{name}' created.[/green]")
        else:
            self.console.print("[red]✗ Failed to save profiles.[/red]")

    def _display_standard_values_menu(self) -> None:
        """Edit DistributedDefaults interactively."""
        defaults = SettingsManager.load_distributed_defaults()

        while True:
            self.console.print()
            self.console.print(
                Panel(
                    f"  Profile:       {defaults.profile_name}\n"
                    f"  Batch size:    {defaults.batch_size}\n"
                    f"  CREST timeout: {defaults.crest_timeout_minutes}m\n"
                    f"  MOPAC timeout: {defaults.mopac_timeout_minutes}m",
                    title="[bold]Standard Values[/bold]",
                    border_style=COLORS["settings"],
                )
            )
            choices = [
                questionary.Choice(
                    f"Set profile (current: {defaults.profile_name})",
                    value="profile",
                ),
                questionary.Choice(
                    f"Set batch size (current: {defaults.batch_size})",
                    value="batch",
                ),
                questionary.Choice(
                    f"Set CREST timeout (current: {defaults.crest_timeout_minutes}m)",
                    value="crest",
                ),
                questionary.Choice(
                    f"Set MOPAC timeout (current: {defaults.mopac_timeout_minutes}m)",
                    value="mopac",
                ),
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
                return
            if choice == "save":
                if SettingsManager.save_distributed_defaults(defaults):
                    self.console.print("[green]✓ Standard values saved.[/green]")
                else:
                    self.console.print("[red]✗ Failed to save.[/red]")
                return
            if choice == "profile":
                profiles = SettingsManager.load_profiles()
                names = [p.name for p in profiles] if profiles else ["default"]
                val = questionary.select(
                    "Choose profile:",
                    choices=names,
                    style=SETTINGS_STYLE,
                ).ask()
                if val:
                    defaults.profile_name = val
            elif choice == "batch":
                raw = questionary.text(
                    f"Batch size [{defaults.batch_size}]:",
                    default=str(defaults.batch_size),
                    style=SETTINGS_STYLE,
                ).ask()
                try:
                    defaults.batch_size = max(1, int(raw or str(defaults.batch_size)))
                except ValueError:
                    self.console.print("[red]✗ Invalid value.[/red]")
            elif choice == "crest":
                raw = questionary.text(
                    f"CREST timeout in minutes [{defaults.crest_timeout_minutes}]:",
                    default=str(defaults.crest_timeout_minutes),
                    style=SETTINGS_STYLE,
                ).ask()
                try:
                    defaults.crest_timeout_minutes = max(
                        1, int(raw or str(defaults.crest_timeout_minutes))
                    )
                except ValueError:
                    self.console.print("[red]✗ Invalid value.[/red]")
            elif choice == "mopac":
                raw = questionary.text(
                    f"MOPAC timeout in minutes [{defaults.mopac_timeout_minutes}]:",
                    default=str(defaults.mopac_timeout_minutes),
                    style=SETTINGS_STYLE,
                ).ask()
                try:
                    defaults.mopac_timeout_minutes = max(
                        1, int(raw or str(defaults.mopac_timeout_minutes))
                    )
                except ValueError:
                    self.console.print("[red]✗ Invalid value.[/red]")
