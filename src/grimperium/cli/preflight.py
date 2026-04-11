"""
Preflight environment check for GRIMPERIUM CLI.

Verifies Python libraries, external binaries, and configuration
before the main menu is displayed. Results are cached to avoid
re-checking on every launch.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import logging
import platform
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rich.console import Console

logger = logging.getLogger(__name__)

# ── Import-name to pip-package mapping (only where they differ) ───────────
_IMPORT_TO_PACKAGE: dict[str, str] = {
    "sklearn": "scikit-learn",
}

# Ordered list of (import_name, pip_package_name)
REQUIRED_PYTHON_LIBS: list[tuple[str, str]] = [
    ("pandas", "pandas"),
    ("numpy", "numpy"),
    ("sklearn", "scikit-learn"),
    ("xgboost", "xgboost"),
    ("rdkit", "rdkit"),
    ("rich", "rich"),
    ("joblib", "joblib"),
]

# Binary checks: (executable_name, display_name, version_flag, install_hint, critical)
REQUIRED_BINARIES: list[dict[str, str | bool]] = [
    {
        "name": "crest",
        "display_name": "CREST",
        "version_flag": "--version",
        "install_hint": (
            "Install CREST: https://crest-lab.github.io/crest-docs/page/installation"
        ),
        "critical": True,
    },
    {
        "name": "mopac",
        "display_name": "MOPAC",
        "version_flag": "-v",
        "install_hint": ("Install MOPAC: https://github.com/openmopac/mopac/releases"),
        "critical": True,
    },
    {
        "name": "xtb",
        "display_name": "xTB",
        "version_flag": "--version",
        "install_hint": (
            "Install xTB: https://github.com/grimme-lab/xtb/releases\n"
            "Note: Many CREST distributions bundle xTB internally."
        ),
        "critical": False,
    },
    {
        "name": "obabel",
        "display_name": "Open Babel",
        "version_flag": "--version",
        "install_hint": (
            "Install Open Babel: https://openbabel.org/docs/Installation\n"
            "Note: RDKit fallback is available for XYZ->SDF conversion."
        ),
        "critical": False,
    },
]


# ── Data classes ──────────────────────────────────────────────────────────


@dataclass
class CheckResult:
    """Result of a single preflight check."""

    name: str
    category: str  # "python_lib" | "binary" | "binary_optional" | "config"
    status: str  # "ok" | "missing" | "wrong_version" | "not_configured"
    detail: str
    fix_available: bool = False
    fix_command: str = ""
    install_hint: str = ""

    @property
    def is_ok(self) -> bool:
        """Check passed."""
        return self.status == "ok"

    @property
    def is_critical(self) -> bool:
        """Check is a critical failure (blocks full functionality)."""
        return not self.is_ok and self.category not in ("config", "binary_optional")


@dataclass
class PreflightCache:
    """Manages the preflight result cache at ~/.grimperium/preflight_cache.json."""

    cache_dir: Path = field(default_factory=lambda: Path.home() / ".grimperium")

    @property
    def cache_path(self) -> Path:
        return self.cache_dir / "preflight_cache.json"

    def load(self) -> dict[str, object] | None:
        """Load and return cache data, or None if invalid/missing."""
        if not self.cache_path.is_file():
            return None
        try:
            data = json.loads(self.cache_path.read_text())
            if not isinstance(data, dict):
                return None
            return data
        except (json.JSONDecodeError, OSError):
            return None

    def is_valid(self, current_version: str) -> bool:
        """Check if cached results are still valid.

        Invalidation triggers:
        - Cache file missing or corrupted
        - grimperium_version changed
        - python_version changed (major.minor.patch)
        - settings_hash changed
        """
        data = self.load()
        if data is None:
            return False

        if data.get("grimperium_version") != current_version:
            return False

        if data.get("python_version") != platform.python_version():
            return False

        # Invalidate if settings file hash changed
        current_hash = _compute_settings_hash()
        cached_hash = data.get("settings_hash", "")
        return current_hash == cached_hash

    def save(
        self,
        checks_passed: list[str],
        grimperium_version: str,
    ) -> None:
        """Write valid cache after all checks pass."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        payload = {
            "grimperium_version": grimperium_version,
            "python_version": platform.python_version(),
            "last_ok_timestamp": datetime.now(timezone.utc).isoformat(),
            "checks_passed": checks_passed,
            "settings_hash": _compute_settings_hash(),
        }
        self.cache_path.write_text(json.dumps(payload, indent=2))

    def invalidate(self) -> None:
        """Delete the cache file."""
        try:
            self.cache_path.unlink()
        except FileNotFoundError:
            pass


# ── Utility helpers ───────────────────────────────────────────────────────


def _compute_settings_hash() -> str:
    """Compute MD5 hash of the settings file, or empty string if absent."""
    from grimperium.cli.settings_manager import SettingsManager

    path = SettingsManager.default_config_path()
    if not path.is_file():
        return ""
    try:
        return hashlib.md5(path.read_bytes()).hexdigest()
    except OSError:
        return ""


def _safe_confirm(message: str, default: bool = False) -> bool:
    """confirm() wrapper that returns *default* when stdin is not a TTY."""
    if not sys.stdin.isatty():
        logger.info(
            "Non-interactive session, using default=%s for: %s", default, message
        )
        return default
    from grimperium.cli.menu import confirm

    return confirm(message, default=default)


def _get_version() -> str:
    """Return the current grimperium version string."""
    from grimperium import __version__

    return __version__


# ── Individual check functions ────────────────────────────────────────────


def check_python_lib(import_name: str, pip_name: str) -> CheckResult:
    """Check if a Python library is importable.

    Args:
        import_name: Module name to import (e.g. "sklearn").
        pip_name: Package name for pip install (e.g. "scikit-learn").
    """
    try:
        mod = importlib.import_module(import_name)
        version = getattr(mod, "__version__", "installed")
        return CheckResult(
            name=pip_name,
            category="python_lib",
            status="ok",
            detail=f"{pip_name} {version}",
        )
    except ImportError:
        hint = f"pip install {pip_name}"
        if pip_name == "rdkit":
            hint += "  (or: conda install -c conda-forge rdkit)"
        return CheckResult(
            name=pip_name,
            category="python_lib",
            status="missing",
            detail=f"{pip_name} not installed",
            fix_available=True,
            fix_command=f"{sys.executable} -m pip install {pip_name}",
            install_hint=hint,
        )


def check_binary(
    name: str,
    display_name: str,
    version_flag: str,
    install_hint: str,
    critical: bool,
) -> CheckResult:
    """Check if an external binary is available.

    Reuses check_executable from grimperium.crest_pm7.validation.
    """
    from grimperium.crest_pm7.validation import check_executable

    found, version = check_executable(name, version_flag)
    category = "binary" if critical else "binary_optional"
    if found:
        return CheckResult(
            name=display_name,
            category=category,
            status="ok",
            detail=f"{display_name}: {version or 'found'}",
        )
    return CheckResult(
        name=display_name,
        category=category,
        status="missing",
        detail=f"{display_name} not found in PATH",
        fix_available=False,
        install_hint=install_hint,
    )


def check_settings_file() -> CheckResult:
    """Check if grimperium_settings.json exists."""
    from grimperium.cli.settings_manager import SettingsManager

    path = SettingsManager.default_config_path()
    if path.is_file():
        return CheckResult(
            name="grimperium_settings.json",
            category="config",
            status="ok",
            detail=f"Found at {path}",
        )
    return CheckResult(
        name="grimperium_settings.json",
        category="config",
        status="not_configured",
        detail="Settings file not found (create via Settings menu)",
        fix_available=False,
        install_hint="Navigate to Settings > Save as Default to create this file",
    )


# ── Status icons ──────────────────────────────────────────────────────────

_STATUS_ICON: dict[str, str] = {
    "ok": "[success]✓[/success]",
    "missing": "[error]✗[/error]",
    "wrong_version": "[warning]⚠[/warning]",
    "not_configured": "[warning]⚠[/warning]",
}


# ── Orchestrator ──────────────────────────────────────────────────────────


class PreflightRunner:
    """Orchestrates the preflight check sequence with Rich display."""

    def __init__(self, console: Console) -> None:
        self.console = console
        self.cache = PreflightCache()
        self.results: list[CheckResult] = []

    def run(self, force: bool = False) -> bool:
        """Execute preflight checks.

        Args:
            force: If True, ignore cache and run all checks.

        Returns:
            True if the application should proceed, False to exit.
        """
        version = _get_version()

        # ── Cache fast-path ───────────────────────────────────────────
        if not force and self.cache.is_valid(version):
            self.console.print(
                "[success]✓ Environment OK[/success] [dim](cached)[/dim]"
            )
            self.console.print()
            return True

        # ── Run checks ────────────────────────────────────────────────
        from rich.panel import Panel

        from grimperium.cli.styles import COLORS

        self.console.print()
        self.console.print(
            Panel(
                "[bold]SYSTEM PREFLIGHT CHECK[/bold]",
                border_style=COLORS["calc"],
                padding=(0, 2),
            )
        )
        self.console.print()

        self.results = self._run_all_checks()
        self._display_results_table(self.results)

        # ── Evaluate results ──────────────────────────────────────────
        if all(r.is_ok for r in self.results):
            self._save_cache(version)
            self.console.print()
            self.console.print("[success]✓ Environment ready[/success]")
            self.console.print()
            return True

        # ── Handle failures ───────────────────────────────────────────
        return self._handle_failures(version)

    def _run_all_checks(self) -> list[CheckResult]:
        """Execute all check functions and return results."""
        results: list[CheckResult] = []

        # Python libraries
        for import_name, pip_name in REQUIRED_PYTHON_LIBS:
            results.append(check_python_lib(import_name, pip_name))

        # External binaries
        for binary in REQUIRED_BINARIES:
            results.append(
                check_binary(
                    name=str(binary["name"]),
                    display_name=str(binary["display_name"]),
                    version_flag=str(binary["version_flag"]),
                    install_hint=str(binary["install_hint"]),
                    critical=bool(binary["critical"]),
                )
            )

        # Configuration
        results.append(check_settings_file())

        return results

    def _display_results_table(self, results: list[CheckResult]) -> None:
        """Display check results as a Rich table."""
        from rich.table import Table

        from grimperium.cli.styles import COLORS

        table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
        table.add_column("Check", style="bold", min_width=25)
        table.add_column("Status", min_width=6, justify="center")
        table.add_column("Detail", style=COLORS["muted"])

        for r in results:
            icon = _STATUS_ICON.get(r.status, "?")
            table.add_row(r.name, icon, r.detail)

        self.console.print(table)

    def _handle_failures(self, version: str) -> bool:
        """Handle check failures: offer fixes, re-check, ask to proceed."""
        self.console.print()

        # ── Offer pip install for missing Python libs ─────────────────
        fixable = [r for r in self.results if r.fix_available and not r.is_ok]
        if fixable:
            pkg_names = [r.name for r in fixable]
            self.console.print(
                f"[warning]Missing Python packages: {', '.join(pkg_names)}[/warning]"
            )
            if _safe_confirm(
                f"Install {len(pkg_names)} missing package(s) via pip?",
                default=False,
            ):
                success = self._attempt_pip_install(pkg_names)
                if success:
                    self.console.print(
                        "[success]✓ Packages installed successfully[/success]"
                    )
                else:
                    self.console.print(
                        "[error]✗ Some packages failed to install[/error]"
                    )

        # ── Show hints for missing binaries ───────────────────────────
        failed_binaries = [
            r
            for r in self.results
            if r.category in ("binary", "binary_optional") and not r.is_ok
        ]
        if failed_binaries:
            self._show_binary_hints(failed_binaries)

        # ── Show config hints ─────────────────────────────────────────
        config_issues = [
            r for r in self.results if r.category == "config" and not r.is_ok
        ]
        for r in config_issues:
            self.console.print(f"[warning]⚠ {r.name}: {r.install_hint}[/warning]")

        # ── Re-run failed checks ──────────────────────────────────────
        self.console.print()
        self.console.print("[muted]Re-checking failed items...[/muted]")
        recheck_results = self._rerun_failed()
        self._display_results_table(recheck_results)

        # Merge recheck into results
        recheck_by_name = {r.name: r for r in recheck_results}
        self.results = [recheck_by_name.get(r.name, r) for r in self.results]

        # ── Final evaluation ──────────────────────────────────────────
        if all(r.is_ok for r in self.results):
            self._save_cache(version)
            self.console.print()
            self.console.print("[success]✓ All issues resolved[/success]")
            self.console.print()
            return True

        critical = [r for r in self.results if r.is_critical]
        if critical:
            self.console.print()
            names = ", ".join(r.name for r in critical)
            self.console.print(f"[error]Critical: {names} still missing[/error]")
            if _safe_confirm("Start anyway with limited functionality?", default=False):
                logger.warning(
                    "User chose to start with missing dependencies: %s", names
                )
                return True
            return False

        # Only non-critical issues remain (config, optional binaries)
        self._save_cache(version)
        self.console.print()
        self.console.print("[warning]⚠ Starting with warnings (non-critical)[/warning]")
        self.console.print()
        return True

    def _rerun_failed(self) -> list[CheckResult]:
        """Re-run only previously failed checks."""
        recheck: list[CheckResult] = []
        for r in self.results:
            if r.is_ok:
                continue
            if r.category == "python_lib":
                # Find original import name
                pip_to_import = {p: i for i, p in REQUIRED_PYTHON_LIBS}
                import_name = pip_to_import.get(r.name, r.name)
                recheck.append(check_python_lib(import_name, r.name))
            elif r.category in ("binary", "binary_optional"):
                # Find original binary spec
                for b in REQUIRED_BINARIES:
                    if b["display_name"] == r.name:
                        recheck.append(
                            check_binary(
                                name=str(b["name"]),
                                display_name=str(b["display_name"]),
                                version_flag=str(b["version_flag"]),
                                install_hint=str(b["install_hint"]),
                                critical=bool(b["critical"]),
                            )
                        )
                        break
            elif r.category == "config":
                recheck.append(check_settings_file())
        return recheck

    def _attempt_pip_install(self, packages: list[str]) -> bool:
        """Run pip install for the given packages.

        Returns True on success, False on failure.
        """
        cmd = [sys.executable, "-m", "pip", "install", *packages]
        self.console.print(f"[muted]Running: {' '.join(cmd)}[/muted]")
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
            )
            if result.returncode != 0:
                logger.error("pip install failed: %s", result.stderr)
                self.console.print(f"[error]{result.stderr.strip()}[/error]")
                return False
            return True
        except subprocess.TimeoutExpired:
            self.console.print("[error]pip install timed out (120s)[/error]")
            return False
        except subprocess.SubprocessError as e:
            self.console.print(f"[error]pip install error: {e}[/error]")
            return False

    def _show_binary_hints(self, failed: list[CheckResult]) -> None:
        """Display installation hints for missing binaries."""
        from rich.panel import Panel

        from grimperium.cli.styles import COLORS

        for r in failed:
            border = COLORS["error"] if r.is_critical else COLORS["warning"]
            level = "REQUIRED" if r.is_critical else "OPTIONAL"
            self.console.print()
            self.console.print(
                Panel(
                    f"[bold]{r.name}[/bold] [{level}]\n\n{r.install_hint}",
                    border_style=border,
                    padding=(1, 2),
                )
            )

    def _save_cache(self, version: str) -> None:
        """Save successful check results to cache."""
        passed = [r.name for r in self.results if r.is_ok]
        self.cache.save(passed, version)


# ── Module-level convenience function ─────────────────────────────────────


def run_preflight(console: Console, force: bool = False) -> bool:
    """Run preflight environment checks.

    Args:
        console: Rich console instance (from controller).
        force: Force re-check ignoring cache.

    Returns:
        True to proceed with CLI startup, False to exit.
    """
    runner = PreflightRunner(console)
    return runner.run(force=force)
