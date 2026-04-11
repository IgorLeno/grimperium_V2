"""Unit tests for preflight environment check."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from rich.console import Console

from grimperium.cli.preflight import (
    CheckResult,
    PreflightCache,
    PreflightRunner,
    check_binary,
    check_python_lib,
    check_settings_file,
)

# ── Fixtures ──────────────────────────────────────────────────────────────


@pytest.fixture
def console() -> Console:
    """Non-interactive console for testing (with CLI theme for style tags)."""
    from grimperium.cli.styles import CLI_THEME

    return Console(force_terminal=False, no_color=True, theme=CLI_THEME)


@pytest.fixture
def cache_dir(tmp_path: Path) -> Path:
    """Temporary cache directory."""
    return tmp_path / ".grimperium"


@pytest.fixture
def cache(cache_dir: Path) -> PreflightCache:
    """PreflightCache pointing at tmp dir."""
    return PreflightCache(cache_dir=cache_dir)


@pytest.fixture
def runner(console: Console, cache: PreflightCache) -> PreflightRunner:
    """PreflightRunner with mocked cache location."""
    r = PreflightRunner(console)
    r.cache = cache
    return r


# ── Helpers ───────────────────────────────────────────────────────────────

_ALL_OK = [
    CheckResult("pandas", "python_lib", "ok", "pandas 2.0"),
    CheckResult("numpy", "python_lib", "ok", "numpy 1.26"),
    CheckResult("CREST", "binary", "ok", "CREST found"),
    CheckResult("MOPAC", "binary", "ok", "MOPAC found"),
    CheckResult("xTB", "binary_optional", "ok", "xTB found"),
    CheckResult("Open Babel", "binary_optional", "ok", "obabel found"),
    CheckResult("grimperium_settings.json", "config", "ok", "found"),
]

_VERSION = "0.2.0"


# ── Required Test 1 ──────────────────────────────────────────────────────


class TestCacheValidSkipsChecks:
    """When cache is valid, no checks are executed."""

    def test_cache_valid_skips_checks(
        self, runner: PreflightRunner, cache: PreflightCache
    ) -> None:
        with patch(
            "grimperium.cli.preflight._compute_settings_hash",
            return_value="",
        ):
            cache.save(
                checks_passed=["pandas", "numpy", "CREST"],
                grimperium_version=_VERSION,
            )

        with (
            patch.object(runner, "_run_all_checks") as mock_checks,
            patch("grimperium.cli.preflight._get_version", return_value=_VERSION),
            patch(
                "grimperium.cli.preflight._compute_settings_hash",
                return_value="",
            ),
        ):
            result = runner.run(force=False)

        assert result is True
        mock_checks.assert_not_called()


# ── Required Test 2 ──────────────────────────────────────────────────────


class TestCacheInvalidOnVersionChange:
    """Cache is invalid when grimperium version changes."""

    def test_cache_invalid_on_version_change(self, cache: PreflightCache) -> None:
        with patch(
            "grimperium.cli.preflight._compute_settings_hash",
            return_value="",
        ):
            cache.save(
                checks_passed=["pandas"],
                grimperium_version="0.1.0",
            )
            assert cache.is_valid("0.2.0") is False


# ── Required Test 3 ──────────────────────────────────────────────────────


class TestAllOkSavesCache:
    """When all checks pass, cache is written."""

    def test_all_ok_saves_cache(
        self, runner: PreflightRunner, cache: PreflightCache
    ) -> None:
        with (
            patch.object(runner, "_run_all_checks", return_value=_ALL_OK),
            patch("grimperium.cli.preflight._get_version", return_value=_VERSION),
            patch(
                "grimperium.cli.preflight._compute_settings_hash",
                return_value="abc123",
            ),
        ):
            result = runner.run(force=True)

        assert result is True
        assert cache.cache_path.exists()
        data = json.loads(cache.cache_path.read_text())
        assert "pandas" in data["checks_passed"]
        assert data["grimperium_version"] == _VERSION
        assert data["settings_hash"] == "abc123"


# ── Required Test 4 ──────────────────────────────────────────────────────


class TestMissingPythonLibOffersFix:
    """Missing Python lib returns fix_available=True with pip command."""

    def test_missing_python_lib_offers_fix(self) -> None:
        with patch("grimperium.cli.preflight.importlib") as mock_importlib:
            mock_importlib.import_module.side_effect = ImportError("not found")
            result = check_python_lib("pandas", "pandas")

        assert result.status == "missing"
        assert result.fix_available is True
        assert "pip" in result.fix_command
        assert "pandas" in result.fix_command

    def test_sklearn_maps_to_scikit_learn(self) -> None:
        with patch("grimperium.cli.preflight.importlib") as mock_importlib:
            mock_importlib.import_module.side_effect = ImportError("not found")
            result = check_python_lib("sklearn", "scikit-learn")

        assert result.name == "scikit-learn"
        assert "scikit-learn" in result.fix_command


# ── Required Test 5 ──────────────────────────────────────────────────────


class TestMissingBinaryShowsHint:
    """Missing binary returns install_hint, fix_available=False."""

    def test_missing_binary_shows_hint(self) -> None:
        with patch(
            "grimperium.crest_pm7.validation.check_executable",
            return_value=(False, None),
        ):
            result = check_binary(
                "crest", "CREST", "--version", "Install CREST: ...", critical=True
            )

        assert result.status == "missing"
        assert result.fix_available is False
        assert result.install_hint != ""
        assert result.category == "binary"

    def test_optional_binary_not_critical(self) -> None:
        with patch(
            "grimperium.crest_pm7.validation.check_executable",
            return_value=(False, None),
        ):
            result = check_binary(
                "xtb", "xTB", "--version", "Install xTB: ...", critical=False
            )

        assert result.category == "binary_optional"
        assert result.is_critical is False


# ── Required Test 6 ──────────────────────────────────────────────────────


class TestCriticalFailureUserExits:
    """User chooses not to continue after critical failure -> returns False."""

    def test_critical_failure_user_exits(self, runner: PreflightRunner) -> None:
        failures = [
            CheckResult("CREST", "binary", "missing", "not found"),
            CheckResult("MOPAC", "binary", "missing", "not found"),
        ]
        with (
            patch.object(runner, "_run_all_checks", return_value=failures),
            patch.object(runner, "_rerun_failed", return_value=failures),
            patch("grimperium.cli.preflight._safe_confirm", return_value=False),
            patch("grimperium.cli.preflight._get_version", return_value=_VERSION),
        ):
            result = runner.run(force=True)

        assert result is False


# ── Required Test 7 ──────────────────────────────────────────────────────


class TestCriticalFailureUserContinues:
    """User chooses to continue despite critical failure -> returns True."""

    def test_critical_failure_user_continues(self, runner: PreflightRunner) -> None:
        failures = [
            CheckResult("CREST", "binary", "missing", "not found"),
        ]
        with (
            patch.object(runner, "_run_all_checks", return_value=failures),
            patch.object(runner, "_rerun_failed", return_value=failures),
            patch("grimperium.cli.preflight._safe_confirm", return_value=True),
            patch("grimperium.cli.preflight._get_version", return_value=_VERSION),
        ):
            result = runner.run(force=True)

        assert result is True


# ── Required Test 8 ──────────────────────────────────────────────────────


class TestSkipPreflightInAppRun:
    """When skip_preflight=True, preflight module is never called."""

    def test_skip_preflight_in_app_run(self) -> None:
        with patch("grimperium.cli.preflight.run_preflight") as mock_pf:
            from grimperium.cli.app import GrimperiumCLI

            app = GrimperiumCLI()
            # Make the loop exit immediately
            app.controller._running = False
            app.run(skip_preflight=True)

        mock_pf.assert_not_called()


# ── Edge-case tests ───────────────────────────────────────────────────────


class TestCorruptedCacheJson:
    """Corrupted cache JSON is treated as invalid."""

    def test_corrupted_json_is_invalid(self, cache: PreflightCache) -> None:
        cache.cache_dir.mkdir(parents=True, exist_ok=True)
        cache.cache_path.write_text("{corrupted json!@#$")
        assert cache.is_valid(_VERSION) is False

    def test_non_dict_json_is_invalid(self, cache: PreflightCache) -> None:
        cache.cache_dir.mkdir(parents=True, exist_ok=True)
        cache.cache_path.write_text('"just a string"')
        assert cache.is_valid(_VERSION) is False


class TestPipInstallFailure:
    """pip install failure is handled gracefully."""

    def test_pip_failure_returns_false(self, runner: PreflightRunner) -> None:
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "ERROR: No matching distribution"

        with patch(
            "grimperium.cli.preflight.subprocess.run",
            return_value=mock_result,
        ):
            success = runner._attempt_pip_install(["nonexistent-pkg"])

        assert success is False

    def test_pip_timeout_returns_false(self, runner: PreflightRunner) -> None:
        import subprocess

        with patch(
            "grimperium.cli.preflight.subprocess.run",
            side_effect=subprocess.TimeoutExpired(cmd="pip", timeout=120),
        ):
            success = runner._attempt_pip_install(["some-pkg"])

        assert success is False


class TestConfigWarningNotCritical:
    """Config issues are not critical."""

    def test_config_not_critical(self) -> None:
        result = CheckResult(
            "grimperium_settings.json",
            "config",
            "not_configured",
            "not found",
        )
        assert result.is_critical is False


class TestBinaryOptionalNotCritical:
    """binary_optional issues are not critical."""

    def test_binary_optional_not_critical(self) -> None:
        result = CheckResult("xTB", "binary_optional", "missing", "not found")
        assert result.is_critical is False


class TestForceIgnoresValidCache:
    """force=True runs checks even when cache is valid."""

    def test_force_ignores_cache(
        self, runner: PreflightRunner, cache: PreflightCache
    ) -> None:
        with patch(
            "grimperium.cli.preflight._compute_settings_hash",
            return_value="",
        ):
            cache.save(checks_passed=["pandas"], grimperium_version=_VERSION)

        with (
            patch.object(runner, "_run_all_checks", return_value=_ALL_OK),
            patch("grimperium.cli.preflight._get_version", return_value=_VERSION),
            patch(
                "grimperium.cli.preflight._compute_settings_hash",
                return_value="",
            ),
        ):
            result = runner.run(force=True)

        assert result is True

    def test_force_causes_check_execution(
        self, runner: PreflightRunner, cache: PreflightCache
    ) -> None:
        with patch(
            "grimperium.cli.preflight._compute_settings_hash",
            return_value="",
        ):
            cache.save(checks_passed=["pandas"], grimperium_version=_VERSION)

        with (
            patch.object(
                runner, "_run_all_checks", return_value=_ALL_OK
            ) as mock_checks,
            patch("grimperium.cli.preflight._get_version", return_value=_VERSION),
            patch(
                "grimperium.cli.preflight._compute_settings_hash",
                return_value="",
            ),
        ):
            runner.run(force=True)

        mock_checks.assert_called_once()


class TestNonTTYSafeConfirm:
    """Non-TTY stdin returns default without hanging."""

    def test_non_tty_returns_default(self) -> None:
        from grimperium.cli.preflight import _safe_confirm

        with patch("grimperium.cli.preflight.sys") as mock_sys:
            mock_sys.stdin.isatty.return_value = False
            assert _safe_confirm("Install?", default=False) is False
            assert _safe_confirm("Continue?", default=True) is True


class TestCacheInvalidOnSettingsHashChange:
    """Cache invalidates when settings file changes."""

    def test_settings_hash_change_invalidates(self, cache: PreflightCache) -> None:
        with patch(
            "grimperium.cli.preflight._compute_settings_hash",
            return_value="hash_v1",
        ):
            cache.save(checks_passed=["pandas"], grimperium_version=_VERSION)

        with patch(
            "grimperium.cli.preflight._compute_settings_hash",
            return_value="hash_v2",
        ):
            assert cache.is_valid(_VERSION) is False


class TestCacheInvalidOnPythonVersionChange:
    """Cache invalidates when Python version changes."""

    def test_python_version_change_invalidates(self, cache: PreflightCache) -> None:
        with patch(
            "grimperium.cli.preflight._compute_settings_hash",
            return_value="",
        ):
            cache.save(checks_passed=["pandas"], grimperium_version=_VERSION)

        with (
            patch(
                "grimperium.cli.preflight.platform.python_version",
                return_value="3.99.0",
            ),
            patch(
                "grimperium.cli.preflight._compute_settings_hash",
                return_value="",
            ),
        ):
            assert cache.is_valid(_VERSION) is False


class TestCheckSettingsFile:
    """Settings file check returns correct results."""

    def test_settings_file_exists(self, tmp_path: Path) -> None:
        settings_path = tmp_path / "grimperium_settings.json"
        settings_path.write_text("{}")

        with patch(
            "grimperium.cli.settings_manager.SettingsManager.default_config_path",
            return_value=settings_path,
        ):
            result = check_settings_file()

        assert result.is_ok
        assert result.category == "config"

    def test_settings_file_missing(self, tmp_path: Path) -> None:
        settings_path = tmp_path / "grimperium_settings.json"

        with patch(
            "grimperium.cli.settings_manager.SettingsManager.default_config_path",
            return_value=settings_path,
        ):
            result = check_settings_file()

        assert result.status == "not_configured"
        assert result.is_critical is False


class TestCheckPythonLibSuccess:
    """Successful Python lib check returns ok."""

    def test_lib_found(self) -> None:
        mock_mod = MagicMock()
        mock_mod.__version__ = "1.2.3"

        with patch("grimperium.cli.preflight.importlib") as mock_importlib:
            mock_importlib.import_module.return_value = mock_mod
            result = check_python_lib("pandas", "pandas")

        assert result.is_ok
        assert "1.2.3" in result.detail


class TestCheckBinarySuccess:
    """Successful binary check returns ok."""

    def test_binary_found(self) -> None:
        with patch(
            "grimperium.crest_pm7.validation.check_executable",
            return_value=(True, "CREST version 3.0"),
        ):
            result = check_binary("crest", "CREST", "--version", "hint", critical=True)

        assert result.is_ok
        assert "3.0" in result.detail
