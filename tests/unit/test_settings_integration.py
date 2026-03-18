"""Integration tests for settings persistence in controller and view.

Tests that:
- Controller auto-loads settings from file on init
- View reset_all deletes the persisted config file
- View exposes "Save as Default" action
"""

import json
from pathlib import Path
from unittest.mock import patch

import pytest
from rich.console import Console

from grimperium.cli.controller import CliController
from grimperium.cli.settings_manager import SettingsManager
from grimperium.cli.views.settings_view import SettingsView


@pytest.fixture
def config_path(tmp_path: Path) -> Path:
    """Provide a temporary config file path."""
    return tmp_path / "grimperium_settings.json"


class TestControllerAutoLoad:
    """Controller must auto-load settings from config file on init."""

    def test_controller_loads_existing_config(
        self, config_path: Path
    ) -> None:
        """If config file exists, controller must load it at init."""
        # Prepare a config file with non-default values
        data = SettingsManager(console=Console(force_terminal=False)).to_dict()
        data["crest_ewin"] = 42.0
        data["mopac_precise"] = True
        config_path.write_text(json.dumps(data))

        with patch.object(
            SettingsManager, "default_config_path", return_value=config_path
        ):
            ctrl = CliController()

        assert ctrl.settings_manager.crest.ewin == 42.0
        assert ctrl.settings_manager.mopac.precise is True

    def test_controller_works_without_config_file(
        self, config_path: Path
    ) -> None:
        """Without config file, controller uses defaults."""
        with patch.object(
            SettingsManager, "default_config_path", return_value=config_path
        ):
            ctrl = CliController()

        assert ctrl.settings_manager.crest.ewin == 5.0
        assert ctrl.settings_manager.mopac.precise is False


class TestViewSaveAsDefault:
    """SettingsView must support 'save_as_default' action."""

    def test_save_as_default_persists_to_file(
        self, config_path: Path
    ) -> None:
        """'save_as_default' action must write settings to disk."""
        with patch.object(
            SettingsManager, "default_config_path", return_value=config_path
        ):
            ctrl = CliController()
            view = SettingsView(ctrl)
            ctrl.settings_manager.crest.ewin = 77.0

            with patch.object(view, "wait_for_enter"):
                result = view.handle_action("save_as_default")

        assert config_path.exists()
        data = json.loads(config_path.read_text())
        assert data["crest_ewin"] == 77.0
        assert result is None  # stays in settings view

    def test_menu_options_include_save_as_default(self) -> None:
        """Menu options must include 'Save as Default'."""
        ctrl = CliController()
        view = SettingsView(ctrl)
        options = view.get_menu_options()
        values = [opt.value for opt in options]
        assert "save_as_default" in values


class TestViewResetDeletesFile:
    """View reset_all must also delete persisted config file."""

    def test_reset_all_deletes_config_file(
        self, config_path: Path
    ) -> None:
        """'reset_all' action must remove the persisted config file."""
        # Create a persisted config
        config_path.write_text(json.dumps({"crest_ewin": 99.0}))

        with patch.object(
            SettingsManager, "default_config_path", return_value=config_path
        ):
            ctrl = CliController()
            view = SettingsView(ctrl)

            # Simulate reset with auto-confirm and skip wait_for_enter
            with (
                patch("grimperium.cli.menu.confirm", return_value=True),
                patch.object(view, "wait_for_enter"),
            ):
                view.handle_action("reset_all")

        assert not config_path.exists()
        assert ctrl.settings_manager.crest.ewin == 5.0
