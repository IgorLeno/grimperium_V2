"""Unit tests for SettingsManager persistence (save/load to JSON).

Tests the file-based persistence layer for settings configuration.
TDD RED phase — these tests should all FAIL before implementation.
"""

import json
from pathlib import Path

import pytest
from rich.console import Console

from grimperium.cli.settings_manager import SettingsManager


@pytest.fixture
def manager() -> SettingsManager:
    """Create a SettingsManager instance for testing."""
    return SettingsManager(console=Console(force_terminal=False))


@pytest.fixture
def config_path(tmp_path: Path) -> Path:
    """Provide a temporary config file path."""
    return tmp_path / "grimperium_settings.json"


class TestSaveToFile:
    """Tests for SettingsManager.save_to_file()."""

    def test_save_creates_json_file(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """save_to_file() must create a JSON file at the given path."""
        manager.save_to_file(config_path)
        assert config_path.exists()

    def test_save_writes_valid_json(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """Saved file must contain valid, parseable JSON."""
        manager.save_to_file(config_path)
        data = json.loads(config_path.read_text())
        assert isinstance(data, dict)

    def test_save_contains_all_keys(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """Saved JSON must contain all keys from to_dict()."""
        manager.save_to_file(config_path)
        data = json.loads(config_path.read_text())
        expected_keys = set(manager.to_dict().keys())
        assert set(data.keys()) == expected_keys

    def test_save_preserves_modified_values(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """Modified settings must be reflected in the saved file."""
        manager.crest.ewin = 12.5
        manager.mopac.precise = True
        manager.xtb.preopt = False

        manager.save_to_file(config_path)
        data = json.loads(config_path.read_text())

        assert data["crest_ewin"] == 12.5
        assert data["mopac_precise"] is True
        assert data["crest_xtb_preopt"] is False

    def test_save_overwrites_existing_file(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """Saving again must overwrite the previous file."""
        manager.crest.ewin = 5.0
        manager.save_to_file(config_path)

        manager.crest.ewin = 99.0
        manager.save_to_file(config_path)

        data = json.loads(config_path.read_text())
        assert data["crest_ewin"] == 99.0

    def test_save_returns_true_on_success(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """save_to_file() must return True on success."""
        result = manager.save_to_file(config_path)
        assert result is True

    def test_save_creates_parent_directories(
        self, manager: SettingsManager, tmp_path: Path
    ) -> None:
        """save_to_file() must create parent dirs if they don't exist."""
        nested_path = tmp_path / "subdir" / "deep" / "grimperium_settings.json"
        result = manager.save_to_file(nested_path)
        assert result is True
        assert nested_path.exists()


class TestLoadFromFile:
    """Tests for SettingsManager.load_from_file()."""

    def test_load_restores_saved_settings(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """Roundtrip: save → fresh manager → load must restore values."""
        manager.crest.ewin = 7.5
        manager.crest.optlev = "vtight"
        manager.mopac.scfcrt = 1.0e-6
        manager.mopac.precise = True
        manager.xtb.preopt = False
        manager.save_to_file(config_path)

        fresh = SettingsManager(console=Console(force_terminal=False))
        result = fresh.load_from_file(config_path)

        assert result is True
        assert fresh.crest.ewin == 7.5
        assert fresh.crest.optlev == "vtight"
        assert fresh.mopac.scfcrt == 1.0e-6
        assert fresh.mopac.precise is True
        assert fresh.xtb.preopt is False

    def test_load_nonexistent_file_returns_false(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """load_from_file() with missing file must return False."""
        result = manager.load_from_file(config_path)
        assert result is False

    def test_load_nonexistent_preserves_defaults(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """Loading a missing file must not alter default settings."""
        defaults = manager.to_dict()
        manager.load_from_file(config_path)
        assert manager.to_dict() == defaults

    def test_load_corrupted_json_returns_false(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """load_from_file() with invalid JSON must return False."""
        config_path.write_text("{invalid json content!!!")
        result = manager.load_from_file(config_path)
        assert result is False

    def test_load_corrupted_json_preserves_defaults(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """Loading corrupted JSON must not alter default settings."""
        defaults = manager.to_dict()
        config_path.write_text("not json")
        manager.load_from_file(config_path)
        assert manager.to_dict() == defaults

    def test_load_partial_json_applies_available_keys(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """Loading JSON with only some keys must apply those keys only."""
        config_path.write_text(json.dumps({"crest_ewin": 42.0}))
        result = manager.load_from_file(config_path)
        assert result is True
        assert manager.crest.ewin == 42.0
        # Other settings remain default
        assert manager.crest.v3 is True
        assert manager.mopac.precise is False


class TestDeleteConfigFile:
    """Tests for SettingsManager.delete_config_file()."""

    def test_delete_removes_existing_file(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """delete_config_file() must remove the file."""
        manager.save_to_file(config_path)
        assert config_path.exists()

        manager.delete_config_file(config_path)
        assert not config_path.exists()

    def test_delete_nonexistent_no_error(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """delete_config_file() with missing file must not raise."""
        assert not config_path.exists()
        manager.delete_config_file(config_path)  # Should not raise

    def test_delete_returns_true_when_file_existed(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """delete_config_file() must return True if file was removed."""
        manager.save_to_file(config_path)
        result = manager.delete_config_file(config_path)
        assert result is True

    def test_delete_returns_false_when_no_file(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """delete_config_file() must return False if file didn't exist."""
        result = manager.delete_config_file(config_path)
        assert result is False


class TestDefaultConfigPath:
    """Tests for SettingsManager.default_config_path()."""

    def test_returns_path_object(self, manager: SettingsManager) -> None:
        """default_config_path() must return a Path, not a string."""
        path = manager.default_config_path()
        assert isinstance(path, Path)

    def test_filename_is_grimperium_settings_json(
        self, manager: SettingsManager
    ) -> None:
        """Config file must be named grimperium_settings.json."""
        path = manager.default_config_path()
        assert path.name == "grimperium_settings.json"


class TestResetAllWithPersistence:
    """Integration: reset_all + delete_config_file."""

    def test_reset_and_delete_restores_defaults(
        self, manager: SettingsManager, config_path: Path
    ) -> None:
        """After reset_all + delete, settings are default and file is gone."""
        manager.crest.ewin = 99.0
        manager.mopac.precise = True
        manager.save_to_file(config_path)

        manager.reset_all()
        manager.delete_config_file(config_path)

        assert manager.crest.ewin == 5.0
        assert manager.mopac.precise is False
        assert not config_path.exists()
