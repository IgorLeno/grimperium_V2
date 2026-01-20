"""Tests for Phase A settings system changes."""

import pytest

from src.grimperium.cli.settings_manager import CRESTSettings, SettingsManager


def test_crest_settings_has_new_fields() -> None:
    """Test that CRESTSettings has crest_method and quick_mode fields."""
    settings = CRESTSettings()

    # New fields should exist with defaults
    assert hasattr(settings, "crest_method")
    assert hasattr(settings, "quick_mode")
    assert settings.crest_method == "gfn2"
    assert settings.quick_mode == "off"

    # Old fields should NOT exist
    assert not hasattr(settings, "gfnff")
    assert not hasattr(settings, "quick")


def test_crest_method_options() -> None:
    """Test that CREST_METHOD_OPTIONS are defined."""
    assert hasattr(CRESTSettings, "CREST_METHOD_OPTIONS")
    assert CRESTSettings.CREST_METHOD_OPTIONS == ["gfn2", "gfnff", "gfn2//gfnff"]


def test_quick_mode_options() -> None:
    """Test that QUICK_MODE_OPTIONS are defined."""
    assert hasattr(CRESTSettings, "QUICK_MODE_OPTIONS")
    assert CRESTSettings.QUICK_MODE_OPTIONS == ["off", "quick", "squick", "mquick"]


def test_settings_to_dict_new_format() -> None:
    """Test that to_dict exports new fields correctly."""
    sm = SettingsManager()
    sm.crest.crest_method = "gfnff"
    sm.crest.quick_mode = "quick"

    result = sm.to_dict()

    # New fields should be present
    assert "crest_method" in result
    assert "crest_quick_mode" in result
    assert result["crest_method"] == "gfnff"
    assert result["crest_quick_mode"] == "quick"

    # Old fields should NOT be present
    assert "crest_gfnff" not in result
    assert "crest_quick" not in result
