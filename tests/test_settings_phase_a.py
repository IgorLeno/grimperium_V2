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


def test_from_dict_backward_compat_gfnff_true() -> None:
    """Test loading old format with gfnff=true converts to crest_method=gfnff."""
    sm = SettingsManager()

    old_data = {"crest_gfnff": "true"}
    sm.from_dict(old_data)

    assert sm.crest.crest_method == "gfnff"


def test_from_dict_backward_compat_gfnff_false() -> None:
    """Test loading old format with gfnff=false converts to crest_method=gfn2."""
    sm = SettingsManager()

    old_data = {"crest_gfnff": "false"}
    sm.from_dict(old_data)

    assert sm.crest.crest_method == "gfn2"


def test_from_dict_backward_compat_quick_true() -> None:
    """Test loading old format with quick=true converts to quick_mode=quick."""
    sm = SettingsManager()

    old_data = {"crest_quick": "true"}
    sm.from_dict(old_data)

    assert sm.crest.quick_mode == "quick"


def test_from_dict_backward_compat_quick_false() -> None:
    """Test loading old format with quick=false converts to quick_mode=off."""
    sm = SettingsManager()

    old_data = {"crest_quick": "false"}
    sm.from_dict(old_data)

    assert sm.crest.quick_mode == "off"


def test_from_dict_new_format() -> None:
    """Test loading new format works directly."""
    sm = SettingsManager()

    new_data = {"crest_method": "gfn2//gfnff", "crest_quick_mode": "squick"}
    sm.from_dict(new_data)

    assert sm.crest.crest_method == "gfn2//gfnff"
    assert sm.crest.quick_mode == "squick"


def test_settings_roundtrip() -> None:
    """Test that to_dict -> from_dict preserves settings."""
    sm1 = SettingsManager()
    sm1.crest.crest_method = "gfnff"
    sm1.crest.quick_mode = "mquick"

    data = sm1.to_dict()

    sm2 = SettingsManager()
    sm2.from_dict(data)

    assert sm2.crest.crest_method == "gfnff"
    assert sm2.crest.quick_mode == "mquick"
