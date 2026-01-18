"""Tests for settings view module."""

from unittest.mock import MagicMock, patch

import pytest

from grimperium.cli.views.settings_view import SettingsView


@pytest.fixture
def mock_controller() -> MagicMock:
    """Create mock controller for SettingsView."""
    controller = MagicMock()
    return controller


@pytest.fixture
def settings_view(mock_controller: MagicMock) -> SettingsView:
    """Create SettingsView instance."""
    return SettingsView(mock_controller)


def test_display_all_settings_no_duplicate_headers(settings_view: SettingsView) -> None:
    """
    Bug #5: Test that _display_all_settings doesn't show duplicate headers.

    The top-level "⚙️ Current Settings" header should be removed.
    Each section should appear with its Panel title only.
    """
    # Mock the console to capture output
    with patch.object(settings_view, "console") as mock_console:
        mock_console.print = MagicMock()

        # Mock settings manager methods
        settings_view._settings_manager = MagicMock()
        settings_view._settings_manager.show_crest_summary.return_value = (
            "CREST summary"
        )
        settings_view._settings_manager.show_mopac_summary.return_value = (
            "MOPAC summary"
        )
        settings_view._settings_manager.show_xtb_summary.return_value = "xTB summary"

        settings_view._display_all_settings()

        # Get all print calls
        print_calls = [str(call) for call in mock_console.print.call_args_list]
        print_output = " ".join(print_calls)

        # Should NOT have the top-level "Current Settings" header
        # Fixed: explicit check without `or` logic that could produce false positives
        assert (
            "⚙️  Current Settings" not in print_output
        ), "Top-level 'Current Settings' header should be removed"

        # Note: Panel titles ("CREST Settings", "MOPAC Settings", "xTB Settings")
        # are Rich objects and not easily checked in string representation.
        # The absence of "Current Settings" header is the primary bug fix validation.


def test_display_all_settings_shows_all_panels(settings_view: SettingsView) -> None:
    """
    Bug #5: Test that all settings panels are displayed exactly once.

    Should show CREST, MOPAC, and xTB panels with titles.
    """
    with patch.object(settings_view, "console") as mock_console:
        settings_view._settings_manager = MagicMock()
        settings_view._settings_manager.show_crest_summary.return_value = (
            "CREST summary"
        )
        settings_view._settings_manager.show_mopac_summary.return_value = (
            "MOPAC summary"
        )
        settings_view._settings_manager.show_xtb_summary.return_value = "xTB summary"

        settings_view._display_all_settings()

        # Verify all three summary methods were called exactly once
        settings_view._settings_manager.show_crest_summary.assert_called_once()
        settings_view._settings_manager.show_mopac_summary.assert_called_once()
        settings_view._settings_manager.show_xtb_summary.assert_called_once()

        # Verify console.print was called (panels were displayed)
        assert (
            mock_console.print.call_count > 0
        ), "Should display panels via console.print"


def test_get_menu_options_includes_back(settings_view: SettingsView) -> None:
    """
    Bug #6: Test that menu options include a Back button.

    The settings menu should have a "Back to Main Menu" option.
    """
    options = settings_view.get_menu_options()
    assert len(options) > 0, "Should have menu options"

    # Check if back option exists
    option_values = [opt.value for opt in options]
    assert "back" in option_values, "Menu should include 'back' option"

    # Find the back option
    back_option = next((opt for opt in options if opt.value == "back"), None)
    assert back_option is not None, "Back option should exist"
    assert (
        "back" in back_option.label.lower() or "main" in back_option.label.lower()
    ), "Back option label should indicate return to main menu"


def test_handle_action_back_returns_main(settings_view: SettingsView) -> None:
    """
    Bug #6: Test that selecting 'back' action returns to main menu.

    Should return "main" to navigate back to main menu.
    """
    result = settings_view.handle_action("back")
    assert result == "main", "Back action should return 'main'"


def test_handle_action_other_stays_in_settings(settings_view: SettingsView) -> None:
    """
    Bug #6: Test that other actions don't exit settings view.

    Actions like "view_all" should return None to stay in settings.
    """
    with patch.object(settings_view, "console"):
        settings_view._settings_manager = MagicMock()
        settings_view._settings_manager.show_crest_summary.return_value = "CREST"
        settings_view._settings_manager.show_mopac_summary.return_value = "MOPAC"
        settings_view._settings_manager.show_xtb_summary.return_value = "xTB"

        result = settings_view.handle_action("view_all")

    assert result is None, "Non-back actions should return None (stay in view)"
