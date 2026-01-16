"""Unit tests for SettingsManager.

Tests the CREST, MOPAC, and xTB settings management functionality.
"""

import pytest
from rich.console import Console

from grimperium.cli.settings_manager import (
    CRESTSettings,
    MOPACSettings,
    SettingsManager,
    xTBSettings,
)


class TestCRESTSettings:
    """Tests for CRESTSettings dataclass."""

    def test_defaults(self) -> None:
        """Verify CREST default values."""
        settings = CRESTSettings()
        assert settings.v3 is True
        assert settings.quick is False
        assert settings.nci is False
        assert settings.gfnff is False
        assert settings.ewin == 5.0
        assert settings.rthr == 0.125
        assert settings.optlev == "normal"
        assert settings.threads == 4

    def test_optlev_choices(self) -> None:
        """Verify OPTLEV_CHOICES contains expected values."""
        expected = ["loose", "normal", "tight", "vtight", "extreme"]
        assert expected == CRESTSettings.OPTLEV_CHOICES

    def test_custom_values(self) -> None:
        """Test setting custom values."""
        settings = CRESTSettings(
            v3=False,
            quick=True,
            ewin=10.0,
            threads=8,
        )
        assert settings.v3 is False
        assert settings.quick is True
        assert settings.ewin == 10.0
        assert settings.threads == 8


class TestMOPACSettings:
    """Tests for MOPACSettings dataclass."""

    def test_defaults(self) -> None:
        """Verify MOPAC default values."""
        settings = MOPACSettings()
        assert settings.precise is False
        assert settings.scfcrt == 1.0e-4
        assert settings.itry == 1000
        assert settings.pulay is False
        assert settings.prtall is False
        assert settings.archive is False

    def test_custom_values(self) -> None:
        """Test setting custom values."""
        settings = MOPACSettings(
            precise=True,
            scfcrt=1.0e-6,
            itry=2000,
        )
        assert settings.precise is True
        assert settings.scfcrt == 1.0e-6
        assert settings.itry == 2000


class TestxTBSettings:
    """Tests for xTBSettings dataclass."""

    def test_defaults(self) -> None:
        """Verify xTB default values (disabled by default)."""
        settings = xTBSettings()
        assert settings.preopt is False
        assert settings.timeout_seconds == 300

    def test_enable_preopt(self) -> None:
        """Test enabling pre-optimization."""
        settings = xTBSettings(preopt=True)
        assert settings.preopt is True


class TestSettingsManager:
    """Tests for SettingsManager class."""

    @pytest.fixture
    def manager(self) -> SettingsManager:
        """Create a SettingsManager instance for testing."""
        return SettingsManager(console=Console(force_terminal=False))

    def test_init_defaults(self, manager: SettingsManager) -> None:
        """Verify manager initializes with default settings."""
        assert manager.crest.v3 is True
        assert manager.mopac.precise is False
        assert manager.xtb.preopt is False

    def test_to_dict_contains_all_keys(self, manager: SettingsManager) -> None:
        """Verify to_dict() includes all expected settings keys."""
        settings_dict = manager.to_dict()

        expected_keys = [
            "crest_v3",
            "crest_quick",
            "crest_nci",
            "crest_gfnff",
            "crest_ewin",
            "crest_rthr",
            "crest_optlev",
            "crest_threads",
            "mopac_precise",
            "mopac_scfcrt",
            "mopac_itry",
            "mopac_pulay",
            "mopac_prtall",
            "mopac_archive",
            "crest_xtb_preopt",
        ]

        for key in expected_keys:
            assert key in settings_dict, f"Missing key: {key}"

    def test_to_dict_default_values(self, manager: SettingsManager) -> None:
        """Verify to_dict() returns correct default values."""
        settings_dict = manager.to_dict()

        assert settings_dict["crest_v3"] is True
        assert settings_dict["crest_quick"] is False
        assert settings_dict["crest_ewin"] == 5.0
        assert settings_dict["crest_threads"] == 4
        assert settings_dict["mopac_precise"] is False
        assert settings_dict["mopac_itry"] == 1000
        assert settings_dict["crest_xtb_preopt"] is False

    def test_to_dict_after_modification(self, manager: SettingsManager) -> None:
        """Verify to_dict() reflects modified values."""
        manager.crest.ewin = 10.0
        manager.mopac.precise = True
        manager.xtb.preopt = True

        settings_dict = manager.to_dict()

        assert settings_dict["crest_ewin"] == 10.0
        assert settings_dict["mopac_precise"] is True
        assert settings_dict["crest_xtb_preopt"] is True

    def test_from_dict_loads_settings(self, manager: SettingsManager) -> None:
        """Verify from_dict() correctly loads settings."""
        input_dict = {
            "crest_v3": False,
            "crest_ewin": 8.5,
            "crest_optlev": "tight",
            "mopac_precise": True,
            "mopac_itry": 1500,
            "crest_xtb_preopt": True,
        }

        manager.from_dict(input_dict)

        assert manager.crest.v3 is False
        assert manager.crest.ewin == 8.5
        assert manager.crest.optlev == "tight"
        assert manager.mopac.precise is True
        assert manager.mopac.itry == 1500
        assert manager.xtb.preopt is True

    def test_from_dict_ignores_invalid_optlev(self, manager: SettingsManager) -> None:
        """Verify from_dict() ignores invalid optlev values."""
        original_optlev = manager.crest.optlev
        manager.from_dict({"crest_optlev": "invalid_value"})
        assert manager.crest.optlev == original_optlev

    def test_from_dict_partial_update(self, manager: SettingsManager) -> None:
        """Verify from_dict() only updates provided keys."""
        original_v3 = manager.crest.v3
        manager.from_dict({"mopac_precise": True})

        assert manager.crest.v3 == original_v3
        assert manager.mopac.precise is True

    def test_reset_crest(self, manager: SettingsManager) -> None:
        """Verify reset_crest() restores defaults."""
        manager.crest.ewin = 999.0
        manager.crest.v3 = False
        manager.reset_crest()

        assert manager.crest.ewin == 5.0
        assert manager.crest.v3 is True

    def test_reset_mopac(self, manager: SettingsManager) -> None:
        """Verify reset_mopac() restores defaults."""
        manager.mopac.precise = True
        manager.mopac.itry = 9999
        manager.reset_mopac()

        assert manager.mopac.precise is False
        assert manager.mopac.itry == 1000

    def test_reset_xtb(self, manager: SettingsManager) -> None:
        """Verify reset_xtb() restores defaults."""
        manager.xtb.preopt = True
        manager.reset_xtb()

        assert manager.xtb.preopt is False

    def test_reset_all(self, manager: SettingsManager) -> None:
        """Verify reset_all() restores all defaults."""
        manager.crest.ewin = 999.0
        manager.mopac.precise = True
        manager.xtb.preopt = True

        manager.reset_all()

        assert manager.crest.ewin == 5.0
        assert manager.mopac.precise is False
        assert manager.xtb.preopt is False

    def test_help_text_coverage(self, manager: SettingsManager) -> None:
        """Verify help text exists for all expected settings."""
        expected_help_keys = [
            "crest_v3",
            "crest_quick",
            "crest_nci",
            "crest_gfnff",
            "crest_ewin",
            "crest_rthr",
            "crest_optlev",
            "crest_threads",
            "mopac_precise",
            "mopac_scfcrt",
            "mopac_itry",
            "mopac_pulay",
            "mopac_prtall",
            "mopac_archive",
            "xtb_preopt",
        ]

        for key in expected_help_keys:
            assert key in manager.HELP_TEXT, f"Missing help text for: {key}"
            assert len(manager.HELP_TEXT[key]) > 0, f"Empty help text for: {key}"

    def test_show_crest_summary_returns_table(self, manager: SettingsManager) -> None:
        """Verify show_crest_summary() returns a Rich Table."""
        from rich.table import Table

        table = manager.show_crest_summary()
        assert isinstance(table, Table)

    def test_show_mopac_summary_returns_table(self, manager: SettingsManager) -> None:
        """Verify show_mopac_summary() returns a Rich Table."""
        from rich.table import Table

        table = manager.show_mopac_summary()
        assert isinstance(table, Table)

    def test_show_xtb_summary_returns_table(self, manager: SettingsManager) -> None:
        """Verify show_xtb_summary() returns a Rich Table."""
        from rich.table import Table

        table = manager.show_xtb_summary()
        assert isinstance(table, Table)

    def test_roundtrip_to_from_dict(self, manager: SettingsManager) -> None:
        """Verify settings survive a to_dict/from_dict roundtrip."""
        manager.crest.ewin = 7.5
        manager.crest.optlev = "vtight"
        manager.mopac.scfcrt = 1.0e-6
        manager.xtb.preopt = True

        saved = manager.to_dict()

        new_manager = SettingsManager(console=Console(force_terminal=False))
        new_manager.from_dict(saved)

        assert new_manager.crest.ewin == 7.5
        assert new_manager.crest.optlev == "vtight"
        assert new_manager.mopac.scfcrt == 1.0e-6
        assert new_manager.xtb.preopt is True
