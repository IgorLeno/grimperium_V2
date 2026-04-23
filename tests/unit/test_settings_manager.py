"""Unit tests for SettingsManager.

Tests the CREST, MOPAC, and xTB settings management functionality.
"""

import pytest
from rich.console import Console

from grimperium.cli.settings_manager import (
    CalculationProfile,
    CRESTSettings,
    DistributedDefaults,
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
        assert settings.nci is False
        assert settings.crest_method == "gfn2"
        assert settings.quick_mode == "off"
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
            quick_mode="quick",
            crest_method="gfnff",
            ewin=10.0,
            threads=8,
        )
        assert settings.v3 is False
        assert settings.quick_mode == "quick"
        assert settings.crest_method == "gfnff"
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
        assert settings.preopt is True
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
        assert manager.xtb.preopt is True

    def test_to_dict_contains_all_keys(self, manager: SettingsManager) -> None:
        """Verify to_dict() includes all expected settings keys."""
        settings_dict = manager.to_dict()

        expected_keys = [
            "crest_v3",
            "crest_nci",
            "crest_method",
            "crest_quick_mode",
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
        assert settings_dict["crest_method"] == "gfn2"
        assert settings_dict["crest_quick_mode"] == "off"
        assert settings_dict["crest_ewin"] == 5.0
        assert settings_dict["crest_threads"] == 4
        assert settings_dict["mopac_precise"] is False
        assert settings_dict["mopac_itry"] == 1000
        assert settings_dict["crest_xtb_preopt"] is True

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

        assert manager.xtb.preopt is True

    def test_reset_all(self, manager: SettingsManager) -> None:
        """Verify reset_all() restores all defaults."""
        manager.crest.ewin = 999.0
        manager.mopac.precise = True
        manager.xtb.preopt = True

        manager.reset_all()

        assert manager.crest.ewin == 5.0
        assert manager.mopac.precise is False
        assert manager.xtb.preopt is True

    def test_help_text_coverage(self, manager: SettingsManager) -> None:
        """Verify help text exists for all expected settings."""
        expected_help_keys = [
            "crest_v3",
            "crest_nci",
            "crest_method",
            "crest_quick_mode",
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


class TestCalculationProfile:
    """Tests for CalculationProfile dataclass."""

    def test_defaults(self) -> None:
        p = CalculationProfile()
        assert p.name == "default"
        assert p.crest_ewin == 6.0
        assert p.crest_rthr == 0.125
        assert p.crest_opt == 2
        assert p.crest_threads == 4
        assert p.crest_timeout_minutes == 60
        assert p.mopac_keywords == "PM7 EF GNORM=0.01"
        assert p.mopac_timeout_minutes == 30
        assert p.is_standard is False
        assert p.created_at != ""

    def test_created_at_auto_set(self) -> None:
        p = CalculationProfile()
        assert "T" in p.created_at  # ISO datetime format

    def test_to_dict_round_trip(self) -> None:
        p = CalculationProfile(name="fast", crest_timeout_minutes=30, is_standard=True)
        d = p.to_dict()
        restored = CalculationProfile.from_dict(d)
        assert restored.name == "fast"
        assert restored.crest_timeout_minutes == 30
        assert restored.is_standard is True

    def test_from_dict_ignores_unknown_keys(self) -> None:
        p = CalculationProfile.from_dict({"name": "x", "unknown_field": 99})
        assert p.name == "x"

    def test_custom_values_preserved(self) -> None:
        p = CalculationProfile(
            name="heavy",
            crest_ewin=10.0,
            crest_threads=8,
            mopac_timeout_minutes=120,
        )
        restored = CalculationProfile.from_dict(p.to_dict())
        assert restored.crest_ewin == 10.0
        assert restored.crest_threads == 8
        assert restored.mopac_timeout_minutes == 120


class TestDistributedDefaults:
    """Tests for DistributedDefaults dataclass."""

    def test_defaults(self) -> None:
        d = DistributedDefaults()
        assert d.profile_name == "default"
        assert d.batch_size == 10
        assert d.crest_timeout_minutes == 60
        assert d.mopac_timeout_minutes == 30

    def test_to_dict_round_trip(self) -> None:
        d = DistributedDefaults(profile_name="fast", batch_size=5)
        restored = DistributedDefaults.from_dict(d.to_dict())
        assert restored.profile_name == "fast"
        assert restored.batch_size == 5

    def test_from_dict_ignores_unknown_keys(self) -> None:
        d = DistributedDefaults.from_dict({"batch_size": 20, "extra": "ignored"})
        assert d.batch_size == 20


class TestDistributedPersistence:
    """Tests for profile and defaults persistence in ~/.grimperium/."""

    def test_save_and_load_profiles(self, tmp_path: pytest.TempPathFactory) -> None:
        profiles_file = tmp_path / "profiles.json"  # type: ignore[operator]
        profiles = [
            CalculationProfile(name="default", is_standard=True),
            CalculationProfile(name="fast", crest_timeout_minutes=30),
        ]
        profiles_file.write_text(  # type: ignore[union-attr]
            __import__("json").dumps({"profiles": [p.to_dict() for p in profiles]},
                                    indent=2),
            encoding="utf-8",
        )
        raw = __import__("json").loads(profiles_file.read_text(encoding="utf-8"))  # type: ignore[union-attr]
        loaded = [CalculationProfile.from_dict(p) for p in raw["profiles"]]
        assert len(loaded) == 2
        assert loaded[0].name == "default"
        assert loaded[1].crest_timeout_minutes == 30

    def test_save_and_load_distributed_defaults(
        self, tmp_path: pytest.TempPathFactory
    ) -> None:
        defaults_file = tmp_path / "distributed_defaults.json"  # type: ignore[operator]
        d = DistributedDefaults(profile_name="fast", batch_size=5)
        defaults_file.write_text(  # type: ignore[union-attr]
            __import__("json").dumps(d.to_dict(), indent=2), encoding="utf-8"
        )
        raw = __import__("json").loads(defaults_file.read_text(encoding="utf-8"))  # type: ignore[union-attr]
        restored = DistributedDefaults.from_dict(raw)
        assert restored.profile_name == "fast"
        assert restored.batch_size == 5

    def test_load_profiles_auto_creates_default(
        self, tmp_path: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            SettingsManager,
            "grimperium_home_dir",
            staticmethod(lambda: tmp_path),  # type: ignore[arg-type]
        )
        profiles = SettingsManager.load_profiles()
        assert any(p.name == "default" for p in profiles)

    def test_load_profiles_default_not_duplicated(
        self, tmp_path: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            SettingsManager,
            "grimperium_home_dir",
            staticmethod(lambda: tmp_path),  # type: ignore[arg-type]
        )
        SettingsManager.save_profiles(
            [CalculationProfile(name="default", is_standard=True)]
        )
        profiles = SettingsManager.load_profiles()
        assert sum(1 for p in profiles if p.name == "default") == 1

    def test_save_profiles_returns_true(
        self, tmp_path: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            SettingsManager,
            "grimperium_home_dir",
            staticmethod(lambda: tmp_path),  # type: ignore[arg-type]
        )
        result = SettingsManager.save_profiles(
            [CalculationProfile(name="default", is_standard=True)]
        )
        assert result is True

    def test_load_distributed_defaults_fallback(
        self, tmp_path: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            SettingsManager,
            "grimperium_home_dir",
            staticmethod(lambda: tmp_path),  # type: ignore[arg-type]
        )
        d = SettingsManager.load_distributed_defaults()
        assert d.profile_name == "default"
        assert d.batch_size == 10

    def test_save_and_load_distributed_defaults_roundtrip(
        self, tmp_path: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            SettingsManager,
            "grimperium_home_dir",
            staticmethod(lambda: tmp_path),  # type: ignore[arg-type]
        )
        d = DistributedDefaults(profile_name="heavy", batch_size=3)
        SettingsManager.save_distributed_defaults(d)
        loaded = SettingsManager.load_distributed_defaults()
        assert loaded.profile_name == "heavy"
        assert loaded.batch_size == 3
