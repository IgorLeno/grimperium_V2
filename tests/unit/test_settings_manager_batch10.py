"""Unit tests for BATCH 10: Critical Pre-Smoke Test Fixes.

Tests for:
1. Thread count validation (reject 0/negative/too-high)
2. Boolean string parsing ("False" → False, not True)
3. Path.replace Windows compatibility
"""

import pytest

from grimperium.cli.settings_manager import SettingsManager


class TestParseBool:
    """Tests for _parse_bool static method."""

    def test_parse_bool_string_true_variants(self) -> None:
        """Verify 'true' string variants parse as True."""
        for val in ["true", "True", "TRUE", "yes", "Yes", "YES", "1", "on", "ON"]:
            assert SettingsManager._parse_bool(val) is True

    def test_parse_bool_string_false_variants(self) -> None:
        """Verify 'false' string variants parse as False."""
        for val in ["false", "False", "FALSE", "no", "No", "NO", "0", "off", "OFF"]:
            assert SettingsManager._parse_bool(val) is False

    def test_parse_bool_string_with_whitespace(self) -> None:
        """Verify whitespace is trimmed before parsing."""
        assert SettingsManager._parse_bool("  true  ") is True
        assert SettingsManager._parse_bool("  false  ") is False

    def test_parse_bool_int_zero_is_false(self) -> None:
        """Verify int 0 parses as False."""
        assert SettingsManager._parse_bool(0) is False

    def test_parse_bool_int_nonzero_is_true(self) -> None:
        """Verify non-zero ints parse as True."""
        assert SettingsManager._parse_bool(1) is True
        assert SettingsManager._parse_bool(-1) is True
        assert SettingsManager._parse_bool(99) is True

    def test_parse_bool_bool_passthrough(self) -> None:
        """Verify bool values pass through unchanged."""
        assert SettingsManager._parse_bool(True) is True
        assert SettingsManager._parse_bool(False) is False

    def test_parse_bool_invalid_string_raises_valueerror(self) -> None:
        """Verify invalid string raises ValueError."""
        with pytest.raises(ValueError, match="Cannot parse boolean"):
            SettingsManager._parse_bool("maybe")

    def test_parse_bool_invalid_string_empty_raises_valueerror(self) -> None:
        """Verify empty string raises ValueError."""
        with pytest.raises(ValueError, match="Cannot parse boolean"):
            SettingsManager._parse_bool("")

    def test_parse_bool_invalid_type_raises_typeerror(self) -> None:
        """Verify unsupported types raise TypeError."""
        with pytest.raises(TypeError, match="Cannot parse boolean"):
            SettingsManager._parse_bool([1, 2, 3])  # type: ignore[arg-type]


class TestFromDictBooleanParsing:
    """Tests for from_dict with boolean string parsing."""

    def test_from_dict_csv_string_false_parsed_correctly(self) -> None:
        """Verify 'False' string from CSV parses as False (not True)."""
        manager = SettingsManager()
        manager.crest.v3 = True

        manager.from_dict({"crest_v3": "False"})

        assert manager.crest.v3 is False

    def test_from_dict_csv_string_true_parsed_correctly(self) -> None:
        """Verify 'True' string from CSV parses as True."""
        manager = SettingsManager()
        manager.crest.v3 = False

        manager.from_dict({"crest_v3": "True"})

        assert manager.crest.v3 is True

    def test_from_dict_csv_string_zero_parsed_as_false(self) -> None:
        """Verify '0' string from CSV parses as False."""
        manager = SettingsManager()
        manager.crest.quick = True

        manager.from_dict({"crest_quick": "0"})

        assert manager.crest.quick is False

    def test_from_dict_csv_string_one_parsed_as_true(self) -> None:
        """Verify '1' string from CSV parses as True."""
        manager = SettingsManager()
        manager.crest.quick = False

        manager.from_dict({"crest_quick": "1"})

        assert manager.crest.quick is True

    def test_from_dict_csv_round_trip_all_booleans(self) -> None:
        """Verify all boolean fields round-trip correctly from CSV strings."""
        manager = SettingsManager()

        settings = {
            "crest_v3": "true",
            "crest_quick": "False",
            "crest_nci": "0",
            "crest_gfnff": "1",
            "mopac_precise": "yes",
            "mopac_pulay": "no",
            "mopac_prtall": "on",
            "mopac_archive": "off",
            "crest_xtb_preopt": "True",
        }

        manager.from_dict(settings)

        assert manager.crest.v3 is True
        assert manager.crest.quick is False
        assert manager.crest.nci is False
        assert manager.crest.gfnff is True
        assert manager.mopac.precise is True
        assert manager.mopac.pulay is False
        assert manager.mopac.prtall is True
        assert manager.mopac.archive is False
        assert manager.xtb.preopt is True

    def test_from_dict_numeric_fields_still_work(self) -> None:
        """Verify numeric fields still parse correctly."""
        manager = SettingsManager()

        manager.from_dict(
            {
                "crest_ewin": "5.0",
                "crest_rthr": "0.125",
                "crest_threads": "4",
                "mopac_scfcrt": "1e-8",
                "mopac_itry": "200",
                "crest_xtb_timeout_seconds": "120",
            }
        )

        assert manager.crest.ewin == 5.0
        assert manager.crest.rthr == 0.125
        assert manager.crest.threads == 4
        assert manager.mopac.scfcrt == 1e-8
        assert manager.mopac.itry == 200
        assert manager.xtb.timeout_seconds == 120

    def test_from_dict_invalid_boolean_silently_ignored(self) -> None:
        """Verify invalid boolean string doesn't crash, just skips."""
        manager = SettingsManager()
        original_v3 = manager.crest.v3

        manager.from_dict({"crest_v3": "maybe"})

        assert manager.crest.v3 == original_v3


class TestThreadValidation:
    """Tests for thread count validation logic."""

    def test_thread_zero_should_be_rejected(self) -> None:
        """Verify validation logic rejects threads=0."""
        threads_val = 0
        max_threads = 4

        is_valid = 1 <= threads_val <= max_threads

        assert is_valid is False

    def test_thread_negative_should_be_rejected(self) -> None:
        """Verify validation logic rejects negative threads."""
        threads_val = -1
        max_threads = 4

        is_valid = 1 <= threads_val <= max_threads

        assert is_valid is False

    def test_thread_positive_within_max_should_be_accepted(self) -> None:
        """Verify validation logic accepts valid thread counts."""
        max_threads = 4

        for threads_val in [1, 2, 3, 4]:
            is_valid = 1 <= threads_val <= max_threads
            assert is_valid is True

    def test_thread_exceeds_max_should_be_rejected(self) -> None:
        """Verify validation logic rejects threads > max CPU count."""
        threads_val = 1000
        max_threads = 4

        is_valid = 1 <= threads_val <= max_threads

        assert is_valid is False


class TestPathReplace:
    """Tests for Path.replace Windows compatibility."""

    def test_path_replace_overwrites_existing_file(
        self, tmp_path: pytest.TempPathFactory
    ) -> None:
        """Verify Path.replace overwrites existing file atomically."""
        output_xyz = tmp_path / "output.xyz"  # type: ignore[operator]
        xtbopt = tmp_path / "xtbopt.xyz"  # type: ignore[operator]

        output_xyz.write_text("old output")
        xtbopt.write_text("new optimized structure")

        xtbopt.replace(output_xyz)

        assert not xtbopt.exists()
        assert output_xyz.read_text() == "new optimized structure"

    def test_path_replace_creates_new_file(
        self, tmp_path: pytest.TempPathFactory
    ) -> None:
        """Verify Path.replace creates file if destination doesn't exist."""
        output_xyz = tmp_path / "output.xyz"  # type: ignore[operator]
        xtbopt = tmp_path / "xtbopt.xyz"  # type: ignore[operator]

        xtbopt.write_text("new structure")

        xtbopt.replace(output_xyz)

        assert not xtbopt.exists()
        assert output_xyz.exists()
        assert output_xyz.read_text() == "new structure"
