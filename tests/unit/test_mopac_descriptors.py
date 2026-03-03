"""Unit tests for MOPAC descriptor parsing.

Tests dipole parsing in both summary and table formats.
"""

import numpy as np
import pytest

from grimperium.crest_pm7.mopac_descriptors import _parse_out_file


class TestDipoleParsing:
    """Tests for dipole moment extraction from MOPAC .out content."""

    def test_summary_format(self) -> None:
        """Summary format 'DIPOLE = X.XXX DEBYE' is parsed correctly."""
        content = """
 IONIZATION POTENTIAL    =    10.123 EV
 DIPOLE           =     2.345 DEBYE
 HOMO LUMO ENERGIES (EV) =   -10.096   -0.351
"""
        result = _parse_out_file(content)
        assert result["mopac_dipole_debye"] == pytest.approx(2.345)

    def test_table_format(self) -> None:
        """Table format with SUM row and TOTAL column is parsed correctly."""
        content = """
          DIPOLE           X         Y         Z       TOTAL
 POINT-CHG.    -0.123     0.456    -0.789     0.935
 HYBRID         0.001     0.002    -0.001     0.002
 SUM           -0.835     1.242     0.003     1.497
"""
        result = _parse_out_file(content)
        assert result["mopac_dipole_debye"] == pytest.approx(1.497)

    def test_summary_format_takes_priority(self) -> None:
        """When both formats are present, summary format wins."""
        content = """
 DIPOLE           =     2.345 DEBYE

          DIPOLE           X         Y         Z       TOTAL
 SUM           -0.835     1.242     0.003     1.497
"""
        result = _parse_out_file(content)
        assert result["mopac_dipole_debye"] == pytest.approx(2.345)

    def test_no_dipole_returns_nan(self) -> None:
        """Missing dipole data returns NaN."""
        content = " IONIZATION POTENTIAL    =    10.123 EV\n"
        result = _parse_out_file(content)
        assert np.isnan(result["mopac_dipole_debye"])

    def test_table_format_negative_components(self) -> None:
        """Table format handles negative component values (TOTAL is always positive)."""
        content = """
          DIPOLE           X         Y         Z       TOTAL
 SUM          -3.210    -2.100    -0.500     3.889
"""
        result = _parse_out_file(content)
        assert result["mopac_dipole_debye"] == pytest.approx(3.889)


class TestHomoLumoParsing:
    """Tests for HOMO/LUMO energy extraction from MOPAC .out content."""

    def test_standard_format(self) -> None:
        """Standard 'HOMO LUMO ENERGIES (EV) = X Y' is parsed correctly."""
        content = """
 IONIZATION POTENTIAL    =    10.123 EV
 HOMO LUMO ENERGIES (EV) =   -8.696   -0.081
"""
        result = _parse_out_file(content)
        assert result["mopac_homo_ev"] == pytest.approx(-8.696)
        assert result["mopac_lumo_ev"] == pytest.approx(-0.081)
        assert result["mopac_gap_ev"] == pytest.approx(8.615)

    def test_various_decimal_places(self) -> None:
        """Values with many decimal places are parsed without truncation."""
        content = """
 HOMO LUMO ENERGIES (EV) =   -9.1234567   -0.7654321
"""
        result = _parse_out_file(content)
        assert result["mopac_homo_ev"] == pytest.approx(-9.1234567)
        assert result["mopac_lumo_ev"] == pytest.approx(-0.7654321)
        assert result["mopac_gap_ev"] == pytest.approx(8.3580246)

    def test_near_zero_values(self) -> None:
        """Near-zero HOMO/LUMO values are handled correctly."""
        content = """
 HOMO LUMO ENERGIES (EV) =   -0.001    0.002
"""
        result = _parse_out_file(content)
        assert result["mopac_homo_ev"] == pytest.approx(-0.001)
        assert result["mopac_lumo_ev"] == pytest.approx(0.002)
        assert result["mopac_gap_ev"] == pytest.approx(0.003)
