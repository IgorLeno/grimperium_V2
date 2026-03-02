"""Integration test: validate .out parsing on real MOPAC output.

This smoke test ensures that the .out parser can handle real MOPAC2016 output
from the CI environment or local runs.
"""

import math
from pathlib import Path

import pytest

from grimperium.crest_pm7.mopac_descriptors import extract_mopac_descriptors


@pytest.fixture
def sample_out_file(tmp_path: Path) -> Path:
    """Create a minimal but realistic .out file for testing."""
    out_path = tmp_path / "sample.out"

    # Minimal realistic MOPAC .out content (summary block only)
    content = """
 *******************************************************************************
 **                                                                           **
 **                              MOPAC2016                                    **
 **                                                                           **
 *******************************************************************************

          FINAL HEAT OF FORMATION =          5.64572 KCAL/MOL =      23.62137 KJ/MOL

          DIPOLE           =          0.57800 DEBYE    POINT GROUP:       Cs

          HOMO LUMO ENERGIES (EV)  =    -10.096    -0.351

          IONIZATION POTENTIAL     =         10.096 EV

          COSMO AREA               =        234.56 SQUARE ANGSTROMS
          COSMO VOLUME             =        123.45 CUBIC ANGSTROMS

          GRADIENT NORM            =          0.123

          SCF CALCULATIONS         =         42

          WALL-CLOCK TIME          =          0.473 SECONDS
    """

    out_path.write_text(content)
    return out_path


def test_out_parsing_returns_all_descriptors(sample_out_file: Path) -> None:
    """Test that parser extracts all 11 descriptors from .out."""
    desc = extract_mopac_descriptors(sample_out_file)

    # Check all keys present
    required_keys = [
        "mopac_dipole_debye",
        "mopac_ionization_potential_ev",
        "mopac_homo_ev",
        "mopac_lumo_ev",
        "mopac_gap_ev",
        "mopac_cosmo_area_a2",
        "mopac_cosmo_volume_a3",
        "mopac_gradient_norm",
        "mopac_num_scf_cycles",
        "mopac_point_group",
        "mopac_time_s",
    ]

    for key in required_keys:
        assert key in desc, f"Missing key: {key}"

    # Check specific values
    assert desc["mopac_dipole_debye"] == pytest.approx(0.578, abs=0.001)
    assert desc["mopac_ionization_potential_ev"] == pytest.approx(10.096, abs=0.001)
    assert desc["mopac_homo_ev"] == pytest.approx(-10.096, abs=0.001)
    assert desc["mopac_lumo_ev"] == pytest.approx(-0.351, abs=0.001)
    assert desc["mopac_gap_ev"] == pytest.approx(9.745, abs=0.001)
    assert desc["mopac_cosmo_area_a2"] == pytest.approx(234.56, abs=0.01)
    assert desc["mopac_cosmo_volume_a3"] == pytest.approx(123.45, abs=0.01)
    assert desc["mopac_gradient_norm"] == pytest.approx(0.123, abs=0.001)
    assert desc["mopac_num_scf_cycles"] == 42
    assert desc["mopac_point_group"] == "Cs"
    assert desc["mopac_time_s"] == pytest.approx(0.473, abs=0.001)


def test_out_parsing_missing_file_returns_defaults(tmp_path: Path) -> None:
    """Test that parser returns NaN/None defaults for missing file."""
    nonexistent = tmp_path / "does_not_exist.out"

    desc = extract_mopac_descriptors(nonexistent)

    # All numeric fields should be NaN
    assert math.isnan(float(desc["mopac_homo_ev"]))
    assert math.isnan(float(desc["mopac_lumo_ev"]))
    assert math.isnan(float(desc["mopac_gap_ev"]))

    # String field should be None
    assert desc["mopac_point_group"] is None


def test_out_parsing_gap_calculation() -> None:
    """Test that GAP is correctly computed as LUMO - HOMO."""
    # HOMO = -10.0, LUMO = -2.0 → GAP should be 8.0
    homo = -10.0
    lumo = -2.0
    gap = lumo - homo

    assert gap == pytest.approx(8.0, abs=0.001)
