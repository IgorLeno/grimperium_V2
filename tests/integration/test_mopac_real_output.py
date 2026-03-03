"""Integration test: full descriptor extraction from realistic MOPAC .out content.

Uses a representative snippet of MOPAC2016 output to verify end-to-end parsing.
"""

import numpy as np
import pytest
from pathlib import Path

from grimperium.crest_pm7.mopac_descriptors import extract_mopac_descriptors

# Representative MOPAC2016 .out content (trimmed to essential blocks)
_MOPAC_OUT_CONTENT = """\
 *******************************************************************************
 **                                                                           **
 **                              MOPAC v22.1.1                                **
 **                                                                           **
 *******************************************************************************

          FINAL HEAT OF FORMATION =       -119.38921 KCAL/MOL =    -499.52285 KJ/MOL

          TOTAL ENERGY            =      -1530.90442 EV
          ELECTRONIC ENERGY       =      -7419.70512 EV
          CORE-CORE REPULSION     =       5888.80070 EV

          GRADIENT NORM           =          0.98743 = 0.28520 PER ATOM

          IONIZATION POTENTIAL    =         10.12300 EV
          HOMO LUMO ENERGIES (EV) =        -10.096   -0.351

          MOLECULAR WEIGHT        =        110.112

          COSMO AREA              =        166.23 SQUARE ANGSTROMS
          COSMO VOLUME            =        161.87 CUBIC ANGSTROMS

          SCF CALCULATIONS        =         42

          WALL-CLOCK TIME         =          3.516 SECONDS

          POINT GROUP:    C2v

          DIPOLE           X         Y         Z       TOTAL
 POINT-CHG.    -0.123     0.456    -0.789     0.935
 HYBRID         0.001     0.002    -0.001     0.002
 SUM           -0.835     1.242     0.003     1.497
"""


class TestRealOutputExtraction:
    """End-to-end descriptor extraction from realistic MOPAC output."""

    @pytest.fixture()
    def mopac_out(self, tmp_path: Path) -> Path:
        """Write realistic output to a temp file."""
        out = tmp_path / "mol_001_conf_000.out"
        out.write_text(_MOPAC_OUT_CONTENT)
        return out

    def test_dipole_extracted(self, mopac_out: Path) -> None:
        """Dipole from table-format SUM line is extracted."""
        desc = extract_mopac_descriptors(mopac_out)
        assert desc["mopac_dipole_debye"] == pytest.approx(1.497)

    def test_homo_lumo_gap(self, mopac_out: Path) -> None:
        """HOMO, LUMO, and GAP are extracted correctly."""
        desc = extract_mopac_descriptors(mopac_out)
        assert desc["mopac_homo_ev"] == pytest.approx(-10.096)
        assert desc["mopac_lumo_ev"] == pytest.approx(-0.351)
        assert desc["mopac_gap_ev"] == pytest.approx(9.745)

    def test_ionization_potential(self, mopac_out: Path) -> None:
        desc = extract_mopac_descriptors(mopac_out)
        assert desc["mopac_ionization_potential_ev"] == pytest.approx(10.123)

    def test_cosmo(self, mopac_out: Path) -> None:
        desc = extract_mopac_descriptors(mopac_out)
        assert desc["mopac_cosmo_area_a2"] == pytest.approx(166.23)
        assert desc["mopac_cosmo_volume_a3"] == pytest.approx(161.87)

    def test_gradient_norm(self, mopac_out: Path) -> None:
        desc = extract_mopac_descriptors(mopac_out)
        assert desc["mopac_gradient_norm"] == pytest.approx(0.98743)

    def test_scf_cycles(self, mopac_out: Path) -> None:
        desc = extract_mopac_descriptors(mopac_out)
        assert desc["mopac_num_scf_cycles"] == 42

    def test_point_group(self, mopac_out: Path) -> None:
        desc = extract_mopac_descriptors(mopac_out)
        assert desc["mopac_point_group"] == "C2v"

    def test_wall_clock_time(self, mopac_out: Path) -> None:
        desc = extract_mopac_descriptors(mopac_out)
        assert desc["mopac_time_s"] == pytest.approx(3.516)

    def test_missing_file_returns_defaults(self, tmp_path: Path) -> None:
        """Non-existent file returns all NaN/None defaults."""
        desc = extract_mopac_descriptors(tmp_path / "nonexistent.out")
        assert np.isnan(desc["mopac_dipole_debye"])
        assert desc["mopac_point_group"] is None
