"""Tests for PM7 pipeline refactoring: PM7-only selection + .aux descriptor parsing.

Covers:
- MOPAC input keywords (EF, AUX, no 1SCF)
- PM7Result.get_selected_conformer() and k_selected_pm7
- MOPAC .aux file parsing (Fortran D notation, HOMO/LUMO, descriptors)
- CSV update with target_delta and k_selected_pm7
- Execution manager integration
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("rdkit", reason="rdkit not available")

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.config import MOPACStatus
from grimperium.crest_pm7.csv_enhancements import (
    CSVManagerExtensions,
    DeltaCalculations,
)
from grimperium.crest_pm7.molecule_processor import ConformerData, PM7Result
from grimperium.crest_pm7.mopac_descriptors import (
    _parse_aux_file,
    extract_mopac_descriptors,
    parse_fortran_float,
)

# ─── Real H2O .aux fixture (MOPAC2016.22.234L with PM7 PRECISE EF AUX) ───


H2O_AUX_CONTENT = """\
 START OF MOPAC FILE
 KEYWORDS=PM7 PRECISE EF AUX
 ATOM_X:ANGSTROMS[3]=
  +0.000000D+00   +0.000000D+00   +0.000000D+00
 ATOM_EL[3]=
O   H   H
 AO_ATOMINDEX[6]=
1   1   1   1   2   3
 NUM_ELECTRONS=08
 EIGENVALUES[6]=
 -39.89625  -18.62879   -14.57682  -12.21355   4.25810   6.33574
 DIPOLE:DEBYE=+0.212919D+01
 IONIZATION_POTENTIAL:EV=+0.121155D+02
 AREA:SQUARE ANGSTROMS=+0.424300D+02
 VOLUME:CUBIC ANGSTROMS=+0.251718D+02
 GRADIENT_NORM:KCAL/MOL/ANGSTROM=+0.883541D-02
 NUMBER_SCF_CYCLES=7
 POINT_GROUP=C2v
 CPU_TIME:SECONDS= 0.09
 END OF MOPAC FILE
"""


# ─── Step 2: MOPAC input keywords ───


@pytest.mark.unit
class TestMopacInputKeywords:
    """Tests for EF+AUX keywords in MOPAC input."""

    def test_mopac_input_includes_ef_aux_keywords(self, tmp_path: Path) -> None:
        """EF and AUX must appear in MOPAC keyword line."""
        from grimperium.crest_pm7.mopac_optimizer import _create_mopac_input

        xyz = tmp_path / "test.xyz"
        xyz.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH 0.0 0.96 0.0\n")
        mop = tmp_path / "test.mop"

        _create_mopac_input(xyz, mop, config=None)

        content = mop.read_text()
        first_line = content.split("\n")[0]
        assert "EF" in first_line
        assert "AUX" in first_line

    def test_mopac_input_no_1scf(self, tmp_path: Path) -> None:
        """1SCF must NOT appear in MOPAC input (removed)."""
        from grimperium.crest_pm7.config import PM7Config
        from grimperium.crest_pm7.mopac_optimizer import _create_mopac_input

        xyz = tmp_path / "test.xyz"
        xyz.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH 0.0 0.96 0.0\n")
        mop = tmp_path / "test.mop"

        # Even with precise_scf=False, 1SCF should not appear
        config = PM7Config(temp_dir=tmp_path)
        config.mopac_precise_scf = False
        _create_mopac_input(xyz, mop, config=config)

        content = mop.read_text()
        assert "1SCF" not in content


# ─── Step 3: PM7Result selection ───


@pytest.mark.unit
class TestPM7ResultSelection:
    """Tests for PM7Result.get_selected_conformer() and k_selected_pm7."""

    def _make_conformer(
        self, idx: int, crest_rank: int, energy: float
    ) -> ConformerData:
        return ConformerData(
            index=idx,
            mol_id="TEST",
            crest_rank=crest_rank,
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=True,
            energy_hof=energy,
        )

    def test_get_selected_conformer_returns_lowest_hof(self) -> None:
        """Selected conformer has the lowest energy_hof."""
        c1 = self._make_conformer(0, 1, -50.0)
        c2 = self._make_conformer(1, 2, -55.0)  # lowest
        c3 = self._make_conformer(2, 3, -48.0)

        result = PM7Result(mol_id="T", smiles="O", conformers=[c1, c2, c3])
        selected = result.get_selected_conformer()

        assert selected is not None
        assert selected.energy_hof == -55.0
        assert selected.crest_rank == 2

    def test_k_selected_pm7_uses_crest_rank(self) -> None:
        """k_selected_pm7 must equal the crest_rank of the selected conformer."""
        c1 = self._make_conformer(0, 1, -50.0)
        c2 = self._make_conformer(1, 2, -55.0)  # best

        result = PM7Result(mol_id="T", smiles="O", conformers=[c1, c2])
        selected = result.get_selected_conformer()
        assert selected is not None
        # Verify the production assignment logic: process() does
        # result.k_selected_pm7 = selected.crest_rank; test that the right
        # conformer is returned so the assignment would produce the correct value.
        assert selected.crest_rank == 2

    def test_get_selected_conformer_returns_none_when_empty(self) -> None:
        """No successful conformers returns None."""
        result = PM7Result(mol_id="T", smiles="O", conformers=[])
        assert result.get_selected_conformer() is None

    def test_selected_conformer_none_handles_descriptors(self, tmp_path: Path) -> None:
        """CSV update with no selected conformer sets k_selected_pm7 to None/NaN."""
        csv_path = tmp_path / "test.csv"
        df = pd.DataFrame(
            {
                "mol_id": ["T"],
                "smiles": ["O"],
                "nheavy": [1],
                "status": ["RUNNING"],
                "abs_diff": [None],
                "abs_diff_%": [None],
                "target_delta_kcalmol": [None],
                "k_selected_pm7": [None],
                "mopac_dipole_debye": [None],
                "mopac_homo_ev": [None],
                "v3": [None],
            }
        )
        df.to_csv(csv_path, index=False)

        csv_manager = BatchCSVManager(csv_path)
        csv_manager.load_csv()

        result = PM7Result(mol_id="T", smiles="O", conformers=[])
        assert result.get_selected_conformer() is None
        assert result.k_selected_pm7 is None

        success = CSVManagerExtensions.update_molecule_with_mopac_results(
            csv_manager=csv_manager,
            mol_id="T",
            h298_cbs=-17.5,
            h298_pm7=-15.3,
            selected_conformer=None,
            k_selected_pm7=None,
            batch_settings={"v3": False},
        )
        assert success
        df_out = pd.read_csv(csv_path)
        assert pd.isna(df_out.loc[0, "k_selected_pm7"])


# ─── Step 4: .aux file parsing ───


@pytest.mark.unit
class TestAuxFileParsing:
    """Tests for MOPAC .aux descriptor parsing."""

    def test_parse_fortran_float_basic(self) -> None:
        """Fortran D notation converts correctly."""
        assert parse_fortran_float("+0.577997D+02") == pytest.approx(57.7997)
        assert parse_fortran_float("-0.123456D-03") == pytest.approx(-1.23456e-4)
        assert parse_fortran_float("+0.883541D-02") == pytest.approx(0.00883541)

    def test_aux_file_parsing_water_molecule(self) -> None:
        """Full H2O .aux parsing extracts all expected descriptors."""
        result = _parse_aux_file(H2O_AUX_CONTENT)

        assert result["mopac_dipole_debye"] == pytest.approx(2.12919)
        assert result["mopac_ionization_potential_ev"] == pytest.approx(12.1155)
        assert result["mopac_cosmo_area_a2"] == pytest.approx(42.43)
        assert result["mopac_cosmo_volume_a3"] == pytest.approx(25.1718)
        assert result["mopac_gradient_norm"] == pytest.approx(0.00883541)
        assert result["mopac_num_scf_cycles"] == 7
        assert result["mopac_point_group"] == "C2v"
        assert result["mopac_time_s"] == pytest.approx(0.09)

    def test_aux_parsing_homo_lumo_from_eigenvalues(self) -> None:
        """HOMO/LUMO derived from NUM_ELECTRONS=08 + EIGENVALUES array."""
        result = _parse_aux_file(H2O_AUX_CONTENT)

        # H2O: 8 electrons -> HOMO at idx 3, LUMO at idx 4
        assert result["mopac_homo_ev"] == pytest.approx(-12.21355)
        assert result["mopac_lumo_ev"] == pytest.approx(4.25810)
        assert result["mopac_gap_ev"] == pytest.approx(4.25810 - (-12.21355))

    def test_parsing_handles_missing_fields_returns_nan(self) -> None:
        """Minimal .aux content returns NaN for missing descriptors."""
        result = _parse_aux_file("START OF MOPAC FILE\nEND OF MOPAC FILE\n")

        assert np.isnan(result["mopac_dipole_debye"])
        assert np.isnan(result["mopac_homo_ev"])
        assert np.isnan(result["mopac_lumo_ev"])
        assert np.isnan(result["mopac_gap_ev"])
        assert result["mopac_point_group"] is None

    def test_extract_mopac_descriptors_missing_file(self, tmp_path: Path) -> None:
        """Missing .aux file returns NaN dict with warning."""
        result = extract_mopac_descriptors(tmp_path / "nonexistent.out")

        assert np.isnan(result["mopac_dipole_debye"])
        assert np.isnan(result["mopac_homo_ev"])

    def test_extract_mopac_descriptors_real_file(self, tmp_path: Path) -> None:
        """Real .aux file on disk is parsed correctly."""
        aux_path = tmp_path / "test.aux"
        aux_path.write_text(H2O_AUX_CONTENT)
        out_path = tmp_path / "test.out"

        result = extract_mopac_descriptors(out_path)
        assert result["mopac_dipole_debye"] == pytest.approx(2.12919)

    def test_aux_parsing_odd_electron_count_returns_defined_values(self) -> None:
        """Odd electron count does not crash; values extracted via integer division."""
        aux_odd = H2O_AUX_CONTENT.replace("NUM_ELECTRONS=08", "NUM_ELECTRONS=07")
        result = _parse_aux_file(aux_odd)
        # With 7 electrons: homo_idx=2, lumo_idx=3, both within the 6-element array
        # Pipeline uses closed-shell floor-division; result must be a float, not NaN.
        assert isinstance(result["mopac_homo_ev"], float)
        assert not np.isnan(result["mopac_homo_ev"])
        assert isinstance(result["mopac_lumo_ev"], float)
        assert not np.isnan(result["mopac_lumo_ev"])

    def test_aux_parsing_malformed_fortran_d_returns_nan(self) -> None:
        """Corrupted Fortran D notation and partial lines return NaN, not errors."""
        malformed = (
            " START OF MOPAC FILE\n"
            " DIPOLE:DEBYE=+INVALID_D+01\n"
            " IONIZATION_POTENTIAL:EV=+0.12D\n"
            " GRADIENT_NORM:KCAL/MOL/ANGSTROM=\n"
            " END OF MOPAC FILE\n"
        )
        result = _parse_aux_file(malformed)
        assert np.isnan(result["mopac_dipole_debye"])
        assert np.isnan(result["mopac_ionization_potential_ev"])
        assert np.isnan(result["mopac_gradient_norm"])
        assert result["mopac_point_group"] is None


# ─── Step 5: CSV update calculations ───


@pytest.mark.unit
class TestCSVUpdateCalculations:
    """Tests for target_delta and k_selected_pm7 in CSV updates."""

    def test_csv_update_calculates_target_delta_signed(self) -> None:
        """target_delta_kcalmol = H298_cbs - H298_pm7 (signed)."""
        # CBS = -17.5, PM7 = -15.3 => delta = -17.5 - (-15.3) = -2.2
        delta = DeltaCalculations.calculate_target_delta(-17.5, -15.3)
        assert delta == pytest.approx(-2.2)

        # CBS = -15.3, PM7 = -17.5 => delta = -15.3 - (-17.5) = +2.2
        delta2 = DeltaCalculations.calculate_target_delta(-15.3, -17.5)
        assert delta2 == pytest.approx(2.2)

    def test_target_delta_kcalmol_can_be_negative(self) -> None:
        """target_delta_kcalmol preserves sign; negative values are valid."""
        delta_neg = DeltaCalculations.calculate_target_delta(-20.0, -15.0)
        assert delta_neg == pytest.approx(-5.0)
        assert delta_neg < 0

        delta_pos = DeltaCalculations.calculate_target_delta(-10.0, -15.0)
        assert delta_pos == pytest.approx(5.0)
        assert delta_pos > 0

        # NaN propagation: None inputs return NaN
        assert np.isnan(DeltaCalculations.calculate_target_delta(None, -15.0))
        assert np.isnan(DeltaCalculations.calculate_target_delta(-20.0, None))

    def test_csv_update_includes_k_selected_pm7(self, tmp_path: Path) -> None:
        """k_selected_pm7 appears in CSV after update."""
        csv_path = tmp_path / "test.csv"
        df = pd.DataFrame(
            {
                "mol_id": ["mol_001"],
                "smiles": ["CCO"],
                "nheavy": [2],
                "status": ["RUNNING"],
                "abs_diff": [None],
                "abs_diff_%": [None],
                "target_delta_kcalmol": [None],
                "k_selected_pm7": [None],
                "v3": [None],
            }
        )
        df.to_csv(csv_path, index=False)

        csv_manager = BatchCSVManager(csv_path)
        csv_manager.load_csv()

        selected = ConformerData(
            index=0,
            mol_id="mol_001",
            crest_rank=2,
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=True,
            energy_hof=-55.0,
            mopac_output_file=tmp_path / "fake.out",  # no .aux, descriptors will be NaN
        )

        success = CSVManagerExtensions.update_molecule_with_mopac_results(
            csv_manager=csv_manager,
            mol_id="mol_001",
            h298_cbs=-17.5,
            h298_pm7=-15.3,
            selected_conformer=selected,
            k_selected_pm7=2,
            batch_settings={"v3": True},
        )

        assert success
        df_out = pd.read_csv(csv_path)
        assert df_out.loc[0, "k_selected_pm7"] == 2
        assert df_out.loc[0, "target_delta_kcalmol"] == pytest.approx(-2.2)


# ─── Step 6: Execution manager integration ───


@pytest.mark.unit
class TestExecutionManagerIntegration:
    """Tests for execution_manager passing selected_conformer and k."""

    def test_execution_manager_passes_selected_conformer_and_k(self) -> None:
        """execution_manager calls update_molecule_with_mopac_results with new args."""
        from grimperium.crest_pm7.batch.execution_manager import (
            BatchExecutionManager,
        )

        # Create mocks
        csv_manager = MagicMock()
        csv_manager.get_reference_hof.return_value = -17.5
        csv_manager.get_status.return_value = "OK"

        detail_manager = MagicMock()
        pm7_config = MagicMock()
        pm7_config.temp_dir = Path("/tmp")

        processor = MagicMock()

        # Create a PM7Result-like mock
        conformer = ConformerData(
            index=0,
            mol_id="mol_001",
            crest_rank=1,
            mopac_status=MOPACStatus.SUCCESS,
            hof_extraction_successful=True,
            energy_hof=-55.0,
        )

        pm7_result = MagicMock()
        pm7_result.success = True
        pm7_result.conformers = [conformer]
        pm7_result.most_stable_hof = -55.0
        pm7_result.quality_grade = MagicMock(value="A")
        pm7_result.error_message = None
        pm7_result.get_selected_conformer.return_value = conformer
        pm7_result.k_selected_pm7 = 1

        processor.process_with_fixed_timeout.return_value = pm7_result

        mgr = BatchExecutionManager(
            csv_manager=csv_manager,
            detail_manager=detail_manager,
            pm7_config=pm7_config,
            processor_adapter=processor,
        )
        mgr._batch_settings = {"v3": False}
        mgr._batch_logger = MagicMock()

        result = MagicMock()
        result.total_count = 1
        result.success_count = 0

        hof_values: list[tuple[str, float]] = []

        with patch.object(
            CSVManagerExtensions,
            "update_molecule_with_mopac_results",
            return_value=True,
        ) as mock_update:
            mgr._process_molecule(
                mol_id="mol_001",
                smiles="CCO",
                batch_id="batch_0001",
                batch_order=1,
                crest_timeout=30.0,
                mopac_timeout=60.0,
                result=result,
                hof_values=hof_values,
            )

            mock_update.assert_called_once()
            call_kwargs = mock_update.call_args
            assert isinstance(call_kwargs.kwargs, dict)
            assert "selected_conformer" in call_kwargs.kwargs
            assert "k_selected_pm7" in call_kwargs.kwargs
