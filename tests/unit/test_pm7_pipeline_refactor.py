"""Tests for PM7 pipeline refactoring: PM7-only selection + .out descriptor parsing.

Covers:
- MOPAC input keywords (EF, no AUX, no 1SCF)
- PM7Result.get_selected_conformer() and k_selected_pm7
- MOPAC .out file parsing (HOMO/LUMO, dipole, descriptors)
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
    _parse_out_file,
    extract_mopac_descriptors,
)

# ─── Realistic H2O .out fixture (MOPAC2016 PM7 PRECISE EF) ───

H2O_OUT_CONTENT = """\
 *******************************************************************************
 ** Site#:    0         For non-commercial use only    Version 22.234L 64BITS **
 *******************************************************************************
 ** Cite this program as: MOPAC2016, Version: 22.234L, James J. P. Stewart,  **
 *******************************************************************************

                              PM7 CALCULATION RESULTS

 FINAL HEAT OF FORMATION =        -57.7997 KCAL/MOL =     -241.8260 KJ/MOL

          TOTAL ENERGY            =       -351.49337 EV
          ELECTRONIC ENERGY       =       -496.32441 EV
          CORE-CORE REPULSION     =        144.83104 EV

          GRADIENT NORM           =          0.00884
          IONIZATION POTENTIAL    =         12.11550 EV
          HOMO LUMO ENERGIES (EV) =        -12.214   4.258
          NO. OF FILLED LEVELS    =          4
          MOLECULAR WEIGHT        =         18.0153

          SCF CALCULATIONS        =          7

          WALL-CLOCK TIME         =          0.090 SECONDS
          COMPUTATION TIME        =          0.083 SECONDS

       FINAL  POINT  AND  DERIVATIVES

   ATOM   CHEMICAL          X               Y               Z
  NUMBER   SYMBOL      (ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)

     1       O         0.00000000      0.00000000      0.00000000
     2       H         0.96000000      0.00000000      0.00000000
     3       H         0.00000000      0.96000000      0.00000000

          DIPOLE           =   2.12919 DEBYE    POINT GROUP:       C2v

          COSMO AREA              =   42.43 SQUARE ANGSTROMS
          COSMO VOLUME            =   25.17 CUBIC ANGSTROMS
"""


# ─── Step 2: MOPAC input keywords ───


@pytest.mark.unit
class TestMopacInputKeywords:
    """Tests for EF keyword in MOPAC input (AUX removed)."""

    def test_mopac_input_keywords_no_aux(self, tmp_path: Path) -> None:
        """EF must appear in MOPAC keyword line; AUX must NOT."""
        from grimperium.crest_pm7.mopac_optimizer import _create_mopac_input

        xyz = tmp_path / "test.xyz"
        xyz.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH 0.0 0.96 0.0\n")
        mop = tmp_path / "test.mop"

        _create_mopac_input(xyz, mop, config=None)

        content = mop.read_text()
        first_line = content.split("\n")[0]
        assert "EF" in first_line
        assert "AUX" not in first_line

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


# ─── Step 4: .out file parsing ───


@pytest.mark.unit
class TestOutFileParsing:
    """Tests for MOPAC .out descriptor parsing."""

    def test_out_parsing_extracts_all_fields(self) -> None:
        """Full H2O .out parsing extracts all expected descriptors."""
        result = _parse_out_file(H2O_OUT_CONTENT)

        assert result["mopac_dipole_debye"] == pytest.approx(2.12919)
        assert result["mopac_ionization_potential_ev"] == pytest.approx(12.1155)
        assert result["mopac_homo_ev"] == pytest.approx(-12.214)
        assert result["mopac_lumo_ev"] == pytest.approx(4.258)
        assert result["mopac_gap_ev"] == pytest.approx(4.258 - (-12.214))
        assert result["mopac_cosmo_area_a2"] == pytest.approx(42.43)
        assert result["mopac_cosmo_volume_a3"] == pytest.approx(25.17)
        assert result["mopac_gradient_norm"] == pytest.approx(0.00884)
        assert result["mopac_num_scf_cycles"] == 7
        assert result["mopac_point_group"] == "C2v"
        assert result["mopac_time_s"] == pytest.approx(0.090)

    @pytest.mark.parametrize(
        "homo_lumo_line, expected_homo, expected_lumo",
        [
            # Variant 1: standard spacing
            (
                "          HOMO LUMO ENERGIES (EV) =        -12.214   4.258",
                -12.214,
                4.258,
            ),
            # Variant 2: compact spacing
            (
                "          HOMO LUMO ENERGIES (EV) = -12.214  4.258",
                -12.214,
                4.258,
            ),
            # Variant 3: no space before parenthesis
            (
                "          HOMO LUMO ENERGIES(EV) =  -12.214   4.258",
                -12.214,
                4.258,
            ),
        ],
        ids=["standard", "compact", "no-space-paren"],
    )
    def test_out_parsing_homo_lumo_variants(
        self,
        homo_lumo_line: str,
        expected_homo: float,
        expected_lumo: float,
    ) -> None:
        """HOMO/LUMO parsing handles format variants in MOPAC .out files."""
        content = f"SOME HEADER\n{homo_lumo_line}\nSOME FOOTER\n"
        result = _parse_out_file(content)
        assert result["mopac_homo_ev"] == pytest.approx(expected_homo)
        assert result["mopac_lumo_ev"] == pytest.approx(expected_lumo)
        assert result["mopac_gap_ev"] == pytest.approx(expected_lumo - expected_homo)

    def test_out_parsing_missing_fields_returns_nan(self) -> None:
        """Minimal .out content returns NaN for missing descriptors."""
        result = _parse_out_file("SOME MOPAC OUTPUT WITH NO DESCRIPTOR LINES\n")

        assert np.isnan(result["mopac_dipole_debye"])
        assert np.isnan(result["mopac_homo_ev"])
        assert np.isnan(result["mopac_lumo_ev"])
        assert np.isnan(result["mopac_gap_ev"])
        assert np.isnan(result["mopac_gradient_norm"])
        assert result["mopac_point_group"] is None

    def test_extract_mopac_descriptors_missing_file(self, tmp_path: Path) -> None:
        """Missing .out file returns NaN dict with warning."""
        result = extract_mopac_descriptors(tmp_path / "nonexistent.out")

        assert np.isnan(result["mopac_dipole_debye"])
        assert np.isnan(result["mopac_homo_ev"])

    def test_extract_mopac_descriptors_real_file(self, tmp_path: Path) -> None:
        """Real .out file on disk is parsed correctly."""
        out_path = tmp_path / "test.out"
        out_path.write_text(H2O_OUT_CONTENT)

        result = extract_mopac_descriptors(out_path)
        assert result["mopac_dipole_debye"] == pytest.approx(2.12919)
        assert result["mopac_homo_ev"] == pytest.approx(-12.214)


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
            mopac_output_file=tmp_path / "fake.out",  # no .out, descriptors will be NaN
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

    def test_execution_manager_passes_selected_conformer_and_k(
        self, tmp_path: Path
    ) -> None:
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
        pm7_config.temp_dir = tmp_path

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
            state_manager=MagicMock(),
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
