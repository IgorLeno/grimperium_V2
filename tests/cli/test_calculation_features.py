"""Tests for reusable PM7+Delta feature assembly."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from grimperium.ml.features import FEATURE_COLUMNS


class _SelectedConformer:
    mopac_output_file = Path("selected.out")


class _PM7Result:
    nheavy = 3
    most_stable_hof = -65.5
    crest_conformers_generated = 7
    num_conformers_selected = 3
    crest_time = 12.0

    def get_selected_conformer(self) -> _SelectedConformer:
        return _SelectedConformer()


def test_build_pm7_delta_feature_frame_uses_catalog_column_order(
    monkeypatch: Any,
) -> None:
    from grimperium.cli import calculation_features

    rdkit_descs = {
        "rdkit_nrotbonds": 1,
        "rdkit_tpsa": 20.0,
        "rdkit_num_rings": 0,
        "rdkit_fsp3": 1.0,
        "rdkit_mol_weight": 46.0,
        "rdkit_hbond_donors": 1,
        "rdkit_hbond_acceptors": 1,
        "rdkit_nC": 2,
        "rdkit_nH": 6,
        "rdkit_nO": 1,
        "rdkit_nN": 0,
        "rdkit_bonds_single": 8,
        "rdkit_bonds_double": 0,
        "rdkit_bonds_triple": 0,
        "rdkit_bonds_aromatic": 0,
    }
    mopac_descs = {
        "mopac_homo_ev": -10.0,
        "mopac_lumo_ev": 1.0,
        "mopac_gap_ev": 11.0,
        "mopac_dipole_debye": 1.7,
        "mopac_ionization_potential_ev": 10.5,
        "mopac_cosmo_area_a2": 80.0,
        "mopac_cosmo_volume_a3": 70.0,
        "mopac_gradient_norm": 0.01,
        "mopac_num_scf_cycles": 12,
    }
    monkeypatch.setattr(
        calculation_features,
        "extract_all_rdkit_descriptors",
        lambda _smiles: rdkit_descs,
    )
    monkeypatch.setattr(
        calculation_features,
        "extract_mopac_descriptors",
        lambda _path: mopac_descs,
    )

    frame = calculation_features.build_pm7_delta_feature_frame("CCO", _PM7Result())

    assert list(frame.columns) == FEATURE_COLUMNS
    assert frame.shape == (1, len(FEATURE_COLUMNS))
    assert frame.loc[0, "H298_pm7"] == -65.5
    assert frame.loc[0, "crest_conformers_generated"] == 7
