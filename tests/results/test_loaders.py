from __future__ import annotations

from pathlib import Path

import pandas as pd

from grimperium.results.loaders import load_analysis_dataframe


def test_load_analysis_dataframe_reads_legacy_wide_csv(tmp_path: Path) -> None:
    csv_path = tmp_path / "legacy.csv"
    pd.DataFrame(
        {
            "mol_id": ["m1", "m2"],
            "smiles": ["C", "CC"],
            "H298_cbs": [-10.0, -20.0],
            "H298_pm7": [-9.0, -18.0],
            "H298_predicted": [-10.5, -19.5],
        }
    ).to_csv(csv_path, index=False)

    df = load_analysis_dataframe(csv_path)

    assert list(df["mol_id"]) == ["m1", "m2"]
    assert list(df["H298_predicted"]) == [-10.5, -19.5]
    assert list(df["H298_cbs"]) == [-10.0, -20.0]


def test_load_analysis_dataframe_adapts_canonical_long_form_csv(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "calculation_results.csv"
    pd.DataFrame(
        [
            {
                "estimate_id": "m1-baseline",
                "run_id": "run_1",
                "molecule_smiles": "C",
                "molecule_name": "m1",
                "method_id": "crest_pm7",
                "method_version": "1.0.0",
                "property_id": "standard_enthalpy_of_formation",
                "role": "baseline",
                "hamiltonian": "PM7",
                "canonical_value": -9.0,
                "canonical_unit": "kcal/mol",
                "value_kcal_mol": -9.0,
            },
            {
                "estimate_id": "m1-final",
                "run_id": "run_1",
                "molecule_smiles": "C",
                "molecule_name": "m1",
                "method_id": "pm7_delta_learning",
                "method_version": "0.1.0",
                "property_id": "standard_enthalpy_of_formation",
                "role": "final",
                "hamiltonian": "",
                "canonical_value": -10.5,
                "canonical_unit": "kcal/mol",
                "value_kcal_mol": -10.5,
            },
        ]
    ).to_csv(csv_path, index=False)

    df = load_analysis_dataframe(csv_path)

    assert list(df["mol_id"]) == ["m1"]
    assert list(df["smiles"]) == ["C"]
    assert list(df["H298_pm7"]) == [-9.0]
    assert list(df["H298_predicted"]) == [-10.5]
