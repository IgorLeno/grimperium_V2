"""Feature-frame assembly for CLI calculation workflows."""

from __future__ import annotations

from typing import Any

import pandas as pd

from grimperium.core.descriptors import extract_all_rdkit_descriptors
from grimperium.crest_pm7.mopac_descriptors import extract_mopac_descriptors
from grimperium.ml.features import FEATURE_COLUMNS


def build_pm7_delta_feature_frame(smiles: str, pm7_result: Any) -> pd.DataFrame:
    """Build one PM7+Delta feature row in the catalog column order."""
    h298_pm7 = pm7_result.most_stable_hof
    if h298_pm7 is None:
        raise ValueError("No stable conformer found - cannot extract HOF")

    selected_conformer = pm7_result.get_selected_conformer()
    if selected_conformer is None or selected_conformer.mopac_output_file is None:
        raise ValueError("No successful MOPAC output found")

    rdkit_descs = extract_all_rdkit_descriptors(smiles)
    mopac_descs = extract_mopac_descriptors(selected_conformer.mopac_output_file)

    features_dict: dict[str, float | int] = {
        "nheavy": pm7_result.nheavy or 0,
        "multiplicity": 1,
        "charge": 0,
        **rdkit_descs,
        "crest_conformers_generated": pm7_result.crest_conformers_generated,
        "num_conformers_selected": pm7_result.num_conformers_selected or 0,
        "crest_time_s": pm7_result.crest_time or 0.0,
        **mopac_descs,
        "H298_pm7": h298_pm7,
    }

    return pd.DataFrame([features_dict], columns=FEATURE_COLUMNS)
