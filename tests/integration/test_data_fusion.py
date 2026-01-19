"""
Integration tests for DataFusion with real dataset.

This module tests the DataFusion class using the real CBS dataset
with synthetic PM7 values (PM7 calculations not yet available).

Tests:
    - Real data fusion with CBS and synthetic PM7
    - Task view creation with real data
"""

import numpy as np

from tests.fixtures.real_data import load_real_subset


import numpy as np
import pandas as pd

from tests.fixtures.real_data import load_real_subset


def _create_synthetic_pm7(df: "pd.DataFrame", random_state: int = 42) -> "pd.DataFrame":
    """Create synthetic PM7 values with realistic noise for testing."""
    rng = np.random.default_rng(random_state)
    pm7_df = df[["smiles"]].copy()
    pm7_df["H298_pm7"] = df["H298_cbs"] + rng.normal(5.0, 10.0, len(df))
    return pm7_df


def test_datafusion_with_real_cbs_and_synthetic_pm7():
    """
    Test DataFusion using real CBS and synthetic PM7 data.

    Note: We use synthetic PM7 values since PM7 calculations
    are not yet available.
    """
    from grimperium.data import DataFusion

    # Load real subset
    df = load_real_subset(n=500, stratified=True)

    # Separate CBS and create synthetic PM7 (simulating two data sources)
    cbs_df = df[["smiles", "charge", "multiplicity", "nheavy", "H298_cbs"]].copy()
    pm7_df = _create_synthetic_pm7(df, random_state=42)

    # Fuse
    fusion = DataFusion()
    merged = fusion.merge(cbs_df, pm7_df, on="smiles")

    assert len(merged) == 500
    assert "H298_cbs" in merged.columns
    assert "H298_pm7" in merged.columns

    # Compute deltas
    merged_with_delta = fusion.compute_deltas(merged)

    assert "delta_pm7" in merged_with_delta.columns
    assert merged_with_delta["delta_pm7"].notna().all()

    # Analyze deltas
    stats = fusion.analyze_deltas(merged_with_delta)

    # Synthetic PM7 has mean offset ~5, std ~10
    assert abs(stats["mean"]) < 20  # Reasonable delta range
    assert stats["std"] > 0


def test_real_data_task_views():
    """Test task view creation with real data."""
    from grimperium.data import DataFusion

    df = load_real_subset(n=200, stratified=True)

    fusion = DataFusion()

    # Enthalpy task
    X, y = fusion.select_task_view(df, task="enthalpy")
    assert len(X) == 200
    assert len(y) == 200
    assert "nheavy" in X.columns
    assert y.name == "H298_cbs"
