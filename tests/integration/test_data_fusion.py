"""
Integration tests for DataFusion with real dataset.

This module tests the DataFusion class using the real CBS/B3LYP dataset,
using H298_b3 as a proxy for PM7 since PM7 calculations are not yet available.

Tests:
    - Real data fusion with CBS and B3LYP
    - Task view creation with real data
"""

from tests.fixtures.real_data import load_real_subset


def test_datafusion_with_real_cbs_and_b3():
    """
    Test DataFusion using real CBS and B3LYP data.

    Note: We use H298_b3 as a proxy for semiempirical PM7
    since PM7 calculations are not yet available.
    """
    from grimperium.data import DataFusion

    # Load real subset
    df = load_real_subset(n=500, stratified=True)

    # Separate CBS and B3LYP columns (simulating two data sources)
    cbs_df = df[["smiles", "charge", "multiplicity", "nheavy", "H298_cbs"]].copy()
    b3_df = df[["smiles", "H298_b3"]].copy()
    b3_df.rename(columns={"H298_b3": "H298_pm7"}, inplace=True)

    # Fuse
    fusion = DataFusion()
    merged = fusion.merge(cbs_df, b3_df, on="smiles")

    assert len(merged) == 500
    assert "H298_cbs" in merged.columns
    assert "H298_pm7" in merged.columns

    # Compute deltas
    merged_with_delta = fusion.compute_deltas(merged)

    assert "delta_pm7" in merged_with_delta.columns
    assert merged_with_delta["delta_pm7"].notna().all()

    # Analyze deltas
    stats = fusion.analyze_deltas(merged_with_delta)

    # B3LYP typically higher than CBS for these molecules
    assert abs(stats["mean"]) < 200  # Reasonable delta range (B3LYP vs CBS)
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
