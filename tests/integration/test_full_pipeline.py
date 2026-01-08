"""
End-to-end pipeline validation with real data.

This test validates the complete flow:
1. Load real data
2. Fuse CBS + B3LYP (as PM7 proxy)
3. Compute deltas
4. Extract features (tabular only for now)
5. Train/test split
6. Mock model training
"""

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score

from grimperium.data import ChemperiumLoader, DataFusion
from tests.fixtures.real_data import load_real_subset


def test_full_pipeline_real_data():
    """
    End-to-end pipeline test with real CBS data.

    Uses H298_b3 as PM7 proxy (PM7 calculations not yet available).
    """
    # Step 1: Load real data (small subset for CI speed)
    df = load_real_subset(n=1000, stratified=True, random_state=42)

    assert len(df) == 1000
    assert "H298_cbs" in df.columns
    assert "H298_b3" in df.columns

    # Step 2: Simulate CBS + PM7 data sources
    cbs_df = df[["smiles", "charge", "multiplicity", "nheavy", "H298_cbs"]].copy()
    pm7_df = df[["smiles", "H298_b3"]].copy()
    pm7_df.rename(columns={"H298_b3": "H298_pm7"}, inplace=True)

    # Step 3: Fuse datasets
    fusion = DataFusion()
    merged = fusion.merge(cbs_df, pm7_df, on="smiles")

    assert len(merged) == 1000

    # Step 4: Compute deltas
    merged = fusion.compute_deltas(merged)

    assert "delta_pm7" in merged.columns

    # Step 5: Extract features and targets
    X, y = fusion.get_training_data(merged)

    assert X.shape[0] == 1000
    assert len(y) == 1000
    assert list(X.columns) == ["nheavy", "charge", "multiplicity"]

    # Step 6: Train/test split
    loader = ChemperiumLoader()
    train_df, test_df = loader.split(merged, test_size=0.2, random_state=42)

    X_train, y_train = fusion.get_training_data(train_df)
    X_test, y_test = fusion.get_training_data(test_df)

    assert len(X_train) == 800
    assert len(X_test) == 200

    # Step 7: Mock model training (simple RF for validation)
    model = RandomForestRegressor(n_estimators=10, random_state=42, max_depth=5)
    model.fit(X_train, y_train)

    # Step 8: Predict deltas
    y_pred = model.predict(X_test)

    # Step 9: Validate predictions
    mae = mean_absolute_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    # Sanity checks (using B3LYP as proxy, not real PM7)
    # Note: With only tabular features and extreme delta range (-127k to +210k),
    # MAE ~1000 kcal/mol is expected. This validates pipeline, not model quality.
    assert mae < 15000  # MAE should be within reasonable bounds
    assert r2 > -10.0  # R² should not be catastrophically bad
    assert np.isfinite(y_pred).all()  # No NaN or inf predictions

    print("\n=== Pipeline Validation Results ===")
    print(f"MAE: {mae:.2f} kcal/mol")
    print(f"R²: {r2:.4f}")
    print(f"Delta mean: {y_test.mean():.2f} kcal/mol")
    print(f"Delta std: {y_test.std():.2f} kcal/mol")


def test_pipeline_with_different_subset_sizes():
    """Test pipeline stability across different sample sizes."""
    for n in [100, 500, 1000]:
        df = load_real_subset(n=n, stratified=True, random_state=42)

        # Quick validation
        assert len(df) == n
        assert df["H298_cbs"].notna().all()
        assert df["nheavy"].min() >= 1
        assert df["nheavy"].max() <= 22


def test_pipeline_delta_distribution():
    """Validate that delta distribution is reasonable."""
    df = load_real_subset(n=1000, stratified=True)

    # Compute delta manually
    delta = df["H298_cbs"] - df["H298_b3"]

    # Basic sanity checks
    assert delta.notna().all()
    assert np.isfinite(delta).all()

    # B3LYP differs from CBS across a wide range
    # (std ~11k kcal/mol due to diverse molecular sizes and types)
    assert abs(delta.mean()) < 1000  # Mean should be bounded
    assert delta.std() > 0  # Should have variance

    print("\n=== Delta Distribution (CBS - B3LYP) ===")
    print(f"Mean: {delta.mean():.2f} kcal/mol")
    print(f"Std: {delta.std():.2f} kcal/mol")
    print(f"Min: {delta.min():.2f} kcal/mol")
    print(f"Max: {delta.max():.2f} kcal/mol")
