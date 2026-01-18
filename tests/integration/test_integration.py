"""
Integration tests for real dataset fixtures.

Tests cover:
    - Real data loading and subsetting
    - Stratified sampling
    - Train/test splits
    - Dataset statistics
"""

from tests.fixtures.real_data import (
    get_dataset_stats,
    load_real_subset,
    load_real_train_test_split,
)


def test_real_fixture_availability():
    """Test that real data fixture loads correctly."""
    df = load_real_subset(n=100)

    assert len(df) == 100
    assert "smiles" in df.columns
    assert "H298_cbs" in df.columns
    assert df["H298_cbs"].notna().all()


def test_real_subset_stratification():
    """Test that stratified sampling preserves nheavy distribution."""
    df = load_real_subset(n=1000, stratified=True)

    # Should have diverse molecular sizes
    nheavy_unique = df["nheavy"].nunique()
    assert nheavy_unique >= 10, "Should sample diverse molecular sizes"


def test_real_train_test_split():
    """Test train/test split on real data."""
    train, test = load_real_train_test_split(test_size=0.2, max_samples=1000)

    assert len(train) == 800
    assert len(test) == 200
    assert set(train.columns) == set(test.columns)

    # No overlap
    train_smiles = set(train["smiles"])
    test_smiles = set(test["smiles"])
    assert len(train_smiles & test_smiles) == 0


def test_dataset_stats():
    """Test dataset statistics calculation."""
    stats = get_dataset_stats()

    assert stats["n_molecules"] == 29568
    assert stats["n_columns"] == 7
    assert "H298_cbs" in stats["columns"]
    assert stats["nheavy_min"] >= 1
    assert stats["nheavy_max"] <= 22
