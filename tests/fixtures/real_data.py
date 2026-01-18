"""
Real dataset fixtures for integration testing.

This module provides functions to load subsets of the real
thermo_cbs_chon.csv dataset (29,568 CHON molecules) for CI testing
without requiring the full dataset.

Example:
    >>> from tests.fixtures.real_data import load_real_subset
    >>> df = load_real_subset(n=1000, stratified=True)
    >>> print(f"Loaded {len(df)} real molecules")

"""

from pathlib import Path

import pandas as pd
from sklearn.model_selection import train_test_split

_DATASET_PATH = (
    Path(__file__).parent.parent.parent / "data" / "thermo_cbs_chon.csv"
)  # Primary CHON dataset


def load_real_subset(
    n: int = 1000,
    stratified: bool = True,
    random_state: int = 42,
    dataset_path: Path | None = None,
) -> pd.DataFrame:
    """
    Load a random subset of the real dataset.

    Args:
        n: Number of molecules to load
        stratified: Whether to stratify by nheavy (ensures diversity)
        random_state: Random seed for reproducibility
        dataset_path: Path to dataset (uses default if None)

    Returns:
        DataFrame with n molecules from real dataset

    Raises:
        FileNotFoundError: If dataset file not found

    Example:
        >>> df = load_real_subset(n=500, stratified=True)
        >>> assert len(df) == 500

    """
    path = dataset_path or _DATASET_PATH

    if not path.exists():
        raise FileNotFoundError(
            f"Real dataset not found at {path}. "
            f"Expected: data/thermo_cbs_chon.csv (29,568 CHON molecules)\n"
            f"Note: thermo_cbs_clean.csv and thermo_cbs_opt.csv are no longer used"
        )

    # Load full dataset
    df_full = pd.read_csv(path)

    if n >= len(df_full):
        return df_full

    # Sample subset
    if stratified:
        # Stratify by nheavy to ensure diverse molecular sizes
        # Handle case where some nheavy values have very few samples
        try:
            _, df_subset = train_test_split(
                df_full,
                test_size=n,
                random_state=random_state,
                stratify=df_full["nheavy"],
            )
        except ValueError:
            # If stratification fails (some classes too small), fall back to random sampling
            df_subset = df_full.sample(n=n, random_state=random_state)
    else:
        df_subset = df_full.sample(n=n, random_state=random_state)

    return df_subset.reset_index(drop=True)


def load_real_train_test_split(
    test_size: float = 0.2,
    max_samples: int | None = None,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load real dataset with train/test split.

    Args:
        test_size: Fraction for test set
        max_samples: Limit total samples (None for full dataset)
        random_state: Random seed

    Returns:
        Tuple of (train_df, test_df)

    Example:
        >>> train, test = load_real_train_test_split(test_size=0.2, max_samples=5000)

    """
    if not _DATASET_PATH.exists():
        raise FileNotFoundError(f"Real dataset not found at {_DATASET_PATH}")

    df = pd.read_csv(_DATASET_PATH)

    # Limit samples if requested
    if max_samples and max_samples < len(df):
        df = df.sample(n=max_samples, random_state=random_state)

    # Split stratified by nheavy
    # Handle case where some nheavy values have very few samples
    try:
        train_df, test_df = train_test_split(
            df,
            test_size=test_size,
            random_state=random_state,
            stratify=df["nheavy"],
        )
    except ValueError:
        # If stratification fails (some classes too small), fall back to non-stratified split
        train_df, test_df = train_test_split(
            df,
            test_size=test_size,
            random_state=random_state,
        )

    return train_df, test_df


def get_dataset_stats() -> dict:
    """
    Get statistics about the full real dataset.

    Returns:
        Dictionary with dataset statistics

    Example:
        >>> stats = get_dataset_stats()
        >>> print(f"Total molecules: {stats['n_molecules']}")

    """
    if not _DATASET_PATH.exists():
        raise FileNotFoundError(f"Real dataset not found at {_DATASET_PATH}")

    df = pd.read_csv(_DATASET_PATH)

    return {
        "n_molecules": len(df),
        "n_columns": len(df.columns),
        "columns": df.columns.tolist(),
        "h298_cbs_mean": float(df["H298_cbs"].mean()),
        "h298_cbs_std": float(df["H298_cbs"].std()),
        "h298_cbs_min": float(df["H298_cbs"].min()),
        "h298_cbs_max": float(df["H298_cbs"].max()),
        "nheavy_min": int(df["nheavy"].min()),
        "nheavy_max": int(df["nheavy"].max()),
        "nheavy_mean": float(df["nheavy"].mean()),
    }
