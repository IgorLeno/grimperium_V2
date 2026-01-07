"""
Test fixtures for Grimperium delta-learning validation.

This module provides fixtures for testing with both synthetic and real data:
    - synthetic_data_1k: Fast synthetic data for unit tests
    - real_data_1k: Real CBS data + mock PM7 for integration tests
    - Mock PM7 generator with realistic error patterns

Mock PM7 Error Model:
    1. Base bias: PM7 overestimates by ~5 kcal/mol
    2. Size-dependent: error scales with sqrt(nheavy)
    3. Magnitude bias: 2% of |y_cbs|
    4. Gaussian noise: ~7 kcal/mol std

"""

import pytest
import numpy as np
from typing import Tuple


# ============================================================
# MOCK PM7 FUNCTIONS
# ============================================================


def create_realistic_mock_pm7(
    y_cbs: np.ndarray, X_basic: np.ndarray, seed: int = 42
) -> np.ndarray:
    """
    Generate H298_PM7 mock with realistic error patterns.

    Components:
    1. Base bias: PM7 overestimates ~5 kcal/mol
    2. Size-dependent: error scales with sqrt(nheavy)
    3. Magnitude bias: 2% of |y_cbs|
    4. Gaussian noise: ~7 kcal/mol std

    Args:
        y_cbs: True CBS energies (kcal/mol)
        X_basic: Basic features [nheavy, charge, mult]
        seed: Random seed for reproducibility

    Returns:
        Mock PM7 energies with realistic error patterns

    Example:
        >>> X_basic = np.array([[10, 0, 1], [20, -1, 2]])
        >>> y_cbs = np.array([-50.0, -100.0])
        >>> y_pm7 = create_realistic_mock_pm7(y_cbs, X_basic)
        >>> delta = y_cbs - y_pm7
        >>> assert np.std(delta) > 5  # Realistic spread

    """
    np.random.seed(seed)
    nheavy = X_basic[:, 0]

    # Component 1: Base bias (PM7 systematically overestimates)
    base_bias = -5.0

    # Component 2: Size-dependent error
    size_error = (1 + np.sqrt(np.maximum(nheavy, 1))) * np.random.normal(
        0, 1.5, len(y_cbs)
    )

    # Component 3: Magnitude-dependent bias
    magnitude_bias = (np.abs(y_cbs) / 100) * np.random.normal(0, 3, len(y_cbs))

    # Component 4: Random Gaussian noise
    gaussian_noise = np.random.normal(0, 7, len(y_cbs))

    # Total error
    total_error = base_bias + size_error + magnitude_bias + gaussian_noise

    return y_cbs - total_error


def create_enriched_features(X_basic: np.ndarray) -> np.ndarray:
    """
    Expand [nheavy, charge, mult] to 10 dimensions.

    Output features:
        1. nheavy (original)
        2. charge (original)
        3. mult (original)
        4. nheavy²
        5. sqrt(nheavy)
        6. nheavy × charge
        7. nheavy × mult
        8. charge²
        9. mult²
        10. bias term (1.0)

    Args:
        X_basic: Basic features of shape (n_samples, 3)

    Returns:
        Enriched features of shape (n_samples, 10)

    Example:
        >>> X_basic = np.array([[10, 0, 1], [20, -1, 2]])
        >>> X_enriched = create_enriched_features(X_basic)
        >>> assert X_enriched.shape == (2, 10)

    """
    nheavy = X_basic[:, 0:1]
    charge = X_basic[:, 1:2]
    mult = X_basic[:, 2:3]

    return np.hstack(
        [
            nheavy,  # 1. Original
            charge,  # 2. Original
            mult,  # 3. Original
            nheavy**2,  # 4. Polynomial
            np.sqrt(np.abs(nheavy) + 1),  # 5. Polynomial
            nheavy * charge,  # 6. Interaction
            nheavy * mult,  # 7. Interaction
            charge**2,  # 8. Polynomial
            mult**2,  # 9. Polynomial
            np.ones_like(nheavy),  # 10. Bias term
        ]
    )


# ============================================================
# FIXTURES
# ============================================================


@pytest.fixture
def synthetic_data_1k() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Synthetic data with enriched features (1000 samples).

    Use for: Fast unit tests (<1s), controlled validation

    Returns:
        Tuple of (X_enriched, y_cbs, y_pm7_mock, y_delta)
            - X_enriched: Features (1000, 10)
            - y_cbs: True CBS energies (1000,)
            - y_pm7_mock: Mock PM7 energies (1000,)
            - y_delta: CBS - PM7 deltas (1000,)

    Example:
        >>> def test_model(synthetic_data_1k):
        ...     X, y_cbs, y_pm7, y_delta = synthetic_data_1k
        ...     model.fit(X, y_delta)
        ...     assert model.is_fitted

    """
    np.random.seed(42)
    n = 1000

    # Generate basic features
    nheavy = np.random.randint(1, 25, size=n)
    charge = np.random.choice([-2, -1, 0, 1, 2], size=n)
    mult = np.random.choice([1, 2, 3], size=n)

    X_basic = np.column_stack([nheavy, charge, mult]).astype(float)

    # Generate synthetic CBS energies
    y_cbs = np.random.normal(loc=-40, scale=50, size=n)

    # Create enriched features
    X = create_enriched_features(X_basic)

    # Generate mock PM7 with realistic errors
    y_pm7 = create_realistic_mock_pm7(y_cbs, X_basic, seed=42)

    # Compute delta
    y_delta = y_cbs - y_pm7

    return X, y_cbs, y_pm7, y_delta


@pytest.fixture
def real_data_1k() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    DEPRECATED: Use fixtures from tests/experiments/conftest.py instead.

    This fixture uses UNFILTERED data which creates severe distribution shift
    between train/test splits (train mean=-594, test mean=+21). This causes
    Direct learning to fail catastrophically (RMSE=1008) which is an artifact,
    not meaningful comparison.

    Recommended alternatives:
    - real_data_1k_filtered: Realistic distribution for hypothesis validation
    - real_data_1k_extreme: Pathological distribution for stress testing

    Returns:
        Tuple of (X_enriched, y_cbs, y_pm7_mock, y_delta)

    WARNING: Contains outliers up to -325407 kcal/mol. Results may be
    misleading if used for hypothesis validation. See BATCH 3 documentation.
    """
    import warnings

    warnings.warn(
        "real_data_1k is DEPRECATED. Use real_data_1k_filtered from "
        "tests/experiments/conftest.py for hypothesis validation, or "
        "real_data_1k_extreme for stress testing.",
        DeprecationWarning,
        stacklevel=2,
    )

    from grimperium.data.loader import ChemperiumLoader

    # Load real CBS data (UNFILTERED - contains extreme outliers!)
    loader = ChemperiumLoader()
    df = loader.load_thermo_cbs_opt()

    # Sample 1000 rows (reproducible but includes outliers)
    df = df.sample(n=1000, random_state=42)

    # Extract basic features
    X_basic = df[["nheavy", "charge", "multiplicity"]].values.astype(float)
    y_cbs = df["H298_cbs"].values.astype(float)

    # Create enriched features
    X = create_enriched_features(X_basic)

    # Generate mock PM7 with realistic errors
    y_pm7 = create_realistic_mock_pm7(y_cbs, X_basic, seed=42)

    # Compute delta
    y_delta = y_cbs - y_pm7

    return X, y_cbs, y_pm7, y_delta
