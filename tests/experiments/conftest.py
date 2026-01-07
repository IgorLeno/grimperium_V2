"""
Fixtures for BATCH 3 hypothesis validation experiments.

DESIGN DECISIONS:
-----------------
1. `real_data_1k_filtered`: Realistic distribution with filter
   [-1000, +1000] kcal/mol
   - Removes extreme outliers (0.9% of data)
   - Produces std ~70 kcal/mol (not 7230!)
   - Safe for Week 2/3 reuse

2. `real_data_1k_extreme`: Pathological distribution for stress testing
   - Keeps all outliers (H298_cbs up to -325407)
   - Expected to cause distribution shift between train/test
   - Use ONLY for stress tests, not hypothesis validation

WHY FILTER?
-----------
Original dataset has H298_cbs range [-325407, +164949] with std=7230.
Random sampling creates severe distribution shift
(train mean=-594, test mean=+21).
Direct learning fails catastrophically (RMSE=1008), but this is an artifact,
not evidence of delta-learning superiority.

For hypothesis validation, we need comparable regimes where both models
have a fair chance. The filter achieves this while keeping 99.1% of data.

MOCK PM7:
---------
Real PM7 data not available. Mock PM7 simulates realistic error patterns:
- Base bias: PM7 overestimates by ~5 kcal/mol
- Size-dependent: error scales with sqrt(nheavy)
- Magnitude bias: 2% of |y_cbs|
- Gaussian noise: ~7 kcal/mol std

This mock is deterministic (seed=42) for reproducibility.
"""

import pytest
import numpy as np
from typing import Optional, Tuple


def create_realistic_mock_pm7(
    y_cbs: np.ndarray,
    X_basic: np.ndarray,
    seed: Optional[int] = 42,
    magnitude_bias_std: float = 0.5,
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
        seed: Random seed for reproducibility. If None, uses global RNG state.
        magnitude_bias_std: Standard deviation for magnitude-dependent bias.

    Returns:
        Mock PM7 energies with realistic error patterns
    """
    if seed is not None:
        np.random.seed(seed)
    nheavy = X_basic[:, 0]

    # Component 1: Base bias (PM7 systematically overestimates)
    base_bias = -5.0

    # Component 2: Size-dependent error
    size_error = (1 + np.sqrt(np.maximum(nheavy, 1))) * np.random.normal(
        0, 1.5, len(y_cbs)
    )

    # Component 3: Magnitude-dependent bias (2% de |y_cbs| por padrão).
    # Nota: em dados filtrados, esse termo tende a ser pequeno porque
    # |y_cbs| < 1000.
    magnitude_bias = (np.abs(y_cbs) / 50) * np.random.normal(
        0, magnitude_bias_std, len(y_cbs)
    )

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
        4. nheavy^2
        5. sqrt(nheavy)
        6. nheavy * charge
        7. nheavy * mult
        8. charge^2
        9. mult^2
        10. bias term (1.0)

    Args:
        X_basic: Basic features of shape (n_samples, 3)

    Returns:
        Enriched features of shape (n_samples, 10)
    """
    nheavy = X_basic[:, 0:1]
    charge = X_basic[:, 1:2]
    mult = X_basic[:, 2:3]

    return np.hstack([
        nheavy,                              # 1. Original
        charge,                              # 2. Original
        mult,                                # 3. Original
        nheavy ** 2,                         # 4. Polynomial
        np.sqrt(np.abs(nheavy) + 1),         # 5. Polynomial
        nheavy * charge,                     # 6. Interaction
        nheavy * mult,                       # 7. Interaction
        charge ** 2,                         # 8. Polynomial
        mult ** 2,                           # 9. Polynomial
        np.ones_like(nheavy)                 # 10. Bias term
    ])


# ============================================================
# FIXTURES
# ============================================================


@pytest.fixture
def real_data_1k_filtered() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load 1000 molecules with REALISTIC distribution (filtered).

    Filter: -1000 <= H298_cbs <= +1000 kcal/mol
    Impact: Removes 0.9% outliers, keeps 99.1% of data

    Justification:
    - Most organic molecules have enthalpies in this range
    - Outliers >|1000| are extreme cases (highly exo/endothermic reactions)
    - Representative of main Grimperium use case

    Statistics after filter:
    - mean: ~-34 kcal/mol
    - std: ~70 kcal/mol (NOT 7230!)
    - No severe distribution shift between train/test

    Returns:
        Tuple of (X_enriched, y_cbs, y_pm7_mock)
        - X_enriched: Features (1000, 10)
        - y_cbs: True CBS energies (1000,)
        - y_pm7_mock: Mock PM7 energies (1000,)

    Example:
        >>> def test_hypothesis(real_data_1k_filtered):
        ...     X, y_cbs, y_pm7 = real_data_1k_filtered
        ...     # Train/test split will have similar distributions
    """
    from grimperium.data.loader import ChemperiumLoader

    # Load all data
    loader = ChemperiumLoader()
    df = loader.load_thermo_cbs_opt()

    # Log original statistics
    print(f"\n[FIXTURE] Original data: {len(df)} molecules")
    h_mean = df["H298_cbs"].mean()
    h_std = df["H298_cbs"].std()
    h_min = df["H298_cbs"].min()
    h_max = df["H298_cbs"].max()
    print(f"  H298_cbs: mean={h_mean:.1f}, std={h_std:.1f}")
    print(f"  H298_cbs: range=[{h_min:.1f}, {h_max:.1f}]")

    # Apply physical filter to remove extreme outliers
    df_filtered = df[(df['H298_cbs'] >= -1000) & (df['H298_cbs'] <= 1000)]

    df_filtered_pct = len(df_filtered) / len(df) * 100
    print(
        f"\n[FIXTURE] After filter [-1000, +1000]: "
        f"{len(df_filtered)} molecules ({df_filtered_pct:.1f}%)"
    )
    h_filtered_mean = df_filtered["H298_cbs"].mean()
    h_filtered_std = df_filtered["H298_cbs"].std()
    print(f"  H298_cbs: mean={h_filtered_mean:.1f}, std={h_filtered_std:.1f}")

    # Sample 1000 with reproducibility
    df_sample = df_filtered.sample(n=1000, random_state=42)

    print(f"\n[FIXTURE] Final sample: {len(df_sample)} molecules")
    h_sample_mean = df_sample["H298_cbs"].mean()
    h_sample_std = df_sample["H298_cbs"].std()
    h_sample_min = df_sample["H298_cbs"].min()
    h_sample_max = df_sample["H298_cbs"].max()
    print(f"  H298_cbs: mean={h_sample_mean:.1f}, std={h_sample_std:.1f}")
    print(f"  H298_cbs: range=[{h_sample_min:.1f}, {h_sample_max:.1f}]")

    # Validate: std should be reasonable (not 7000+)
    assert h_sample_std < 500, (
        f"Filtered data should have std < 500, got {h_sample_std:.1f}"
    )

    # Extract basic features
    X_basic = (
        df_sample[["nheavy", "charge", "multiplicity"]].values.astype(float)
    )
    y_cbs = df_sample['H298_cbs'].values.astype(float)

    # Create enriched features
    X = create_enriched_features(X_basic)

    # Generate mock PM7 with realistic errors
    y_pm7 = create_realistic_mock_pm7(y_cbs, X_basic, seed=42)

    # Log delta statistics
    y_delta = y_cbs - y_pm7
    print("\n[FIXTURE] Delta (y_cbs - y_pm7):")
    print(f"  mean={y_delta.mean():.2f}, std={y_delta.std():.2f}")
    print(f"  range=[{y_delta.min():.2f}, {y_delta.max():.2f}]")

    return X, y_cbs, y_pm7


@pytest.fixture
def real_data_1k_extreme() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load 1000 molecules with EXTREME distribution (unfiltered).

    WARNING: This fixture creates pathological data for STRESS TESTING ONLY!

    Characteristics:
    - H298_cbs range: [-325407, +164949] kcal/mol
    - std: ~7230 kcal/mol
    - Random sample creates severe distribution shift
    - Train/test means can differ by 500+ kcal/mol

    Purpose:
    - Test delta-learning robustness to outliers
    - NOT for hypothesis validation
    - Demonstrates PM7 baseline stability

    Returns:
        Tuple of (X_enriched, y_cbs, y_pm7_mock)

    Example:
        >>> def test_stress(real_data_1k_extreme):
        ...     X, y_cbs, y_pm7 = real_data_1k_extreme
        ...     # Expect RMSE ~1000 for direct model!
    """
    from grimperium.data.loader import ChemperiumLoader

    # Load all data (no filter!)
    loader = ChemperiumLoader()
    df = loader.load_thermo_cbs_opt()

    print("\n[STRESS FIXTURE] Loading EXTREME distribution (unfiltered)")
    print(f"  Total molecules: {len(df)}")
    h_mean = df["H298_cbs"].mean()
    h_std = df["H298_cbs"].std()
    h_min = df["H298_cbs"].min()
    h_max = df["H298_cbs"].max()
    print(f"  H298_cbs: mean={h_mean:.1f}, std={h_std:.1f}")
    print(f"  H298_cbs: range=[{h_min:.1f}, {h_max:.1f}]")
    print("  WARNING: Distribution shift expected between train/test!")

    # Sample 1000 (includes outliers)
    df_sample = df.sample(n=1000, random_state=42)

    print(f"\n[STRESS FIXTURE] Sample: {len(df_sample)} molecules")
    h_sample_mean = df_sample["H298_cbs"].mean()
    h_sample_std = df_sample["H298_cbs"].std()
    print(f"  H298_cbs: mean={h_sample_mean:.1f}, std={h_sample_std:.1f}")

    # Extract basic features
    X_basic = (
        df_sample[["nheavy", "charge", "multiplicity"]].values.astype(float)
    )
    y_cbs = df_sample['H298_cbs'].values.astype(float)

    # Create enriched features
    X = create_enriched_features(X_basic)

    # Generate mock PM7 (note: magnitude bias is large for outliers)
    y_pm7 = create_realistic_mock_pm7(
        y_cbs,
        X_basic,
        seed=42,
        magnitude_bias_std=3,
    )

    return X, y_cbs, y_pm7


@pytest.fixture
def synthetic_data_1k() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Synthetic data with controlled distribution (1000 samples).

    Use for: Fast unit tests (<1s), edge case testing

    Returns:
        Tuple of (X_enriched, y_cbs, y_pm7_mock)
    """
    np.random.seed(42)
    n = 1000

    # Generate basic features
    nheavy = np.random.randint(1, 25, size=n)
    charge = np.random.choice([-2, -1, 0, 1, 2], size=n)
    mult = np.random.choice([1, 2, 3], size=n)

    X_basic = np.column_stack([nheavy, charge, mult]).astype(float)

    # Generate synthetic CBS energies (realistic range)
    y_cbs = np.random.normal(loc=-40, scale=50, size=n)

    # Create enriched features
    X = create_enriched_features(X_basic)

    # Generate mock PM7
    y_pm7 = create_realistic_mock_pm7(y_cbs, X_basic, seed=42)

    return X, y_cbs, y_pm7
