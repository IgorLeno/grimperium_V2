"""
Pytest configuration and shared fixtures for Grimperium tests.

This module provides:
    - Mock data fixtures for Chemperium dataset
    - Mock PM7 calculation results
    - Common test utilities

All fixtures use in-memory data to enable fast testing
without external dependencies.

"""

import numpy as np
import pandas as pd
import pytest

# =============================================================================
# Mock Chemperium Data
# =============================================================================


@pytest.fixture
def sample_smiles() -> list[str]:
    """Sample SMILES strings for testing."""
    return [
        "C",  # Methane
        "CC",  # Ethane
        "CCC",  # Propane
        "CCO",  # Ethanol
        "CC(=O)O",  # Acetic acid
        "c1ccccc1",  # Benzene
        "CC(C)C",  # Isobutane
        "CCCO",  # 1-Propanol
        "CC=O",  # Acetaldehyde
        "COC",  # Dimethyl ether
    ]


@pytest.fixture
def mock_chemperium_df(sample_smiles: list[str]) -> pd.DataFrame:
    """
    Mock Chemperium DataFrame with realistic thermodynamic values.

    Contains:
        - smiles: SMILES strings
        - charge: Molecular charge (mostly 0)
        - multiplicity: Spin multiplicity (mostly 1)
        - nheavy: Number of heavy atoms
        - H298_cbs: CBS enthalpy (reference, mock values)
        - H298_b3: B3LYP enthalpy (mock values)
        - S298: Entropy at 298K (mock values)
        - cp_1 to cp_5: Heat capacity columns (mock, subset)

    """
    n_molecules = len(sample_smiles)
    np.random.seed(42)

    # Realistic mock thermodynamic values
    # H298 in kcal/mol, typical range for small organics
    h298_cbs = np.array(
        [
            -17.9,  # Methane
            -20.0,  # Ethane
            -25.0,  # Propane
            -56.2,  # Ethanol
            -103.4,  # Acetic acid
            19.8,  # Benzene
            -32.1,  # Isobutane
            -61.0,  # 1-Propanol
            -39.7,  # Acetaldehyde
            -44.0,  # Dimethyl ether
        ]
    )

    # B3LYP values (slightly different from CBS)
    h298_b3 = h298_cbs + np.random.normal(0, 2, n_molecules)

    # Entropy values (cal/mol·K)
    s298 = np.array([44.5, 54.9, 64.5, 67.5, 68.0, 64.3, 70.4, 77.0, 63.0, 63.7])

    # Heat capacity values (cal/mol·K) - just 5 columns for mock
    cp_base = np.array([8.5, 12.6, 17.6, 15.7, 16.0, 19.5, 23.2, 21.5, 13.2, 15.3])

    df = pd.DataFrame(
        {
            "smiles": sample_smiles,
            "charge": [0] * n_molecules,
            "multiplicity": [1] * n_molecules,
            "nheavy": [1, 2, 3, 3, 4, 6, 4, 4, 3, 3],
            "H298_cbs": h298_cbs,
            "H298_b3": h298_b3,
            "S298": s298,
        }
    )

    # Add heat capacity columns
    for i in range(1, 6):
        df[f"cp_{i}"] = cp_base + i * 0.5 + np.random.normal(0, 0.1, n_molecules)

    return df


@pytest.fixture
def mock_pm7_df(sample_smiles: list[str]) -> pd.DataFrame:
    """
    Mock PM7 calculation results.

    PM7 values have systematic error compared to CBS,
    which is what delta-learning aims to correct.

    """
    n_molecules = len(sample_smiles)

    # PM7 values with systematic underestimation
    h298_pm7 = np.array(
        [
            -15.2,  # Methane (error: ~2.7 kcal/mol)
            -17.8,  # Ethane (error: ~2.2 kcal/mol)
            -22.5,  # Propane (error: ~2.5 kcal/mol)
            -52.0,  # Ethanol (error: ~4.2 kcal/mol)
            -98.0,  # Acetic acid (error: ~5.4 kcal/mol)
            22.5,  # Benzene (error: ~2.7 kcal/mol)
            -29.0,  # Isobutane (error: ~3.1 kcal/mol)
            -57.0,  # 1-Propanol (error: ~4.0 kcal/mol)
            -36.5,  # Acetaldehyde (error: ~3.2 kcal/mol)
            -41.0,  # Dimethyl ether (error: ~3.0 kcal/mol)
        ]
    )

    # Garantia de consistência: o fixture `sample_smiles` e os valores mock devem ter o mesmo tamanho.
    assert h298_pm7.shape[0] == n_molecules

    return pd.DataFrame(
        {
            "smiles": sample_smiles,
            "H298_pm7": h298_pm7,
        }
    )


@pytest.fixture
def mock_merged_df(
    mock_chemperium_df: pd.DataFrame,
    mock_pm7_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merged DataFrame with delta column."""
    merged = mock_chemperium_df.merge(mock_pm7_df, on="smiles")
    merged["delta_pm7"] = merged["H298_cbs"] - merged["H298_pm7"]
    return merged


# =============================================================================
# Mock Features
# =============================================================================


@pytest.fixture
def mock_features(sample_smiles: list[str]) -> np.ndarray:
    """
    Mock feature matrix for testing.

    Shape: (n_samples, n_features) where n_features includes:
        - 3 tabular features
        - 256 Morgan FP bits (mock: reduced to 32 for speed)
        - 10 RDKit descriptors

    Total: 45 features (reduced for testing)

    """
    n_samples = len(sample_smiles)
    np.random.seed(42)

    # Tabular features: nheavy, charge, multiplicity
    tabular = np.array(
        [
            [1, 0, 1],
            [2, 0, 1],
            [3, 0, 1],
            [3, 0, 1],
            [4, 0, 1],
            [6, 0, 1],
            [4, 0, 1],
            [4, 0, 1],
            [3, 0, 1],
            [3, 0, 1],
        ]
    )

    # Mock Morgan fingerprints (32 bits for testing, normally 256)
    morgan_fp = np.random.randint(0, 2, (n_samples, 32))

    # Mock RDKit descriptors (10 features)
    rdkit_desc = np.random.randn(n_samples, 10)

    return np.hstack([tabular, morgan_fp, rdkit_desc])


@pytest.fixture
def mock_targets(mock_merged_df: pd.DataFrame) -> np.ndarray:
    """Mock delta targets for training."""
    return mock_merged_df["delta_pm7"].values


# =============================================================================
# Train/Test Splits
# =============================================================================


@pytest.fixture
def train_test_split(
    mock_features: np.ndarray,
    mock_targets: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Split mock data into train/test sets (80/20)."""
    n_train = 8
    X_train = mock_features[:n_train]
    X_test = mock_features[n_train:]
    y_train = mock_targets[:n_train]
    y_test = mock_targets[n_train:]
    return X_train, X_test, y_train, y_test


# =============================================================================
# Configuration Fixtures
# =============================================================================


@pytest.fixture
def default_config():
    """Default Grimperium configuration."""
    from grimperium.config import GrimperiumConfig

    return GrimperiumConfig()


# =============================================================================
# Utility Functions
# =============================================================================


def assert_array_almost_equal(
    actual: np.ndarray,
    expected: np.ndarray,
    decimal: int = 5,
) -> None:
    """Assert arrays are almost equal."""
    np.testing.assert_array_almost_equal(actual, expected, decimal=decimal)


def assert_valid_predictions(
    predictions: np.ndarray,
    n_samples: int,
) -> None:
    """Assert predictions are valid."""
    assert predictions.shape == (n_samples,)
    assert not np.any(np.isnan(predictions))
    assert not np.any(np.isinf(predictions))
