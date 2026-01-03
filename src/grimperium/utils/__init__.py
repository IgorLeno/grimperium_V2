"""
Utilities module for Grimperium.

This module provides utility functions for:
    - Logging: Configurable logging with multiple handlers
    - Validation: Input validation and data quality checks
    - Feature Engineering: Morgan FP + RDKit descriptor computation

"""

from grimperium.utils.feature_engineering import (
    FeatureEngineer,
    compute_morgan_fingerprints,
    compute_rdkit_descriptors,
)
from grimperium.utils.logging import get_logger, setup_logging
from grimperium.utils.validation import (
    validate_dataframe,
    validate_features,
    validate_smiles,
)

__all__ = [
    "setup_logging",
    "get_logger",
    "validate_smiles",
    "validate_dataframe",
    "validate_features",
    "FeatureEngineer",
    "compute_morgan_fingerprints",
    "compute_rdkit_descriptors",
]
