"""
Global configuration for Grimperium.

This module provides configuration management via dataclasses,
supporting both programmatic and file-based configuration.

Example:
    >>> from grimperium.config import GrimperiumConfig
    >>> config = GrimperiumConfig()
    >>> config.model.krr_alpha = 0.1
    >>> config.features.morgan_bits = 512

"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class ModelConfig:
    """Configuration for ML models."""

    # Kernel Ridge Regression
    krr_alpha: float = 1.0
    krr_kernel: str = "rbf"
    krr_gamma: Optional[float] = None

    # XGBoost
    xgb_n_estimators: int = 100
    xgb_max_depth: int = 6
    xgb_learning_rate: float = 0.1
    xgb_subsample: float = 0.8

    # Ensemble
    ensemble_weights: tuple[float, float] = (0.5, 0.5)


@dataclass
class FeatureConfig:
    """Configuration for feature engineering."""

    # Morgan Fingerprints
    morgan_bits: int = 256
    morgan_radius: int = 2

    # RDKit descriptors
    use_rdkit_descriptors: bool = True
    rdkit_descriptors: list[str] = field(
        default_factory=lambda: ["MolWt", "TPSA", "MolLogP", "NumRotatableBonds"]
    )

    # Tabular features
    tabular_features: list[str] = field(
        default_factory=lambda: ["nheavy", "charge", "multiplicity"]
    )


@dataclass
class DataConfig:
    """Configuration for data loading and processing."""

    # Paths
    chemperium_path: Optional[Path] = None
    pm7_cache_dir: Optional[Path] = None

    # Processing
    test_size: float = 0.2
    random_state: int = 42
    n_cv_folds: int = 5


@dataclass
class GrimperiumConfig:
    """
    Main configuration class for Grimperium.

    Aggregates all sub-configurations for models, features, and data.

    Attributes:
        model: Model hyperparameters configuration
        features: Feature engineering configuration
        data: Data loading and processing configuration
        verbose: Enable verbose logging
        n_jobs: Number of parallel jobs (-1 for all cores)

    Example:
        >>> config = GrimperiumConfig()
        >>> config.model.krr_alpha = 0.5
        >>> config.features.morgan_bits = 512
        >>> config.data.test_size = 0.15

    """

    model: ModelConfig = field(default_factory=ModelConfig)
    features: FeatureConfig = field(default_factory=FeatureConfig)
    data: DataConfig = field(default_factory=DataConfig)

    verbose: bool = True
    n_jobs: int = -1

    def to_dict(self) -> dict:
        """Convert configuration to dictionary."""
        raise NotImplementedError("Will be implemented in future batch")

    @classmethod
    def from_dict(cls, config_dict: dict) -> "GrimperiumConfig":
        """Create configuration from dictionary."""
        raise NotImplementedError("Will be implemented in future batch")

    @classmethod
    def from_yaml(cls, path: Path) -> "GrimperiumConfig":
        """Load configuration from YAML file."""
        raise NotImplementedError("Will be implemented in future batch")
