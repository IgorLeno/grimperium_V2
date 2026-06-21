"""
ML models module for Grimperium.

This module provides the ensemble models for delta-learning:
    - BaseModel: Abstract base class for all models
    - KernelRidgeRegressor: KRR with RBF kernel for smooth corrections
    - XGBoostRegressor: Gradient boosting for complex patterns
    - DeltaLearningEnsemble: Combined ensemble with weighted averaging

Architecture:
    BaseModel (ABC)
    ├── KernelRidgeRegressor
    ├── XGBoostRegressor
    └── DeltaLearningEnsemble (combines KRR + XGB)

"""

from grimperium.models.base import BaseModel
from grimperium.models.kernel_ridge import KernelRidgeRegressor
from grimperium.models.xgboost_model import XGBoostRegressor

__all__ = [
    "BaseModel",
    "KernelRidgeRegressor",
    "XGBoostRegressor",
    "DeltaLearningEnsemble",
]


def __getattr__(name: str) -> object:
    """Lazily resolve DeltaLearningEnsemble to avoid a circular import.

    grimperium.models.delta_ensemble re-exports DeltaLearningEnsemble from
    grimperium.ml.delta_ensemble, which itself imports KernelRidgeRegressor
    and XGBoostRegressor from this package. Importing the shim eagerly here
    creates an import cycle when grimperium.ml.delta_ensemble is imported
    first.
    """
    if name == "DeltaLearningEnsemble":
        from grimperium.models.delta_ensemble import DeltaLearningEnsemble

        return DeltaLearningEnsemble
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
