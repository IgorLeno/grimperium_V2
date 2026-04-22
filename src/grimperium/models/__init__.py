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
from grimperium.models.delta_ensemble import DeltaLearningEnsemble

__all__ = [
    "BaseModel",
    "KernelRidgeRegressor",
    "XGBoostRegressor",
    "DeltaLearningEnsemble",
]
