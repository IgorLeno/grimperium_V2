"""
Core algorithms module for Grimperium.

This module contains the core delta-learning logic and metrics:
    - DeltaLearner: Orchestrates delta computation and model training
    - Metrics: RMSE, MAE, R² and other evaluation metrics

The delta-learning approach:
    delta = H298_CBS - H298_PM7
    prediction = H298_PM7 + model.predict(features)

"""

from grimperium.core.delta_learning import DeltaLearner
from grimperium.core.metrics import (
    mae,
    mean_absolute_percentage_error,
    r2_score,
    rmse,
)

__all__ = [
    "DeltaLearner",
    "rmse",
    "mae",
    "r2_score",
    "mean_absolute_percentage_error",
]
