"""ML pipeline for Grimperium delta-learning.

Modules:
    data_loader — Load and filter CSV for ML training
    features    — Feature engineering with imputation
    trainer     — Train/test split and training orchestration
    evaluator   — Full-dataset evaluation
"""

from grimperium.ml.data_loader import load_ml_data
from grimperium.ml.evaluator import evaluate
from grimperium.ml.features import FeaturePipeline
from grimperium.ml.trainer import train

__all__ = ["FeaturePipeline", "evaluate", "load_ml_data", "train"]
