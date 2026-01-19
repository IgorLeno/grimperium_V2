"""
Core algorithms module for Grimperium.

This module contains:
    - Delta-learning logic (DeltaLearner, metrics)
    - Batch processing components (v2.2)
    - CSV data loading with validation
    - Value conversion utilities

Delta-learning approach:
    delta = H298_CBS - H298_PM7
    prediction = H298_PM7 + model.predict(features)

Batch processing components (v2.2):
    - BatchOrchestrator: End-to-end batch coordination
    - BatchScheduler: Molecule prioritization
    - BatchDataManager: CSV loading and validation
    - Molecule: Type-safe molecule representation
    - MoleculeValueConverter: CSV value conversion
"""

# Delta-learning components
from grimperium.core.delta_learning import DeltaLearner
from grimperium.core.metrics import (
    mae,
    mean_absolute_percentage_error,
    r2_score,
    rmse,
)

# Batch processing components (v2.2)
from grimperium.core.batch_orchestrator import (
    BatchOrchestrator,
    BatchOrchestratorError,
    BatchSummary,
    CalculationSettings,
)
from grimperium.core.batch_scheduler import BatchScheduler
from grimperium.core.csv_data_loader import (
    BatchDataManager,
    CSVDataLoader,
    CSVDataLoaderError,
    ValidationReport,
)
from grimperium.core.molecule import (
    Molecule,
    MoleculeIdentity,
    MoleculeMeta,
    MoleculeProperties,
    MoleculeResults,
    MoleculeStatus,
)
from grimperium.core.value_converter import (
    ConversionError,
    MoleculeValueConverter,
)

__all__ = [
    # Delta-learning
    "DeltaLearner",
    "rmse",
    "mae",
    "r2_score",
    "mean_absolute_percentage_error",
    # Batch processing (v2.2)
    "BatchOrchestrator",
    "BatchOrchestratorError",
    "BatchSummary",
    "CalculationSettings",
    "BatchScheduler",
    "BatchDataManager",
    "CSVDataLoader",
    "CSVDataLoaderError",
    "ValidationReport",
    "Molecule",
    "MoleculeIdentity",
    "MoleculeMeta",
    "MoleculeProperties",
    "MoleculeResults",
    "MoleculeStatus",
    "MoleculeValueConverter",
    "ConversionError",
]
