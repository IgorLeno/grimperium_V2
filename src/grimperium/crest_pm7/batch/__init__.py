"""Batch processing module for CREST PM7 pipeline.

This module provides functionality for processing large batches of molecules
through the CREST PM7 pipeline with:
- CSV-based status tracking
- JSON detail files per molecule
- Fixed timeouts per batch
- Retry logic with configurable policies
"""

from grimperium.crest_pm7.batch.artifact_manager import ArtifactManager, ArtifactPaths
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.detail_manager import ConformerDetailManager
from grimperium.crest_pm7.batch.enums import (
    BatchFailurePolicy,
    BatchSortingStrategy,
    MoleculeStatus,
)
from grimperium.crest_pm7.batch.execution_manager import (
    BatchExecutionManager,
    create_execution_manager,
)
from grimperium.crest_pm7.batch.models import (
    Batch,
    BatchMolecule,
    BatchResult,
    BatchRowCSV,
    ConformerDetail,
    MoleculeDetail,
)
from grimperium.crest_pm7.batch.processor_adapter import (
    FixedTimeoutPredictor,
    FixedTimeoutProcessor,
)

__all__ = [
    # Enums
    "MoleculeStatus",
    "BatchSortingStrategy",
    "BatchFailurePolicy",
    # Models
    "Batch",
    "BatchMolecule",
    "BatchResult",
    "BatchRowCSV",
    "ConformerDetail",
    "MoleculeDetail",
    "ArtifactPaths",
    # Managers
    "BatchCSVManager",
    "ConformerDetailManager",
    "BatchExecutionManager",
    "ArtifactManager",
    # Processor Adapter
    "FixedTimeoutPredictor",
    "FixedTimeoutProcessor",
    # Factory
    "create_execution_manager",
]
