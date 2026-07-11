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
    WorkerStatus,
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
from grimperium.crest_pm7.batch.output_contracts import (
    BATCH_STATE_COLUMNS,
    CALCULATION_RESULTS_CSV_FILENAME,
    CALCULATION_RESULTS_XLSX_FILENAME,
    SCIENTIFIC_RESULT_COLUMNS,
    BatchOutputLayout,
    BatchResultWriteReport,
    write_batch_calculation_results,
)
from grimperium.crest_pm7.batch.processor_adapter import (
    FixedTimeoutPredictor,
    FixedTimeoutProcessor,
)
from grimperium.crest_pm7.batch.result_applier import (
    BatchResultApplier,
    ResultApplicationDecision,
)
from grimperium.crest_pm7.batch.result_ledger import (
    LedgerDecision,
    LedgerStatus,
    ResultLedger,
)
from grimperium.crest_pm7.batch.state_manager import BatchStateManager

__all__ = [
    # Enums
    "MoleculeStatus",
    "BatchSortingStrategy",
    "BatchFailurePolicy",
    "WorkerStatus",
    # Models
    "Batch",
    "BatchMolecule",
    "BatchResult",
    "BatchRowCSV",
    "ConformerDetail",
    "MoleculeDetail",
    "ArtifactPaths",
    "BatchOutputLayout",
    "BatchResultWriteReport",
    # Managers
    "BatchCSVManager",
    "BatchStateManager",
    "ConformerDetailManager",
    "BatchExecutionManager",
    "ArtifactManager",
    "BatchResultApplier",
    "ResultApplicationDecision",
    "ResultLedger",
    "LedgerDecision",
    "LedgerStatus",
    # Processor Adapter
    "FixedTimeoutPredictor",
    "FixedTimeoutProcessor",
    # Factory
    "create_execution_manager",
    # PR6 split-output contracts
    "BATCH_STATE_COLUMNS",
    "SCIENTIFIC_RESULT_COLUMNS",
    "CALCULATION_RESULTS_CSV_FILENAME",
    "CALCULATION_RESULTS_XLSX_FILENAME",
    "write_batch_calculation_results",
]
