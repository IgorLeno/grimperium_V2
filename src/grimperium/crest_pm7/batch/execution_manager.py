"""Batch execution manager for CREST PM7 pipeline.

This module provides BatchExecutionManager for:
- Orchestrating batch execution
- Coordinating CSV tracking, detail files, and processing
- Managing failure policies (PARTIAL_OK, ALL_OR_NOTHING)
- Generating BatchResult with statistics
"""

import logging
import time
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path

from grimperium.crest_pm7.batch.artifact_manager import ArtifactManager
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.detail_manager import ConformerDetailManager
from grimperium.crest_pm7.batch.enums import BatchFailurePolicy, MoleculeStatus
from grimperium.crest_pm7.batch.models import Batch, BatchResult
from grimperium.crest_pm7.batch.processor_adapter import FixedTimeoutProcessor
from grimperium.crest_pm7.config import PM7Config

LOG = logging.getLogger("grimperium.crest_pm7.batch.execution_manager")


class BatchExecutionManager:
    """Orchestrates batch execution of CREST PM7 pipeline.

    Coordinates:
    - BatchCSVManager: Status tracking in CSV
    - ConformerDetailManager: Per-molecule JSON files
    - FixedTimeoutProcessor: Molecule processing with fixed timeouts

    Implements two failure policies:
    - PARTIAL_OK: Failed molecules marked RERUN/SKIP individually
    - ALL_OR_NOTHING: If any fail, reset entire batch to PENDING

    Attributes:
        csv_manager: Manager for CSV status tracking
        detail_manager: Manager for JSON detail files
        pm7_config: PM7 configuration
        processor_adapter: Fixed timeout processor
    """

    def __init__(
        self,
        csv_manager: BatchCSVManager,
        detail_manager: ConformerDetailManager,
        pm7_config: PM7Config,
        processor_adapter: FixedTimeoutProcessor | None = None,
        artifact_manager: ArtifactManager | None = None,
    ) -> None:
        """Initialize batch execution manager.

        Args:
            csv_manager: Manager for CSV tracking
            detail_manager: Manager for JSON detail files
            pm7_config: PM7 configuration
            processor_adapter: Optional processor adapter (created if None)
            artifact_manager: Optional artifact manager for debug/audit files
        """
        self.csv_manager = csv_manager
        self.detail_manager = detail_manager
        self.pm7_config = pm7_config
        self.processor_adapter = processor_adapter or FixedTimeoutProcessor(pm7_config)
        self.artifact_manager = artifact_manager

        LOG.info("BatchExecutionManager initialized")

    def execute_batch(
        self,
        batch: Batch,
        progress_callback: Callable[[str, int, int], None] | None = None,
    ) -> BatchResult:
        """Execute a batch of molecules.

        Processes each molecule sequentially, updates CSV status,
        saves detail files, and handles failures according to policy.

        Args:
            batch: Batch to execute
            progress_callback: Optional callback(mol_id, current, total) for progress

        Returns:
            BatchResult with execution statistics
        """
        if batch.is_empty:
            LOG.warning(f"Batch {batch.batch_id} is empty, nothing to execute")
            return BatchResult(
                batch_id=batch.batch_id,
                total_count=0,
            )

        LOG.info(
            f"Starting batch {batch.batch_id}: {batch.size} molecules, "
            f"policy={batch.failure_policy.value}"
        )

        # Update processor timeouts
        self.processor_adapter.update_timeouts(
            crest_timeout_minutes=batch.crest_timeout_minutes,
            mopac_timeout_minutes=batch.mopac_timeout_minutes,
        )

        # Initialize result
        result = BatchResult(
            batch_id=batch.batch_id,
            total_count=batch.size,
            timestamp_start=datetime.now(timezone.utc),
        )

        # Track HOF values for statistics
        hof_values: list[tuple[str, float]] = []  # (mol_id, hof)
        start_time = time.time()

        # Process each molecule
        for i, mol in enumerate(batch.molecules, start=1):
            try:
                self._process_molecule(
                    mol_id=mol.mol_id,
                    smiles=mol.smiles,
                    batch_id=batch.batch_id,
                    batch_order=mol.batch_order,
                    crest_timeout=batch.crest_timeout_minutes,
                    mopac_timeout=batch.mopac_timeout_minutes,
                    result=result,
                    hof_values=hof_values,
                )
            except Exception as e:
                LOG.error(
                    f"Unexpected error processing {mol.mol_id}: {e}",
                    exc_info=True,
                )
                result.failed_count += 1
                result.failed_mol_ids.append(mol.mol_id)

            # Progress callback
            if progress_callback:
                try:
                    progress_callback(mol.mol_id, i, batch.size)
                except Exception as e:
                    LOG.warning(f"Progress callback error: {e}")

        # Finalize timing
        result.total_time = time.time() - start_time
        result.timestamp_end = datetime.now(timezone.utc)

        # Calculate HOF statistics
        if hof_values:
            # Find min and max simultaneously to avoid floating-point equality issues
            min_mol_id, min_hof_val = min(hof_values, key=lambda x: x[1])
            result.min_hof = round(min_hof_val, 2)
            result.min_hof_mol_id = min_mol_id

            max_mol_id, max_hof_val = max(hof_values, key=lambda x: x[1])
            result.max_hof = round(max_hof_val, 2)
            result.max_hof_mol_id = max_mol_id

        # Handle ALL_OR_NOTHING policy
        if batch.failure_policy == BatchFailurePolicy.ALL_OR_NOTHING and (
            result.failed_count > 0 or result.rerun_count > 0
        ):
            LOG.warning(
                f"ALL_OR_NOTHING: Resetting batch {batch.batch_id} due to failures"
            )
            reset_count = self.csv_manager.reset_batch(batch.batch_id)
            LOG.info(f"Reset {reset_count} molecules from batch {batch.batch_id}")

        LOG.info(
            f"Batch {batch.batch_id} complete: "
            f"{result.success_count}/{result.total_count} OK, "
            f"{result.rerun_count} rerun, {result.skip_count} skip, "
            f"{result.failed_count} failed, "
            f"time={result.total_time:.1f}s"
        )

        return result

    def _process_molecule(
        self,
        mol_id: str,
        smiles: str,
        batch_id: str,
        batch_order: int,
        crest_timeout: float,
        mopac_timeout: float,
        result: BatchResult,
        hof_values: list[tuple[str, float]],
    ) -> None:
        """Process a single molecule within batch context.

        Updates CSV status, saves detail file, and updates result statistics.

        Args:
            mol_id: Molecule identifier
            smiles: SMILES string
            batch_id: Current batch ID
            batch_order: Position in batch
            crest_timeout: CREST timeout (minutes)
            mopac_timeout: MOPAC timeout (minutes)
            result: BatchResult to update
            hof_values: List to append HOF values
        """
        LOG.info(f"Processing {mol_id} ({batch_order}/{result.total_count})")

        # Mark as running
        self.csv_manager.mark_running(mol_id)

        try:
            # Process molecule
            pm7_result = self.processor_adapter.process_with_fixed_timeout(
                mol_id=mol_id,
                smiles=smiles,
            )

            # Create CSV update dict
            csv_update = self.csv_manager.pm7result_to_csv_update(
                mol_id=mol_id,
                result=pm7_result,
                batch_id=batch_id,
                batch_order=batch_order,
                crest_timeout_used=crest_timeout,
                mopac_timeout_used=mopac_timeout,
            )

            # Save detail file
            detail = self.detail_manager.pm7result_to_detail(
                mol_id=mol_id,
                smiles=smiles,
                result=pm7_result,
                batch_id=batch_id,
            )
            self.detail_manager.save_detail(detail)

            # Save artifacts for debug/audit if artifact manager is configured
            if self.artifact_manager is not None:
                try:
                    # Work directories are typically temp_dir/{mol_id}/crest and /mopac
                    mol_work_dir = self.pm7_config.temp_dir / mol_id
                    crest_work_dir = mol_work_dir / "crest" if mol_work_dir.exists() else None
                    mopac_work_dir = mol_work_dir / "mopac" if mol_work_dir.exists() else None

                    self.artifact_manager.save_artifacts(
                        mol_id=mol_id,
                        batch_id=batch_id,
                        crest_work_dir=crest_work_dir,
                        mopac_work_dir=mopac_work_dir,
                        success=pm7_result.success,
                        extra_metadata={
                            "quality_grade": pm7_result.quality_grade.value,
                            "most_stable_hof": pm7_result.most_stable_hof,
                            "error_message": pm7_result.error_message,
                        },
                    )
                except Exception as e:
                    LOG.warning(f"Failed to save artifacts for {mol_id}: {e}")

            # Update status based on success
            if pm7_result.success:
                self.csv_manager.mark_success(mol_id, csv_update)
                result.success_count += 1

                # Track HOF for statistics
                if pm7_result.most_stable_hof is not None:
                    hof_values.append((mol_id, pm7_result.most_stable_hof))

                LOG.info(
                    f"{mol_id}: OK (HOF={pm7_result.most_stable_hof}, "
                    f"grade={pm7_result.quality_grade.value})"
                )
            else:
                # Mark for rerun (or skip if max retries)
                error_msg = pm7_result.error_message or "Unknown error"
                self.csv_manager.mark_rerun(mol_id, error_msg, csv_update)

                # Check if it became SKIP using public method
                new_status = self.csv_manager.get_status(mol_id)

                if new_status == MoleculeStatus.SKIP.value:
                    result.skip_count += 1
                else:
                    result.rerun_count += 1
                    result.rerun_mol_ids.append(mol_id)

                LOG.warning(f"{mol_id}: {new_status} - {error_msg}")

        except Exception as e:
            # Unexpected error - mark for rerun
            error_msg = f"Processing exception: {str(e)}"
            LOG.error(f"{mol_id}: {error_msg}", exc_info=True)

            self.csv_manager.mark_rerun(mol_id, error_msg)

            # Check final status to update correct counter
            new_status = self.csv_manager.get_status(mol_id)
            if new_status == MoleculeStatus.SKIP.value:
                result.skip_count += 1
            else:
                result.rerun_count += 1
                result.rerun_mol_ids.append(mol_id)

    def get_status_summary(self) -> dict[str, int]:
        """Get current status counts from CSV.

        Returns:
            Dict mapping status to count
        """
        return self.csv_manager.get_status_counts()

    def get_pending_count(self) -> int:
        """Get count of molecules still pending.

        Returns:
            Number of PENDING + RERUN molecules
        """
        counts = self.get_status_summary()
        return counts.get(MoleculeStatus.PENDING.value, 0) + counts.get(
            MoleculeStatus.RERUN.value, 0
        )

    def is_complete(self) -> bool:
        """Check if all molecules have been processed.

        Returns:
            True if no PENDING or RERUN molecules remain
        """
        return self.get_pending_count() == 0


def create_execution_manager(
    csv_path: Path,
    detail_dir: Path,
    pm7_config: PM7Config,
    crest_timeout_minutes: float = 30.0,
    mopac_timeout_minutes: float = 60.0,
    artifact_dir: Path | None = None,
    preserve_artifacts_on_success: bool = True,
    preserve_artifacts_on_failure: bool = True,
) -> BatchExecutionManager:
    """Factory function to create BatchExecutionManager with defaults.

    Args:
        csv_path: Path to CSV tracking file
        detail_dir: Directory for JSON detail files
        pm7_config: PM7 configuration
        crest_timeout_minutes: Default CREST timeout
        mopac_timeout_minutes: Default MOPAC timeout
        artifact_dir: Directory for debug/audit artifacts (None to disable)
        preserve_artifacts_on_success: Save artifacts for successful molecules
        preserve_artifacts_on_failure: Save artifacts for failed molecules

    Returns:
        Configured BatchExecutionManager
    """
    csv_manager = BatchCSVManager(csv_path)
    detail_manager = ConformerDetailManager(detail_dir)
    processor = FixedTimeoutProcessor(
        config=pm7_config,
        crest_timeout_minutes=crest_timeout_minutes,
        mopac_timeout_minutes=mopac_timeout_minutes,
    )

    # Create artifact manager if directory specified
    artifact_manager = None
    if artifact_dir is not None:
        artifact_manager = ArtifactManager(
            artifact_dir=artifact_dir,
            preserve_on_success=preserve_artifacts_on_success,
            preserve_on_failure=preserve_artifacts_on_failure,
        )

    return BatchExecutionManager(
        csv_manager=csv_manager,
        detail_manager=detail_manager,
        pm7_config=pm7_config,
        processor_adapter=processor,
        artifact_manager=artifact_manager,
    )
