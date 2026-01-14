"""Conformer detail manager for per-molecule JSON files.

This module provides ConformerDetailManager for:
- Creating and saving detailed JSON files per molecule
- Converting PM7Result to MoleculeDetail
- Loading existing detail files

Thread-safety:
- SAFE for concurrent access to DIFFERENT mol_ids
- NOT safe for concurrent access to the SAME mol_id
- Per-mol locking not implemented; serialize access per mol_id in caller
"""

import json
import logging
import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional

import pydantic

from grimperium.crest_pm7.batch.models import ConformerDetail, MoleculeDetail

if TYPE_CHECKING:
    from grimperium.crest_pm7.molecule_processor import PM7Result

LOG = logging.getLogger("grimperium.crest_pm7.batch.detail_manager")


class ConformerDetailManager:
    """Manages JSON detail files for processed molecules.

    Each molecule gets a separate JSON file containing:
    - Full CREST results
    - Per-conformer MOPAC results
    - Energy values and quality grades

    Thread-safety:
    - SAFE for concurrent access to DIFFERENT mol_ids
    - NOT safe for concurrent access to the SAME mol_id
    - Per-mol locking not implemented; serialize access per mol_id in caller

    Attributes:
        detail_dir: Directory for JSON detail files
    """

    def __init__(self, detail_dir: Path) -> None:
        """Initialize detail manager.

        Args:
            detail_dir: Directory for JSON files (created if needed)
        """
        self.detail_dir = Path(detail_dir)
        self.detail_dir.mkdir(parents=True, exist_ok=True)
        LOG.debug(f"ConformerDetailManager initialized at {self.detail_dir}")

    def get_detail_path(self, mol_id: str) -> Path:
        """Get path for molecule detail file.

        Args:
            mol_id: Molecule identifier

        Returns:
            Path to JSON file
        """
        # Sanitize mol_id for filesystem (replace special chars)
        safe_id = mol_id.replace("/", "_").replace("\\", "_")
        return self.detail_dir / f"{safe_id}.json"

    def save_detail(self, detail: MoleculeDetail) -> Path:
        """Save molecule detail to JSON file atomically.

        Uses atomic write pattern: write to temp file, fsync, then rename.
        This prevents corruption if interrupted during write.
        Ensures file descriptor is always properly closed, even on exception.

        Process:
        1. Create temp file with mkstemp
        2. Try to open as file object (may fail)
        3. If fdopen fails: close raw FD and re-raise
        4. If fdopen succeeds: use with-block to handle closing
        5. After successful write: atomic rename
        6. On any error: cleanup temp file

        Args:
            detail: MoleculeDetail to save

        Returns:
            Path to saved file

        Raises:
            Exception: If write fails (temp file is cleaned up)
        """
        detail_path = self.get_detail_path(detail.mol_id)

        # Write to temp file in same directory (ensure same filesystem)
        temp_fd, temp_path = tempfile.mkstemp(
            dir=detail_path.parent, text=True, prefix=".tmp_"
        )

        try:
            # Step 1: Try to open FD as file object
            # This is in its own try/except so fdopen failure doesn't trigger
            # the double-close issue
            f = None
            try:
                f = os.fdopen(temp_fd, "w", encoding="utf-8")
            except Exception as e:
                # fdopen failed — close raw FD before re-raising
                os.close(temp_fd)
                LOG.error(f"Failed to open temp file for {detail.mol_id}: {e}")
                raise

            # Step 2: fdopen succeeded — use with-block to handle closing
            # Do NOT call os.close(temp_fd) anywhere in this path
            try:
                with f:
                    json.dump(detail.model_dump(mode="json"), f, indent=2)
                    f.flush()
                    os.fsync(f.fileno())  # Force to disk
            except Exception as e:
                LOG.error(f"Failed to write detail for {detail.mol_id}: {e}")
                raise

            # Step 3: Atomic rename
            os.replace(temp_path, detail_path)
            LOG.debug(f"Saved detail for {detail.mol_id}")
            return detail_path

        except Exception as e:
            # Cleanup: remove temp file if it exists
            try:
                os.unlink(temp_path)
            except FileNotFoundError:
                pass
            LOG.error(f"Failed to save detail for {detail.mol_id}: {e}")
            raise

    def load_detail(self, mol_id: str) -> Optional[MoleculeDetail]:
        """Load molecule detail from JSON file.

        Args:
            mol_id: Molecule identifier

        Returns:
            MoleculeDetail if file exists and valid, None otherwise
        """
        path = self.get_detail_path(mol_id)

        if not path.exists():
            return None

        try:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
            return MoleculeDetail.model_validate(data)

        except (json.JSONDecodeError, pydantic.ValidationError) as e:
            LOG.error(f"Failed to load detail for {mol_id} at {path}: {e}")
            return None

    def exists(self, mol_id: str) -> bool:
        """Check if detail file exists for molecule.

        Args:
            mol_id: Molecule identifier

        Returns:
            True if file exists
        """
        return self.get_detail_path(mol_id).exists()

    def pm7result_to_detail(
        self,
        mol_id: str,
        smiles: str,
        result: "PM7Result",
        batch_id: str,
    ) -> MoleculeDetail:
        """Convert PM7Result to MoleculeDetail.

        Args:
            mol_id: Molecule identifier
            smiles: SMILES string
            result: PM7Result from processing
            batch_id: Batch that processed this molecule

        Returns:
            MoleculeDetail ready to save
        """
        # Convert conformers to ConformerDetail list
        conformer_details: list[ConformerDetail] = []
        for i, conf in enumerate(result.conformers):
            detail = ConformerDetail(
                conformer_index=i,
                energy_hof=conf.energy_hof,
                energy_total=conf.energy_total,
                mopac_status=(
                    conf.mopac_status.value if conf.mopac_status else "NOT_ATTEMPTED"
                ),
                mopac_time=conf.mopac_time,
                mopac_error=conf.mopac_error,
                geometry_file=str(conf.output_path) if conf.output_path else None,
            )
            conformer_details.append(detail)

        # Count successful conformers
        num_successful = len(result.successful_conformers)

        return MoleculeDetail(
            mol_id=mol_id,
            smiles=smiles,
            batch_id=batch_id,
            timestamp=result.timestamp,
            crest_status=result.crest_status.value,
            crest_conformers_generated=result.crest_conformers_generated,
            crest_time=result.crest_time,
            crest_error=result.crest_error,
            conformers=conformer_details,
            num_conformers_selected=result.num_conformers_selected or 0,
            num_conformers_successful=num_successful,
            most_stable_hof=result.most_stable_hof,
            quality_grade=result.quality_grade.value,
            issues=result.issues,
            success=result.success,
            error_message=result.error_message,
            total_execution_time=result.total_execution_time,
        )

    def list_details(self) -> list[str]:
        """List all mol_ids with detail files.

        Returns:
            List of mol_ids
        """
        return [p.stem for p in self.detail_dir.glob("*.json")]

    def delete_detail(self, mol_id: str) -> bool:
        """Delete detail file for molecule.

        Args:
            mol_id: Molecule identifier

        Returns:
            True if file was deleted, False if didn't exist
        """
        path = self.get_detail_path(mol_id)

        try:
            path.unlink()  # Raises FileNotFoundError if not exists
            LOG.debug(f"Deleted detail for {mol_id}")
            return True
        except FileNotFoundError:
            LOG.debug(f"No detail file found for {mol_id}, nothing to delete")
            return False

    def get_summary_stats(self) -> dict[str, Any]:
        """Get summary statistics from all detail files.

        Returns:
            Dict with counts and aggregates
        """
        mol_ids = self.list_details()

        if not mol_ids:
            return {
                "total_files": 0,
                "success_count": 0,
                "failed_count": 0,
            }

        success_count = 0
        failed_count = 0
        total_conformers = 0
        hof_values: list[float] = []

        for mol_id in mol_ids:
            detail = self.load_detail(mol_id)
            if detail is None:
                continue

            if detail.success:
                success_count += 1
                if detail.most_stable_hof is not None:
                    hof_values.append(detail.most_stable_hof)
            else:
                failed_count += 1

            total_conformers += detail.num_conformers_successful

        return {
            "total_files": len(mol_ids),
            "success_count": success_count,
            "failed_count": failed_count,
            "total_conformers": total_conformers,
            "min_hof": min(hof_values) if hof_values else None,
            "max_hof": max(hof_values) if hof_values else None,
        }
