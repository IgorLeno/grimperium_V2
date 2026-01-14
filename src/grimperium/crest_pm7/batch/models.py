"""Pydantic models for batch processing of CREST PM7 pipeline.

This module defines data models for:
- BatchResult: Summary statistics from batch execution
- Batch: Selected molecules for processing
- BatchRowCSV: Schema for CSV tracking file
- ConformerDetail: Per-molecule JSON detail file
"""

from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field, computed_field, field_serializer

from grimperium.crest_pm7.batch.enums import (
    BatchFailurePolicy,
    BatchSortingStrategy,
    MoleculeStatus,
)


class BatchResult(BaseModel):
    """Result of batch execution.

    Contains summary statistics and timing for display_batch_summary().
    Statistics are computed from THIS BATCH only, not the entire dataset.
    """

    # Identification
    batch_id: str = Field(..., description="Batch identifier")

    # Counts (from batch ONLY, not entire dataset)
    total_count: int = Field(..., description="Total molecules in this batch")
    success_count: int = Field(default=0, description="Molecules with status OK")
    rerun_count: int = Field(default=0, description="Molecules marked for Rerun")
    skip_count: int = Field(
        default=0, description="Molecules marked Skip (max retries)"
    )
    failed_count: int = Field(
        default=0, description="Molecules that failed processing"
    )

    # Timing (in seconds)
    total_time: float = Field(default=0.0, description="Total execution time")
    timestamp_start: datetime = Field(
        default_factory=lambda: datetime.now(timezone.utc)
    )
    timestamp_end: Optional[datetime] = Field(default=None)

    # Energy statistics (from successful molecules in THIS BATCH only)
    min_hof: Optional[float] = Field(default=None, description="Minimum HOF in batch")
    max_hof: Optional[float] = Field(default=None, description="Maximum HOF in batch")
    min_hof_mol_id: Optional[str] = Field(
        default=None, description="mol_id with min HOF"
    )
    max_hof_mol_id: Optional[str] = Field(
        default=None, description="mol_id with max HOF"
    )

    # Issues tracking (for debugging)
    failed_mol_ids: list[str] = Field(
        default_factory=list, description="mol_ids that failed"
    )
    rerun_mol_ids: list[str] = Field(
        default_factory=list, description="mol_ids marked rerun"
    )

    @computed_field
    @property
    def success_rate(self) -> float:
        """Success rate as percentage."""
        if self.total_count == 0:
            return 0.0
        return 100 * self.success_count / self.total_count

    @computed_field
    @property
    def avg_time_per_molecule(self) -> float:
        """Average processing time per molecule in seconds."""
        if self.total_count == 0:
            return 0.0
        return self.total_time / self.total_count

    @computed_field
    @property
    def hof_range(self) -> Optional[float]:
        """Energy range (max - min) in kcal/mol."""
        if self.min_hof is None or self.max_hof is None:
            return None
        return self.max_hof - self.min_hof

    @field_serializer("timestamp_start", mode="plain")
    def serialize_timestamp_start(self, v: datetime) -> str:
        """Serialize non-optional timestamp_start to ISO format."""
        return v.isoformat()

    @field_serializer("timestamp_end", mode="plain")
    def serialize_timestamp_end(self, v: Optional[datetime]) -> Optional[str]:
        """Serialize optional timestamp_end to ISO format or None."""
        return v.isoformat() if v is not None else None


class BatchMolecule(BaseModel):
    """A molecule selected for batch processing.

    Lightweight representation used during batch execution.
    """

    mol_id: str = Field(..., description="Molecule identifier")
    smiles: str = Field(..., description="SMILES string")
    batch_order: int = Field(..., description="Position in batch (1-indexed)")

    # Molecular descriptors (for timeout prediction)
    nheavy: int = Field(..., description="Number of heavy atoms")
    nrotbonds: int = Field(default=0, description="Number of rotatable bonds")

    # Retry tracking
    retry_count: int = Field(default=0, description="Number of previous attempts")


class Batch(BaseModel):
    """A batch of molecules selected for processing.

    Created by BatchCSVManager.select_batch() and consumed by BatchExecutionManager.
    """

    # Identification
    batch_id: str = Field(..., description="Unique batch identifier")

    # Molecules
    molecules: list[BatchMolecule] = Field(
        default_factory=list, description="Molecules in this batch"
    )

    # Configuration (applied to ALL molecules in batch)
    crest_timeout_minutes: int = Field(
        ..., description="CREST timeout in minutes for all molecules"
    )
    mopac_timeout_minutes: int = Field(
        ..., description="MOPAC timeout in minutes for all molecules"
    )

    # Strategy
    sorting_strategy: BatchSortingStrategy = Field(
        default=BatchSortingStrategy.RERUN_FIRST_THEN_EASY,
        description="Strategy used to select molecules",
    )
    failure_policy: BatchFailurePolicy = Field(
        default=BatchFailurePolicy.PARTIAL_OK,
        description="How to handle failures in this batch",
    )

    @computed_field
    @property
    def size(self) -> int:
        """Number of molecules in batch."""
        return len(self.molecules)

    @computed_field
    @property
    def is_empty(self) -> bool:
        """Whether batch has no molecules."""
        return len(self.molecules) == 0


class BatchRowCSV(BaseModel):
    """Schema for a row in the batch tracking CSV file.

    This model defines all 36 columns in the CSV.
    Used for validation and documentation.
    """

    # === Identification (from input) ===
    mol_id: str = Field(..., description="Unique molecule identifier")
    smiles: str = Field(..., description="SMILES string")

    # === Molecular Descriptors (from input) ===
    nheavy: int = Field(..., description="Number of heavy atoms")
    nrotbonds: int = Field(default=0, description="Number of rotatable bonds")
    tpsa: Optional[float] = Field(default=None, description="Topological PSA")
    aromatic_rings: Optional[int] = Field(
        default=None, description="Number of aromatic rings"
    )
    has_heteroatoms: Optional[bool] = Field(
        default=None, description="Whether has heteroatoms"
    )

    # === Reference Data (from input) ===
    reference_hof: Optional[float] = Field(
        default=None, description="CBS-QB3 reference HOF (kcal/mol)"
    )

    # === Batch Status ===
    status: MoleculeStatus = Field(
        default=MoleculeStatus.PENDING, description="Current processing status"
    )
    retry_count: int = Field(default=0, description="Number of retry attempts")
    max_retries: int = Field(default=3, description="Max retries before Skip")

    # === Batch Assignment ===
    batch_id: Optional[str] = Field(default=None, description="Current batch ID")
    batch_order: Optional[int] = Field(
        default=None, description="Position in batch (1-indexed)"
    )
    batch_failure_policy: Optional[BatchFailurePolicy] = Field(
        default=None, description="Failure policy for current batch"
    )

    # === Timeout Configuration (per batch) ===
    assigned_crest_timeout: Optional[float] = Field(
        default=None, description="Assigned CREST timeout (minutes)"
    )
    assigned_mopac_timeout: Optional[float] = Field(
        default=None, description="Assigned MOPAC timeout (minutes)"
    )

    # === CREST Execution Results ===
    crest_status: Optional[str] = Field(default=None, description="CREST status")
    crest_conformers_generated: Optional[int] = Field(
        default=None, description="Conformers from CREST"
    )
    crest_time: Optional[float] = Field(
        default=None, description="CREST execution time (s)"
    )
    crest_error: Optional[str] = Field(default=None, description="CREST error message")

    # === MOPAC Execution Results ===
    num_conformers_selected: Optional[int] = Field(
        default=None, description="Conformers sent to MOPAC"
    )
    most_stable_hof: Optional[float] = Field(
        default=None, description="Best HOF (kcal/mol)"
    )
    quality_grade: Optional[str] = Field(default=None, description="Quality grade")

    # === Delta-E (Energy Spread) ===
    delta_e_12: Optional[float] = Field(
        default=None, description="Energy diff conf1-conf2 (kcal/mol)"
    )
    delta_e_13: Optional[float] = Field(
        default=None, description="Energy diff conf1-conf3 (kcal/mol)"
    )
    delta_e_15: Optional[float] = Field(
        default=None, description="Energy diff conf1-conf5 (kcal/mol)"
    )

    # === Final Status ===
    success: Optional[bool] = Field(default=None, description="Processing succeeded")
    error_message: Optional[str] = Field(
        default=None, description="Error if failed"
    )
    total_execution_time: Optional[float] = Field(
        default=None, description="Total time (s)"
    )

    # === Actual Timeouts Used ===
    actual_crest_timeout_used: Optional[float] = Field(
        default=None, description="Actual CREST timeout (min)"
    )
    actual_mopac_timeout_used: Optional[float] = Field(
        default=None, description="Actual MOPAC timeout (min)"
    )

    # === Timestamps ===
    timestamp: Optional[datetime] = Field(
        default=None, description="Last processing timestamp"
    )

    # === Error Tracking ===
    last_error_message: Optional[str] = Field(
        default=None, description="Most recent error (kept for audit)"
    )

    @field_serializer("timestamp", mode="plain")
    def serialize_timestamp(self, v: Optional[datetime]) -> Optional[str]:
        """Serialize datetime to ISO format."""
        return v.isoformat() if v is not None else None

    model_config = {
        "use_enum_values": True,
    }


class ConformerDetail(BaseModel):
    """Detail information for a single conformer.

    Stored in JSON files per molecule for detailed analysis.
    """

    conformer_index: int = Field(..., description="Conformer index (0-based)")
    energy_hof: Optional[float] = Field(
        default=None, description="Heat of formation (kcal/mol)"
    )
    energy_total: Optional[float] = Field(
        default=None, description="Total energy (eV)"
    )
    mopac_status: str = Field(default="NOT_ATTEMPTED", description="MOPAC status")
    mopac_time: Optional[float] = Field(
        default=None, description="MOPAC execution time (s)"
    )
    mopac_error: Optional[str] = Field(default=None, description="MOPAC error message")
    geometry_file: Optional[str] = Field(
        default=None, description="Path to optimized geometry"
    )


class MoleculeDetail(BaseModel):
    """Complete detail file for a molecule.

    Stored as JSON in data/conformer_details/{mol_id}.json
    """

    # Identification
    mol_id: str = Field(..., description="Molecule identifier")
    smiles: str = Field(..., description="SMILES string")

    # Processing metadata
    batch_id: str = Field(..., description="Batch that processed this molecule")
    timestamp: datetime = Field(
        default_factory=lambda: datetime.now(timezone.utc),
        description="Processing timestamp",
    )

    # CREST results
    crest_status: str = Field(..., description="CREST status")
    crest_conformers_generated: int = Field(
        default=0, description="Total conformers from CREST"
    )
    crest_time: Optional[float] = Field(default=None, description="CREST time (s)")
    crest_error: Optional[str] = Field(default=None, description="CREST error")

    # Conformer details
    conformers: list[ConformerDetail] = Field(
        default_factory=list, description="Per-conformer results"
    )

    # Summary
    num_conformers_selected: int = Field(
        default=0, description="Conformers sent to MOPAC"
    )
    num_conformers_successful: int = Field(
        default=0, description="Conformers with valid HOF"
    )
    most_stable_hof: Optional[float] = Field(
        default=None, description="Best HOF (kcal/mol)"
    )

    # Quality
    quality_grade: str = Field(default="FAILED", description="Quality grade")
    issues: list[str] = Field(default_factory=list, description="Issues found")

    # Final
    success: bool = Field(default=False, description="Processing succeeded")
    error_message: Optional[str] = Field(default=None, description="Error if failed")
    total_execution_time: Optional[float] = Field(
        default=None, description="Total time (s)"
    )

    def to_json_path(self, base_dir: Path) -> Path:
        """Get path for JSON detail file."""
        return base_dir / f"{self.mol_id}.json"

    @field_serializer("timestamp", mode="plain")
    def serialize_timestamp(self, v: datetime) -> str:
        """Serialize datetime to ISO format."""
        return v.isoformat()
