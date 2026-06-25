"""Pydantic models for batch processing of CREST PM7 pipeline.

This module defines data models for:
- BatchResult: Summary statistics from batch execution
- Batch: Selected molecules for processing
- BatchRowCSV: Schema for CSV tracking file (58 columns)
- ConformerDetail: Per-molecule JSON detail file
"""

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

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
    failed_count: int = Field(default=0, description="Molecules that failed processing")

    # Timing (in seconds)
    total_time: float = Field(default=0.0, description="Total execution time")
    timestamp_start: datetime = Field(
        default_factory=lambda: datetime.now(timezone.utc)
    )
    timestamp_end: datetime | None = Field(default=None)

    # Energy statistics (from successful molecules in THIS BATCH only)
    min_hof: float | None = Field(default=None, description="Minimum HOF in batch")
    max_hof: float | None = Field(default=None, description="Maximum HOF in batch")
    min_hof_mol_id: str | None = Field(default=None, description="mol_id with min HOF")
    max_hof_mol_id: str | None = Field(default=None, description="mol_id with max HOF")

    # Issues tracking (for debugging)
    failed_mol_ids: list[str] = Field(
        default_factory=list, description="mol_ids that failed"
    )
    rerun_mol_ids: list[str] = Field(
        default_factory=list, description="mol_ids marked rerun"
    )

    @computed_field  # type: ignore[prop-decorator]
    @property
    def success_rate(self) -> float:
        """Success rate as percentage."""
        if self.total_count == 0:
            return 0.0
        return 100 * self.success_count / self.total_count

    @computed_field  # type: ignore[prop-decorator]
    @property
    def avg_time_per_molecule(self) -> float:
        """Average processing time per molecule in seconds."""
        if self.total_count == 0:
            return 0.0
        return self.total_time / self.total_count

    @computed_field  # type: ignore[prop-decorator]
    @property
    def hof_range(self) -> float | None:
        """Energy range (max - min) in kcal/mol."""
        if self.min_hof is None or self.max_hof is None:
            return None
        return self.max_hof - self.min_hof

    @field_serializer("timestamp_start", mode="plain")
    def serialize_timestamp_start(self, v: datetime) -> str:
        """Serialize non-optional timestamp_start to ISO format."""
        return v.isoformat()

    @field_serializer("timestamp_end", mode="plain")
    def serialize_timestamp_end(self, v: datetime | None) -> str | None:
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
    reruns: int = Field(default=0, description="Number of previous attempts")


class Batch(BaseModel):
    """A batch of molecules selected for processing.

    Created by BatchCSVManager.select_batch() and consumed by
    BatchExecutionManager.
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
    method_id: str = Field(
        default="pm7_delta_learning",
        description="Calculation method registry identifier used for the batch",
    )
    method_version: str = Field(
        default="1.0.0",
        description="Calculation method version used for the batch",
    )
    method_snapshot: dict[str, Any] = Field(
        default_factory=dict,
        description="Validated calculation method definition snapshot",
    )

    @computed_field  # type: ignore[prop-decorator]
    @property
    def size(self) -> int:
        """Number of molecules in batch."""
        return len(self.molecules)

    @computed_field  # type: ignore[prop-decorator]
    @property
    def is_empty(self) -> bool:
        """Whether batch has no molecules."""
        return len(self.molecules) == 0


class BatchRowCSV(BaseModel):
    """Schema for a row in the batch tracking CSV file (58 columns).

    This model defines the CSV columns matching the target schema.
    Used for validation and documentation.

    Column groups:
    - Identity (3): mol_id, status, smiles
    - Molecular Properties (7): multiplicity, charge, nheavy, H298_cbs,
      H298_pm7, target_delta_kcalmol, abs_diff_%
    - Batch Info (4): batch_id, timestamp, total_time, reruns
    - RDKit Descriptors (15): rdkit_nrotbonds, rdkit_tpsa, rdkit_num_rings,
      rdkit_fsp3, rdkit_mol_weight, rdkit_hbond_donors, rdkit_hbond_acceptors,
      rdkit_nC, rdkit_nH, rdkit_nO, rdkit_nN, rdkit_bonds_single,
      rdkit_bonds_double, rdkit_bonds_triple, rdkit_bonds_aromatic
    - CREST (12): crest_status, xtb, v3, qm, nci, c_method, energy_window,
      rmsd_threshold, opt_lvl, crest_conformers_generated,
      num_conformers_selected, crest_time_s
    - MOPAC (13): mopac_status, k_selected_pm7, mopac_dipole_debye,
      mopac_ionization_potential_ev, mopac_homo_ev, mopac_lumo_ev,
      mopac_gap_ev, mopac_cosmo_area_a2, mopac_cosmo_volume_a3,
      mopac_gradient_norm, mopac_num_scf_cycles, mopac_point_group,
      mopac_time_s
    - Batch Config (4): batch_order, batch_failure_policy,
      assigned_crest_timeout, assigned_mopac_timeout
    Total: 3+7+4+15+12+13+4 = 58 columns
    """

    # === Identity (3) ===
    mol_id: str = Field(..., description="Unique molecule identifier")
    status: MoleculeStatus = Field(
        default=MoleculeStatus.PENDING, description="Current processing status"
    )
    smiles: str = Field(..., description="SMILES string")

    # === Molecular Properties (7) ===
    multiplicity: int = Field(default=1, description="Spin multiplicity")
    charge: int = Field(default=0, description="Molecular charge")
    nheavy: int = Field(..., description="Number of heavy atoms")
    H298_cbs: float | None = Field(
        default=None, description="CBS-QB3 reference enthalpy (kcal/mol)"
    )
    H298_pm7: float | None = Field(default=None, description="PM7 enthalpy (kcal/mol)")
    target_delta_kcalmol: float | None = Field(
        default=None, description="Signed delta: H298_cbs - H298_pm7"
    )
    abs_diff_pct: float | None = Field(
        default=None, alias="abs_diff_%", description="Percentage difference"
    )
    cbs_quality_flag: str = Field(
        default="", description='CBS reference quality: "OK", "SUSPECT", or ""'
    )

    # === Batch Info (4) ===
    batch_id: str | None = Field(default=None, description="Current batch ID")
    timestamp: str | None = Field(default=None, description="Processing timestamp")
    total_time: float | None = Field(default=None, description="Total time (min)")
    reruns: int = Field(default=0, description="Number of reruns")

    # === RDKit Descriptors (15) ===
    rdkit_nrotbonds: int | None = Field(default=None, description="Rotatable bonds")
    rdkit_tpsa: float | None = Field(default=None, description="Topological PSA")
    rdkit_num_rings: int | None = Field(default=None, description="Aromatic rings")
    rdkit_fsp3: float | None = Field(default=None, description="Fraction sp3 carbons")
    rdkit_mol_weight: float | None = Field(default=None, description="Molecular weight")
    rdkit_hbond_donors: int | None = Field(default=None, description="H-bond donors")
    rdkit_hbond_acceptors: int | None = Field(
        default=None, description="H-bond acceptors"
    )
    rdkit_nC: int | None = Field(default=None, description="Carbon atom count")
    rdkit_nH: int | None = Field(default=None, description="Hydrogen atom count")
    rdkit_nO: int | None = Field(default=None, description="Oxygen atom count")
    rdkit_nN: int | None = Field(default=None, description="Nitrogen atom count")
    rdkit_bonds_single: int | None = Field(default=None, description="Single bonds")
    rdkit_bonds_double: int | None = Field(default=None, description="Double bonds")
    rdkit_bonds_triple: int | None = Field(default=None, description="Triple bonds")
    rdkit_bonds_aromatic: int | None = Field(default=None, description="Aromatic bonds")

    # === CREST (12) ===
    crest_status: str | None = Field(default=None, description="CREST status")
    xtb: bool | None = Field(default=None, description="xTB pre-optimization")
    v3: bool | None = Field(default=None, description="CREST v3 algorithm")
    qm: bool | None = Field(default=None, description="CREST quick mode enabled")
    nci: bool | None = Field(default=None, description="CREST NCI mode")
    c_method: str | None = Field(default=None, description="CREST method")
    energy_window: float | None = Field(
        default=None, description="CREST energy window (kcal/mol)"
    )
    rmsd_threshold: float | None = Field(
        default=None, description="CREST RMSD threshold (A)"
    )
    opt_lvl: int | None = Field(
        default=None, description="CREST optimization level (0-2)"
    )
    crest_conformers_generated: int | None = Field(
        default=None, description="Conformers from CREST"
    )
    num_conformers_selected: int | None = Field(
        default=None, description="Conformers sent to MOPAC"
    )
    crest_time_s: float | None = Field(
        default=None, description="CREST execution time (s)"
    )

    # === MOPAC (13) ===
    mopac_status: str | None = Field(default=None, description="MOPAC status")
    k_selected_pm7: int | None = Field(
        default=None, description="CREST rank of PM7-selected conformer"
    )
    mopac_dipole_debye: float | None = Field(default=None, description="Dipole (Debye)")
    mopac_ionization_potential_ev: float | None = Field(
        default=None, description="Ionization potential (eV)"
    )
    mopac_homo_ev: float | None = Field(default=None, description="HOMO energy (eV)")
    mopac_lumo_ev: float | None = Field(default=None, description="LUMO energy (eV)")
    mopac_gap_ev: float | None = Field(default=None, description="HOMO-LUMO gap (eV)")
    mopac_cosmo_area_a2: float | None = Field(
        default=None, description="COSMO area (A^2)"
    )
    mopac_cosmo_volume_a3: float | None = Field(
        default=None, description="COSMO volume (A^3)"
    )
    mopac_gradient_norm: float | None = Field(
        default=None, description="Gradient norm (kcal/mol/A)"
    )
    mopac_num_scf_cycles: int | None = Field(
        default=None, description="SCF cycle count"
    )
    mopac_point_group: str | None = Field(default=None, description="Point group")
    mopac_time_s: float | None = Field(
        default=None, description="MOPAC execution time (s)"
    )

    # === Batch Config (4) ===
    batch_order: int | None = Field(
        default=None, description="Position in batch (1-indexed)"
    )
    batch_failure_policy: BatchFailurePolicy | None = Field(
        default=None, description="Failure policy for current batch"
    )
    assigned_crest_timeout: float | None = Field(
        default=None, description="Assigned CREST timeout (minutes)"
    )
    assigned_mopac_timeout: float | None = Field(
        default=None, description="Assigned MOPAC timeout (minutes)"
    )

    model_config = {
        "use_enum_values": True,
        "populate_by_name": True,
    }


class ConformerDetail(BaseModel):
    """Detail information for a single conformer.

    Stored in JSON files per molecule for detailed analysis.
    """

    conformer_index: int = Field(..., description="Conformer index (0-based)")
    energy_hof: float | None = Field(
        default=None, description="Heat of formation (kcal/mol)"
    )
    energy_total: float | None = Field(default=None, description="Total energy (eV)")
    mopac_status: str = Field(default="NOT_ATTEMPTED", description="MOPAC status")
    mopac_time: float | None = Field(
        default=None, description="MOPAC execution time (s)"
    )
    mopac_error: str | None = Field(default=None, description="MOPAC error message")
    geometry_file: str | None = Field(
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
    crest_time: float | None = Field(default=None, description="CREST time (s)")
    crest_error: str | None = Field(default=None, description="CREST error")

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
    most_stable_hof: float | None = Field(
        default=None, description="Best HOF (kcal/mol)"
    )

    # Quality
    quality_grade: str = Field(default="FAILED", description="Quality grade")
    issues: list[str] = Field(default_factory=list, description="Issues found")

    # Final
    success: bool = Field(default=False, description="Processing succeeded")
    error_message: str | None = Field(default=None, description="Error if failed")
    total_execution_time: float | None = Field(
        default=None, description="Total time (s)"
    )

    def to_json_path(self, base_dir: Path) -> Path:
        """Get path for JSON detail file."""
        return base_dir / f"{self.mol_id}.json"

    @field_serializer("timestamp", mode="plain")
    def serialize_timestamp(self, v: datetime) -> str:
        """Serialize datetime to ISO format."""
        return v.isoformat()
