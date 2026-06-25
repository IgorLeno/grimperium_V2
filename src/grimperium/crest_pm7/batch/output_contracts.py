"""Output contracts for the PR6 batch state/result split."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from grimperium.calculation.contracts.models import MoleculeCalculationResult
from grimperium.calculation.output.csv_writer import (
    UnitSelection,
    estimate_rows,
    write_canonical_csv,
)

BATCH_STATE_FILENAME = "batch_state.csv"
CALCULATION_RESULTS_CSV_FILENAME = "calculation_results.csv"
CALCULATION_RESULTS_XLSX_FILENAME = "calculation_results.xlsx"

BATCH_STATE_COLUMNS = [
    "mol_id",
    "status",
    "smiles",
    "multiplicity",
    "charge",
    "nheavy",
    "batch_id",
    "batch_order",
    "batch_failure_policy",
    "timestamp",
    "total_time",
    "reruns",
    "crest_status",
    "xtb_status",
    "mopac_status",
    "assigned_crest_timeout",
    "assigned_mopac_timeout",
    "assigned_worker",
    "worker_status",
    "assigned_at",
    "method_id",
    "method_version",
    "method_definition_snapshot",
]

SCIENTIFIC_RESULT_COLUMNS = [
    "H298_cbs",
    "H298_pm7",
    "H298_predicted",
    "target_delta_kcalmol",
    "delta_correction",
    "abs_diff_%",
    "k_selected_pm7",
    "rdkit_nrotbonds",
    "rdkit_tpsa",
    "rdkit_num_rings",
    "rdkit_fsp3",
    "rdkit_mol_weight",
    "rdkit_hbond_donors",
    "rdkit_hbond_acceptors",
    "rdkit_nC",
    "rdkit_nH",
    "rdkit_nO",
    "rdkit_nN",
    "rdkit_bonds_single",
    "rdkit_bonds_double",
    "rdkit_bonds_triple",
    "rdkit_bonds_aromatic",
    "mopac_dipole_debye",
    "mopac_ionization_potential_ev",
    "mopac_homo_ev",
    "mopac_lumo_ev",
    "mopac_gap_ev",
    "mopac_cosmo_area_a2",
    "mopac_cosmo_volume_a3",
    "mopac_gradient_norm",
    "mopac_num_scf_cycles",
    "mopac_point_group",
    "mopac_time_s",
]


@dataclass(frozen=True)
class BatchOutputLayout:
    """Filesystem layout for split batch outputs."""

    output_dir: Path

    def __post_init__(self) -> None:
        object.__setattr__(self, "output_dir", Path(self.output_dir))

    @property
    def batch_state_csv(self) -> Path:
        """Path for operational batch state."""
        return self.output_dir / BATCH_STATE_FILENAME

    @property
    def calculation_results_csv(self) -> Path:
        """Path for canonical scientific results in CSV format."""
        return self.output_dir / CALCULATION_RESULTS_CSV_FILENAME

    @property
    def calculation_results_xlsx(self) -> Path:
        """Path for canonical scientific results in XLSX format."""
        return self.output_dir / CALCULATION_RESULTS_XLSX_FILENAME


@dataclass(frozen=True)
class BatchResultWriteReport:
    """Summary of canonical batch result files written."""

    csv_path: Path
    xlsx_path: Path | None
    result_count: int
    estimate_row_count: int


def write_batch_calculation_results(
    results: list[MoleculeCalculationResult],
    layout: BatchOutputLayout,
    *,
    units: UnitSelection = "kcal/mol",
    include_xlsx: bool = True,
) -> BatchResultWriteReport:
    """Write canonical scientific batch results to the split result files."""
    write_canonical_csv(results, layout.calculation_results_csv, units=units)

    xlsx_path = None
    if include_xlsx:
        from grimperium.calculation.output.xlsx_writer import write_canonical_xlsx

        write_canonical_xlsx(results, layout.calculation_results_xlsx, units=units)
        xlsx_path = layout.calculation_results_xlsx

    return BatchResultWriteReport(
        csv_path=layout.calculation_results_csv,
        xlsx_path=xlsx_path,
        result_count=len(results),
        estimate_row_count=len(estimate_rows(results, units=units)),
    )
