"""Small data models used by the mini pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

SemiempiricalMethod = Literal["AM1", "PM3", "PM7"]

KCAL_TO_KJ = 4.184


@dataclass
class MoleculeInput:
    mol_id: str
    molecule: str
    smiles: str
    formula: str = ""
    group: str = ""
    hf_exp_kJmol: float | None = None
    hf_tcc_am1_kJmol: float | None = None
    hf_tcc_pm3_kJmol: float | None = None
    hf_tcc_pm7_kJmol: float | None = None
    hf_joback_kJmol: float | None = None
    hf_constantinou_gani_kJmol: float | None = None


@dataclass
class PipelineTask:
    mol_id: str
    molecule: str
    smiles: str
    method: SemiempiricalMethod
    run_crest: bool = True
    charge: int = 0
    multiplicity: int = 1
    formula: str = ""
    group: str = ""
    max_conformers_to_process: int | None = None
    hf_exp_kJmol: float | None = None
    hf_tcc_am1_kJmol: float | None = None
    hf_tcc_pm3_kJmol: float | None = None
    hf_tcc_pm7_kJmol: float | None = None
    hf_joback_kJmol: float | None = None
    hf_constantinou_gani_kJmol: float | None = None


@dataclass
class CrestResult:
    status: str
    conformer_files: list[Path]
    n_conformers_found: int
    work_dir: Path
    error_message: str = ""
    command: list[str] | None = None


@dataclass
class MopacResult:
    status: str
    hof_kcalmol: float | None = None
    hof_kJmol: float | None = None
    output_file: Path | None = None
    input_file: Path | None = None
    error_message: str = ""
    command: list[str] | None = None


@dataclass
class ConformerResult:
    mol_id: str
    molecule: str
    method: SemiempiricalMethod
    conformer_rank: int
    conformer_file: Path
    mopac_status: str
    hof_kcalmol: float | None
    hof_kJmol: float | None
    mopac_output_file: Path | None
    error_message: str = ""
