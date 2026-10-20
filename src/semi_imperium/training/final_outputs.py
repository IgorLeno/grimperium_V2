"""Durable outputs of the conformer that was finally chosen for a molecule.

The MOPAC working directory is disposable; what survives is, per molecule,
one folder holding only the evidence of the chosen conformer: its optimized
geometry (``.xyz``), the MOPAC optimization and verification outputs
(``.out`` / ``.freq.out``) and a small ``.json`` naming the conformer, its
heat of formation and the lowest frequency. Other conformers are not kept.
"""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

from semi_imperium.conformers import ConformerGeometry
from semi_imperium.domain import CalculationResultData
from semi_imperium.training.geometry_store import GeometryRecord

#: Folder inside the store that holds the kept outputs, one subfolder per molecule.
OUTPUT_DIRNAME = "geometrias"

_SAFE_MOL_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class FinalOutputError(RuntimeError):
    """The outputs of the chosen conformer could not be kept."""


def export_final_outputs(
    *,
    store_root: Path,
    mol_id: str,
    smiles: str,
    hamiltonian: str,
    protocol: str,
    calculation_id: str,
    signature: str,
    geometry: ConformerGeometry,
    result: CalculationResultData,
) -> Path:
    """Write the kept files of one molecule and return their folder.

    ``result.artifact_paths`` holds, relative to ``store_root``, the MOPAC
    optimization output first and the verification output second. The folder
    is built beside its final place and renamed into it, so a crash never
    leaves a half-written molecule; a second export replaces the first.

    Raises:
        FinalOutputError: If ``mol_id`` is unsafe as a folder name, a MOPAC
            output is missing, or a file cannot be written.
    """
    if not _SAFE_MOL_ID.fullmatch(mol_id):
        raise FinalOutputError(f"unsafe mol_id for a folder name: {mol_id!r}")
    if len(result.artifact_paths) < 2:
        raise FinalOutputError(
            "expected the MOPAC optimization and verification outputs, got "
            f"{len(result.artifact_paths)} path(s)"
        )
    sources = {
        f"{mol_id}.out": store_root / result.artifact_paths[0],
        f"{mol_id}.freq.out": store_root / result.artifact_paths[1],
    }
    missing = [str(path) for path in sources.values() if not path.is_file()]
    if missing:
        raise FinalOutputError(f"MOPAC output not found: {', '.join(missing)}")

    root = store_root / OUTPUT_DIRNAME
    final = root / mol_id
    staging = root / f"{mol_id}.tmp"
    try:
        shutil.rmtree(staging, ignore_errors=True)
        staging.mkdir(parents=True)
        for name, source in sources.items():
            shutil.copyfile(source, staging / name)
        xyz = GeometryRecord(
            mol_id=mol_id,
            smiles=smiles,
            hamiltonian=hamiltonian,
            protocol=protocol,
            signature=signature,
            calculation_id=calculation_id,
            conformer_index=result.conformer_index,
            elements=geometry.elements,
            coordinates=geometry.coordinates,
            hof_kcal_mol=result.energy_hof_kcal_mol,
            lowest_frequency_cm1=result.lowest_frequency_cm1,
            minimum_confirmed=True,
        ).to_xyz()
        (staging / f"{mol_id}.xyz").write_text(xyz, encoding="utf-8")
        summary = {
            "mol_id": mol_id,
            "smiles": smiles,
            "hamiltonian": hamiltonian,
            "protocol": protocol,
            "calculation_id": calculation_id,
            "signature": signature,
            "conformer_index": result.conformer_index,
            "conformers_evaluated": result.conformers_evaluated,
            "hof_kcal_mol": result.energy_hof_kcal_mol,
            "lowest_frequency_cm1": result.lowest_frequency_cm1,
            "minimum_confirmed": True,
            "files": sorted(path.name for path in staging.iterdir())
            + [f"{mol_id}.json"],
        }
        (staging / f"{mol_id}.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        shutil.rmtree(final, ignore_errors=True)
        staging.rename(final)
    except OSError as exc:
        shutil.rmtree(staging, ignore_errors=True)
        raise FinalOutputError(f"could not write outputs of {mol_id}: {exc}") from exc
    return final


__all__ = ["OUTPUT_DIRNAME", "FinalOutputError", "export_final_outputs"]
