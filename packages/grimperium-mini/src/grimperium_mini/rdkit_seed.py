"""Initial 3D geometry generation from SMILES."""

from __future__ import annotations

from pathlib import Path


def generate_initial_xyz(smiles: str, output_xyz: Path, name: str = "molecule") -> Path:
    """Generate an initial 3D XYZ geometry using RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:  # pragma: no cover - depends on optional environment
        raise RuntimeError(
            "RDKit is required for non-dry-run geometry generation"
        ) from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol, randomSeed=0xF00D)  # type: ignore[attr-defined]
    if status != 0:
        status = AllChem.EmbedMolecule(  # type: ignore[attr-defined]
            mol, useRandomCoords=True, randomSeed=0xF00D
        )
    if status != 0:
        raise ValueError(f"RDKit failed to embed SMILES: {smiles}")
    AllChem.UFFOptimizeMolecule(mol, maxIters=500)  # type: ignore[attr-defined]
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), name]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol()} {pos.x:.8f} {pos.y:.8f} {pos.z:.8f}")
    output_xyz.parent.mkdir(parents=True, exist_ok=True)
    output_xyz.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_xyz
