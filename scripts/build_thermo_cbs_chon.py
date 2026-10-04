#!/usr/bin/env python3
"""Build the CHON reference cut of ThermoCBS (Chemperium) with geometries.

Source: ThermoCBS, Dobbelaere et al., J. Cheminf. 16:99 (2024),
Zenodo 10.5281/zenodo.11409710, CC-BY 4.0.

Why this script exists:

* the previous cut (``thermo_cbs_chon.csv``) filtered on SMILES text and
  let 1808 molecules with aromatic sulfur, selenium or boron through;
  elements are now read from the reference ``xyz`` block instead;
* the previous cut dropped the reference geometry, entropy and heat
  capacities, which the conformer-selection validation and the planned
  S/Cp/G properties need.

Rows are never deduplicated here. Some SMILES appear more than once
(distinct reference conformers); they are flagged so the training split
can decide explicitly instead of inheriting a silent choice.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ALLOWED_ELEMENTS = frozenset({"C", "H", "O", "N"})
SOURCE_INDEX_COLUMN = "source_index"


def xyz_elements(xyz_block: str) -> frozenset[str]:
    """Return element symbols of an XYZ block (count line, comment, atoms)."""
    lines = xyz_block.split("\n")
    return frozenset(line.split()[0] for line in lines[2:] if line.strip())


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build(source: Path, legacy: Path | None) -> pd.DataFrame:
    raw = pd.read_csv(source)
    raw = raw.rename(columns={raw.columns[0]: SOURCE_INDEX_COLUMN})
    if not raw[SOURCE_INDEX_COLUMN].is_unique:
        raise ValueError(f"{source}: first column is not a unique row index")

    elements = raw["xyz"].map(xyz_elements)
    if (elements.map(len) == 0).any():
        raise ValueError(f"{source}: rows with empty xyz geometry")
    chon = raw[elements.map(lambda found: found <= ALLOWED_ELEMENTS)].copy()

    # Stable ID tied to the source row, not to the filter order.
    chon.insert(0, "mol_id", chon[SOURCE_INDEX_COLUMN].map(lambda i: f"cbs_{i:05d}"))

    group = chon.groupby("smiles")["H298_cbs"]
    chon["smiles_count"] = group.transform("size").astype(int)
    chon["conformer_rank_by_h298"] = group.rank(method="first", ascending=True).astype(
        int
    )

    if legacy is not None:
        legacy_ids = (
            pd.read_csv(legacy, usecols=["mol_id", "smiles"])
            .drop_duplicates("smiles")
            .rename(columns={"mol_id": "legacy_mol_id"})
        )
        chon = chon.merge(legacy_ids, on="smiles", how="left")

    return chon


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, required=True, help="thermo_cbs.csv")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--legacy",
        type=Path,
        default=None,
        help="previous thermo_cbs_chon.csv, to record legacy_mol_id",
    )
    args = parser.parse_args(argv)

    if args.output.exists():
        print(f"refusing to overwrite {args.output}", file=sys.stderr)
        return 1

    table = build(args.source, args.legacy)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": str(args.source),
        "source_sha256": sha256_of(args.source),
        "source_citation": "Dobbelaere et al., J. Cheminf. 16:99 (2024); "
        "Zenodo 10.5281/zenodo.11409710; CC-BY 4.0",
        "filter": "elements from xyz subset of {C,H,O,N}",
        "rows": len(table),
        "unique_smiles": int(table["smiles"].nunique()),
        "multiplicity": table["multiplicity"].value_counts().to_dict(),
        "charge": table["charge"].value_counts().to_dict(),
        "units": {
            "H298_cbs": "kcal/mol",
            "H298_b3": "kcal/mol",
            "cbs_b3": "kcal/mol (H298_cbs - H298_b3)",
            "S298": "cal/(mol K)",
            "cp_*": "cal/(mol K)",
        },
        "units_note": "inferred from CH4/CH3OH sanity values; the source "
        "README does not state units",
        "output_sha256": sha256_of(args.output),
    }
    manifest_path = args.output.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    print(json.dumps({k: manifest[k] for k in ("rows", "unique_smiles")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
