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
  S/Cp/G properties need;
* the source mixes units in the enthalpy columns: about 1% of the
  ``H298_cbs``/``H298_b3``/``cbs_b3`` values are written as plain integers
  in cal/mol (e.g. ``-79523`` for 3-butenoic acid) while all others carry
  decimals in kcal/mol. Integer-formatted values are divided by 1000 and
  the repair is recorded per row in ``h298_unit_repair``. The identity
  ``H298_cbs = H298_b3 + cbs_b3`` holds for every source row only after
  this repair, which is checked here, so a wrong reading fails the build.

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
ENTHALPY_COLUMNS = ("H298_cbs", "H298_b3", "cbs_b3")
UNIT_REPAIR_COLUMN = "h298_unit_repair"
CAL_PER_KCAL = 1000.0
#: Tolerance of H298_cbs = H298_b3 + cbs_b3 (kcal/mol); the integer
#: cal/mol values are rounded to 1 cal/mol, i.e. 0.001 kcal/mol each.
IDENTITY_TOLERANCE = 0.005


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


def is_integer_text(text: str) -> bool:
    """Whether a numeric field was written without decimals or exponent."""
    return not any(mark in text for mark in ".eE")


def repair_enthalpy_units(raw: pd.DataFrame) -> pd.DataFrame:
    """Convert integer-formatted cal/mol enthalpies to kcal/mol.

    ``raw`` must hold the enthalpy columns as text. Returns a copy with
    float columns and ``h298_unit_repair`` naming the repaired columns
    (``;``-separated, empty when none).

    Raises:
        ValueError: If H298_cbs = H298_b3 + cbs_b3 fails after the repair.
    """
    out = raw.copy()
    repaired = pd.DataFrame(index=raw.index)
    for column in ENTHALPY_COLUMNS:
        text = raw[column].astype(str).str.strip()
        integer = text.map(is_integer_text)
        values = text.astype(float)
        out[column] = values.where(~integer, values / CAL_PER_KCAL)
        repaired[column] = integer
    out[UNIT_REPAIR_COLUMN] = repaired.apply(
        lambda row: ";".join(column for column in ENTHALPY_COLUMNS if row[column]),
        axis=1,
    )
    residual = (out["H298_cbs"] - out["H298_b3"] - out["cbs_b3"]).abs()
    broken = residual > IDENTITY_TOLERANCE
    if broken.any():
        sample = out.loc[broken, list(ENTHALPY_COLUMNS)].head(5).to_dict("records")
        raise ValueError(
            f"{int(broken.sum())} rows violate H298_cbs = H298_b3 + cbs_b3 after "
            f"the unit repair, e.g. {sample}"
        )
    return out


def build(source: Path, legacy: Path | None) -> pd.DataFrame:
    raw = pd.read_csv(source, dtype={column: str for column in ENTHALPY_COLUMNS})
    raw = raw.rename(columns={raw.columns[0]: SOURCE_INDEX_COLUMN})
    raw = repair_enthalpy_units(raw)
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
        "h298_unit_repair": {
            column: int(table[UNIT_REPAIR_COLUMN].str.contains(column).sum())
            for column in ENTHALPY_COLUMNS
        },
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
        "README does not state units. Integer-formatted enthalpies in the "
        "source are cal/mol and were divided by 1000 (see h298_unit_repair).",
        "output_sha256": sha256_of(args.output),
    }
    manifest_path = args.output.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    print(json.dumps({k: manifest[k] for k in ("rows", "unique_smiles")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
