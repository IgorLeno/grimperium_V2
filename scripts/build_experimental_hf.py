#!/usr/bin/env python3
"""Build the experimental gas-phase enthalpy-of-formation table (298.15 K).

Sources, in selection priority:

1. ATcT, Active Thermochemical Tables (Ruscic & Bross, Argonne), main table
   of a TN version saved as HTML. The SMILES of each species is read from
   the ``alt`` text of its structure image; the phase from the formula label
   (``(g)``, ``(cr,l)``, ``(aq)``...). Gas species with an ATcT ID suffix
   other than ``*0`` are conformer or spin-state variants (``g, syn``,
   ``g, triplet``) and are rejected as ``gas_variant``; ``*0`` is the
   reference state. No explicit licence is stated on
   the site, so the raw page stays out of git.
2. RMG-database reference set (``input/reference_sets/main``, MIT header),
   one YAML per species with H298, uncertainty and the source tag
   (ATcT, CCCBDB, ...).
3. Bains, Petkowski, Zhan & Seager, Data 7:33 (2022), Zenodo
   10.5281/zenodo.4661783 (CC-BY-SA), ``Measured_Enthalpy_V2.7.xlsx``:
   one column per original compilation, in kJ/mol, no uncertainty.

Why the rules below:

* the scope matches the CBS reference: C-containing, only C/H/N/O, neutral,
  closed-shell (no radical electrons), no isotopic labels, gas phase;
* every candidate value is kept in a long table (one row per source value),
  and one value per InChIKey is selected by source priority, then lowest
  reported uncertainty, so the choice is traceable;
* the Bains ``Stewart`` column is the reference set used to parametrise
  PM6/PM7, and ``Yaws`` mixes measured and estimated values: both get the
  lowest priority, and ``stewart_reference`` flags molecules whose PM7 error
  is not an independent test;
* ATcT ions often carry a neutral SMILES in the image ``alt`` text; the net
  charge is read from the row id and charged species are rejected;
* species of one source that collapse to one InChIKey (stereo missing from
  the SMILES, tautomers merged by the standard InChI) keep only the lowest
  H298 of that source, counted in ``isomers_collapsed``;
* Bains rows with ``FILTER == 1`` are kept in the long table but never
  selected;
* a selected value conflicts with another source when they differ by more
  than ``max(2 * combined uncertainty, 1 kcal/mol)``; unknown uncertainty
  counts as zero.

Units: sources are in kJ/mol; outputs are kcal/mol (thermochemical calorie,
``KJ_PER_KCAL`` = 4.184, the same value as ``KCAL_TO_KJ`` in
``grimperium.cli.views.calc_view``).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections.abc import Iterable, Iterator
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

KJ_PER_KCAL = 4.184
ALLOWED_ELEMENTS = frozenset({"C", "H", "N", "O"})
CONFLICT_FLOOR_KCAL = 1.0
GAS_REFERENCE_SUFFIX = "*0"
#: Formula button text ends with the phase label, e.g. ``CH4  (g)``.
ATCT_PHASE_LABEL = re.compile(r"\(([^()]*)\)\s*$")

#: Lower rank wins. RMG values are ranked by the source they cite.
SOURCE_RANK: dict[str, int] = {
    "atct": 0,
    "rmg:ATcT": 1,
    "rmg:CCCBDB": 2,
    "rmg:Cioslowski": 3,
    "bains:ATcT": 10,
    "bains:Cioslowski": 11,
    "bains:Pedley": 12,
    **{f"bains:Cox-{i}": 13 for i in range(1, 6)},
    **{f"bains:NIST-{i}": 14 for i in range(1, 6)},
    "bains:JANAF": 15,
    "bains:Sandler": 16,
    "bains:Jorgensen": 16,
    **{f"bains:Add Lit {i}": 17 for i in range(1, 6)},
    "bains:Winget": 18,
    "bains:Yaws": 90,
    "bains:Stewart": 91,
}
LOW_PRIORITY_RANK = 90
UNRANKED = 50

BAINS_SHEET = "Measured Enthalpy"
BAINS_HEADER_ROW = 2  # zero-based row with Name, SMILES, source columns
BAINS_FILTER_COLUMN = "FILTER"

#: Row id, e.g. ``s1_6n4_1c0 i23 CAS74-82-8``: the ``c<n>`` suffix of the first
#: token is the net charge (the image ``alt`` SMILES of ions is often neutral).
ATCT_ROW_ID = re.compile(
    r"^s\S*?c(?P<charge>-?\d+)\s+i(?P<number>\d+)\s+CAS(?P<cas>\S*)$"
)

CITATIONS = {
    "atct": "Ruscic, B.; Bross, D. H. Active Thermochemical Tables (ATcT) "
    "values based on ver. {version} of the Thermochemical Network, Argonne "
    "National Laboratory (atct.anl.gov)",
    "rmg": "RMG-database input/reference_sets/main, commit {commit} "
    "(github.com/ReactionMechanismGenerator/RMG-database)",
    "bains": "Bains, Petkowski, Zhan, Seager, Data 7:33 (2022); "
    "Zenodo 10.5281/zenodo.4661783 (CC-BY-SA)",
}


@dataclass(frozen=True)
class SourceValue:
    """One experimental value as written in a source, before selection."""

    source: str
    source_id: str
    name: str
    raw_smiles: str
    h298_kj: float
    uncertainty_kj: float | None
    phase: str
    selectable: bool = True
    declared_charge: int | None = None


@dataclass(frozen=True)
class Structure:
    smiles: str
    inchikey: str
    formula: str
    charge: int
    multiplicity: int
    nheavy: int


def kj_to_kcal(value: float) -> float:
    return value / KJ_PER_KCAL


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inchikey_of(mol: Chem.Mol) -> str:
    key: str = Chem.MolToInchiKey(mol)  # type: ignore[no-untyped-call]
    return key


def parse_uncertainty(text: str) -> float | None:
    """ATcT uncertainty cell: ``± 0.63`` or ``exact`` (elements)."""
    cleaned = text.replace("±", "").strip()
    if not cleaned:
        return None
    if cleaned.lower() == "exact":
        return 0.0
    return float(cleaned)


def standardize(smiles: str) -> tuple[Structure | None, str]:
    """Parse and screen one SMILES; return the structure or a reject reason.

    Accepted: RDKit-parsable, contains C, only C/H/N/O, no isotope label,
    total formal charge 0 and no radical electrons (multiplicity 1).
    """
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None, "unparsable_smiles"
    atoms = list(mol.GetAtoms())
    elements = {atom.GetSymbol() for atom in atoms}
    if "C" not in elements or not elements <= ALLOWED_ELEMENTS:
        return None, "elements"
    if any(atom.GetIsotope() for atom in atoms):
        return None, "isotope"
    charge = Chem.GetFormalCharge(mol)
    if charge != 0:
        return None, "charged"
    if sum(atom.GetNumRadicalElectrons() for atom in atoms) != 0:
        return None, "open_shell"
    return (
        Structure(
            smiles=Chem.MolToSmiles(mol),
            inchikey=inchikey_of(mol),
            formula=CalcMolFormula(mol),
            charge=charge,
            multiplicity=1,
            nheavy=mol.GetNumHeavyAtoms(),
        ),
        "",
    )


# --------------------------------------------------------------------- ATcT


class _AtctTableParser(HTMLParser):
    """Collect species rows (``<tr id="s.. i<number> CAS..">``) of the table."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.rows: list[dict[str, str]] = []
        self._row: dict[str, str] | None = None
        self._field: str | None = None
        self._in_button = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = {key: value or "" for key, value in attrs}
        if tag == "tr":
            match = ATCT_ROW_ID.match(attributes.get("id", "").strip())
            self._row = None
            if match:
                self._row = {
                    "species_number": match["number"],
                    "charge": match["charge"],
                    "cas": match["cas"],
                }
                self.rows.append(self._row)
            return
        if self._row is None:
            return
        if tag == "img":
            self._row["smiles"] = attributes.get("alt", "").strip()
        elif tag == "button":
            self._in_button = True
        elif tag == "span":
            css = attributes.get("class", "")
            if css in {"DHf298", "Uncert", "Units", "ATcTID"}:
                self._field = css
        elif tag == "a" and "species_number=" in attributes.get("href", ""):
            self._field = "name"

    def handle_endtag(self, tag: str) -> None:
        if tag in {"span", "a"}:
            self._field = None
        elif tag == "button":
            self._in_button = False
        elif tag == "tr":
            self._row = None

    def handle_data(self, data: str) -> None:
        if self._row is None:
            return
        if self._in_button:
            self._row["formula"] = self._row.get("formula", "") + data
        elif self._field is not None:
            self._row[self._field] = self._row.get(self._field, "") + data


def atct_phase(formula_text: str, atct_id: str) -> str:
    """``g``, ``g_variant`` (non-reference gas form) or ``condensed``."""
    match = ATCT_PHASE_LABEL.search(formula_text.strip())
    label = match.group(1).strip() if match else ""
    if not label.startswith("g"):
        return "condensed"
    return "g" if atct_id.endswith(GAS_REFERENCE_SUFFIX) else "g_variant"


def atct_version(html: str) -> str:
    match = re.search(r"based on version\s+([\d.]+\w*)", html)
    return match.group(1) if match else "unknown"


def read_atct(path: Path) -> tuple[list[SourceValue], dict[str, Any]]:
    html = path.read_text(encoding="utf-8", errors="replace")
    parser = _AtctTableParser()
    parser.feed(html)
    parser.close()

    values: list[SourceValue] = []
    seen: set[str] = set()
    skipped_units = 0
    for row in parser.rows:
        number = row["species_number"]
        if number in seen:  # tolerate repeated rows in a saved page
            continue
        seen.add(number)
        atct_id = row.get("ATcTID", "").strip()
        units = row.get("Units", "").strip()
        h298_text = row.get("DHf298", "").strip()
        if not h298_text:
            continue
        if units and units != "kJ/mol":
            skipped_units += 1
            continue
        values.append(
            SourceValue(
                source="atct",
                source_id=atct_id or f"species_{number}",
                name=row.get("name", "").strip(),
                raw_smiles=row.get("smiles", ""),
                h298_kj=float(h298_text),
                uncertainty_kj=parse_uncertainty(row.get("Uncert", "")),
                phase=atct_phase(row.get("formula", ""), atct_id),
                declared_charge=int(row["charge"]),
            )
        )
    info = {
        "version": atct_version(html),
        "species_rows": len(seen),
        "non_kj_rows_skipped": skipped_units,
    }
    return values, info


# ---------------------------------------------------------------------- RMG


def read_rmg(directory: Path) -> tuple[list[SourceValue], dict[str, Any]]:
    import yaml  # type: ignore[import-untyped]

    values: list[SourceValue] = []
    files = sorted(directory.glob("*.yml"))
    for path in files:
        with path.open(encoding="utf-8") as handle:
            entry = yaml.safe_load(handle)
        for tag, data in (entry.get("reference_data") or {}).items():
            h298 = (data.get("thermo_data") or {}).get("H298")
            if not h298:
                continue
            if h298.get("units") != "kJ/mol":
                raise ValueError(f"{path}: unexpected H298 units {h298.get('units')}")
            uncertainty = h298.get("uncertainty")
            values.append(
                SourceValue(
                    source=f"rmg:{tag}",
                    source_id=path.stem,
                    name=path.stem,
                    raw_smiles=entry.get("smiles", ""),
                    h298_kj=float(h298["value"]),
                    uncertainty_kj=(
                        float(uncertainty) if uncertainty not in (None, "") else None
                    ),
                    phase="g",
                )
            )
    return values, {"files": len(files)}


# -------------------------------------------------------------------- Bains


def _bains_rows(path: Path) -> Iterator[tuple[Any, ...]]:
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        yield from workbook[BAINS_SHEET].iter_rows(values_only=True)
    finally:
        workbook.close()


def read_bains(path: Path) -> tuple[list[SourceValue], dict[str, Any]]:
    rows = _bains_rows(path)
    header: tuple[Any, ...] = ()
    for index, row in enumerate(rows):
        if index == BAINS_HEADER_ROW:
            header = row
            break
    names = [str(cell).strip() if cell is not None else "" for cell in header]
    if names[:2] != ["Name", "SMILES"] or BAINS_FILTER_COLUMN not in names:
        raise ValueError(f"{path}: unexpected header {names[:4]}")
    filter_index = names.index(BAINS_FILTER_COLUMN)
    source_columns = list(range(2, filter_index))

    values: list[SourceValue] = []
    data_rows = 0
    filtered_rows = 0
    for row in rows:
        smiles = str(row[1]).strip() if row[1] is not None else ""
        if not smiles:
            continue
        data_rows += 1
        excluded = row[filter_index] not in (None, "", 0)
        filtered_rows += int(excluded)
        for column in source_columns:
            cell = row[column]
            if cell is None or (isinstance(cell, str) and not cell.strip()):
                continue
            values.append(
                SourceValue(
                    source=f"bains:{names[column]}",
                    source_id=f"row{data_rows}",
                    name=str(row[0] or "").strip(),
                    raw_smiles=smiles,
                    h298_kj=float(cell),
                    uncertainty_kj=None,
                    phase="g",
                    selectable=not excluded,
                )
            )
    return values, {"rows": data_rows, "filter_flagged_rows": filtered_rows}


# ---------------------------------------------------------------- selection


def long_table(values: Iterable[SourceValue]) -> tuple[pd.DataFrame, dict[str, int]]:
    """Screen every source value; return accepted rows and reject counts."""
    cache: dict[str, tuple[Structure | None, str]] = {}
    records: list[dict[str, Any]] = []
    rejected: dict[str, int] = {}
    for value in values:
        if value.phase != "g":
            reason = "gas_variant" if value.phase == "g_variant" else "phase"
            rejected[reason] = rejected.get(reason, 0) + 1
            continue
        if value.declared_charge not in (None, 0):
            rejected["charged"] = rejected.get("charged", 0) + 1
            continue
        if value.raw_smiles not in cache:
            cache[value.raw_smiles] = standardize(value.raw_smiles)
        structure, reason = cache[value.raw_smiles]
        if structure is None:
            rejected[reason] = rejected.get(reason, 0) + 1
            continue
        uncertainty = value.uncertainty_kj
        records.append(
            {
                **asdict(structure),
                "source": value.source,
                "source_rank": SOURCE_RANK.get(value.source, UNRANKED),
                "source_id": value.source_id,
                "name": value.name,
                "raw_smiles": value.raw_smiles,
                "H298_exp": kj_to_kcal(value.h298_kj),
                "uncertainty": (
                    kj_to_kcal(uncertainty) if uncertainty is not None else math.nan
                ),
                "selectable": value.selectable,
            }
        )
    return pd.DataFrame.from_records(records), rejected


def _conflicts(chosen: pd.Series, others: pd.DataFrame) -> int:
    chosen_u = 0.0 if pd.isna(chosen["uncertainty"]) else chosen["uncertainty"]
    count = 0
    for _, other in others.iterrows():
        other_u = 0.0 if pd.isna(other["uncertainty"]) else other["uncertainty"]
        limit = max(2.0 * math.hypot(chosen_u, other_u), CONFLICT_FLOOR_KCAL)
        if abs(other["H298_exp"] - chosen["H298_exp"]) > limit:
            count += 1
    return count


def _collapse_isomers(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep the lowest H298 per source; return (kept, dropped).

    Distinct species of one source can share an InChIKey when the source
    SMILES omits stereo or the standard InChI merges tautomers (ATcT lists
    cycloheptene and trans-cycloheptene, urea and isourea, E/Z imines with
    the same SMILES). A key without that distinction stands for the most
    stable form, so the lowest enthalpy is kept.
    """
    ordered = candidates.sort_values(["source", "H298_exp", "source_id"])
    first_id = ordered.groupby("source")["source_id"].transform("first")
    same_species = ordered["source_id"] == first_id
    return ordered[same_species], ordered[~same_species]


def select(long: pd.DataFrame) -> pd.DataFrame:
    """One row per InChIKey: best rank, then lowest uncertainty, then source id."""
    rows: list[dict[str, Any]] = []
    for inchikey, group in long.groupby("inchikey", sort=True):
        candidates, collapsed = _collapse_isomers(group[group["selectable"]])
        if candidates.empty:
            continue
        ordered = candidates.assign(
            _unc=candidates["uncertainty"].fillna(math.inf)
        ).sort_values(["source_rank", "_unc", "source", "source_id"])
        chosen = ordered.iloc[0]
        compared = group.drop(index=collapsed.index)
        others = compared.drop(index=chosen.name)
        rows.append(
            {
                "inchikey": inchikey,
                "smiles": chosen["smiles"],
                "formula": chosen["formula"],
                "charge": int(chosen["charge"]),
                "multiplicity": int(chosen["multiplicity"]),
                "nheavy": int(chosen["nheavy"]),
                "phase": "g",
                "H298_exp": chosen["H298_exp"],
                "uncertainty": chosen["uncertainty"],
                "uncertainty_kind": (
                    "reported" if pd.notna(chosen["uncertainty"]) else "none"
                ),
                "source": chosen["source"],
                "source_id": chosen["source_id"],
                "name": chosen["name"],
                "n_values": len(compared),
                "n_sources": compared["source"].nunique(),
                "source_spread": compared["H298_exp"].max()
                - compared["H298_exp"].min(),
                "n_conflicts": _conflicts(chosen, others),
                "isomers_collapsed": collapsed["source_id"].nunique(),
                "low_priority_only": bool(
                    (candidates["source_rank"] >= LOW_PRIORITY_RANK).all()
                ),
                "stewart_reference": bool((group["source"] == "bains:Stewart").any()),
            }
        )
    table = pd.DataFrame.from_records(rows)
    table["conflict_flag"] = table["n_conflicts"] > 0
    table.insert(0, "exp_id", [f"exp_{i:05d}" for i in range(len(table))])
    return table


# ----------------------------------------------------------- CBS comparison


def compare_with_cbs(table: pd.DataFrame, cbs_csv: Path) -> dict[str, Any]:
    """Overlap and CBS - experimental statistics on shared InChIKeys.

    The CBS file can hold several reference conformers per molecule; the
    lowest H298_cbs per InChIKey is used (most stable conformer).
    """
    cbs = pd.read_csv(
        cbs_csv, usecols=["mol_id", "smiles", "multiplicity", "charge", "H298_cbs"]
    )
    cbs = cbs[(cbs["multiplicity"] == 1) & (cbs["charge"] == 0)].copy()
    keys = {s: inchikey_of(Chem.MolFromSmiles(s)) for s in cbs["smiles"].unique()}
    cbs["inchikey"] = cbs["smiles"].map(keys)
    cbs = cbs[cbs["inchikey"].astype(bool)]
    best = cbs.sort_values("H298_cbs").drop_duplicates("inchikey")

    merged = table.merge(best[["inchikey", "mol_id", "H298_cbs"]], on="inchikey")
    skeleton = set(best["inchikey"].str[:14])
    skeleton_only = int(
        (~table["inchikey"].isin(best["inchikey"]))
        .where(table["inchikey"].str[:14].isin(skeleton), False)
        .sum()
    )
    report: dict[str, Any] = {
        "cbs_csv": str(cbs_csv),
        "cbs_sha256": sha256_of(cbs_csv),
        "overlap_inchikey": len(merged),
        "overlap_skeleton_only": skeleton_only,
    }
    if merged.empty:
        return report
    error = merged["H298_cbs"] - merged["H298_exp"]
    merged = merged.assign(cbs_minus_exp=error, abs_error=error.abs())

    def stats(frame: pd.DataFrame) -> dict[str, float | int]:
        err = frame["cbs_minus_exp"]
        return {
            "n": len(frame),
            "mae": float(err.abs().mean()),
            "bias": float(err.mean()),
            "rmse": float(math.sqrt((err**2).mean())),
            "max_abs": float(err.abs().max()),
        }

    report["units"] = "kcal/mol, CBS - experimental"
    report["all"] = stats(merged)
    report["by_source_family"] = {
        family: stats(group)
        for family, group in merged.groupby(merged["source"].str.split(":").str[0])
    }
    report["excluding_conflicts_and_low_priority"] = stats(
        merged[~merged["conflict_flag"] & ~merged["low_priority_only"]]
    )
    report["outliers_top20"] = (
        merged.sort_values("abs_error", ascending=False)
        .head(20)[
            [
                "exp_id",
                "mol_id",
                "smiles",
                "source",
                "H298_exp",
                "H298_cbs",
                "cbs_minus_exp",
                "conflict_flag",
            ]
        ]
        .round(3)
        .to_dict("records")
    )
    return report


# --------------------------------------------------------------------- main


def build(
    atct_html: Path | None, rmg_dir: Path | None, bains_xlsx: Path | None
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    values: list[SourceValue] = []
    info: dict[str, Any] = {}
    readers = (
        ("atct", atct_html, read_atct),
        ("rmg", rmg_dir, read_rmg),
        ("bains", bains_xlsx, read_bains),
    )
    for name, path, reader in readers:
        if path is None:
            continue
        found, details = reader(path)
        values.extend(found)
        details["path"] = str(path)
        details["values"] = len(found)
        if path.is_file():
            details["sha256"] = sha256_of(path)
        info[name] = details
    if not values:
        raise ValueError("no source given")
    long, rejected = long_table(values)
    if long.empty:
        raise ValueError("no source value passed the CHON closed-shell screen")
    table = select(long)
    info["rejected_values"] = rejected
    return table, long, info


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--atct-html", type=Path, default=None)
    parser.add_argument("--atct-version", default=None, help="override parsed version")
    parser.add_argument("--rmg-dir", type=Path, default=None)
    parser.add_argument("--rmg-commit", default="unknown")
    parser.add_argument("--bains-xlsx", type=Path, default=None)
    parser.add_argument("--cbs-csv", type=Path, default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    all_sources = args.output.with_name(args.output.stem + "_all_sources.csv")
    manifest_path = args.output.with_suffix(".manifest.json")
    for target in (args.output, all_sources, manifest_path):
        if target.exists():
            print(f"refusing to overwrite {target}", file=sys.stderr)
            return 1

    RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]
    table, long, info = build(args.atct_html, args.rmg_dir, args.bains_xlsx)
    if "atct" in info and args.atct_version:
        info["atct"]["version"] = args.atct_version

    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False, float_format="%.4f")
    long.drop(columns=["charge", "multiplicity"]).to_csv(
        all_sources, index=False, float_format="%.4f"
    )

    citations = []
    if "atct" in info:
        citations.append(CITATIONS["atct"].format(version=info["atct"]["version"]))
    if "rmg" in info:
        citations.append(CITATIONS["rmg"].format(commit=args.rmg_commit))
    if "bains" in info:
        citations.append(CITATIONS["bains"])

    manifest: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "sources": info,
        "citations": citations,
        "license": (
            "CC-BY-SA (share-alike inherited from Bains et al. 2022)"
            if "bains" in info
            else "see citations"
        ),
        "license_note": "ATcT states no explicit licence; its values are "
        "redistributed with citation. RMG-database carries an MIT header.",
        "filter": "contains C; elements subset of {C,H,N,O}; no isotopes; "
        "formal charge 0; no radical electrons; gas phase",
        "selection": "per InChIKey: within one source keep the lowest H298 "
        "species (stereo/tautomer collapse), then lowest source_rank, then "
        "lowest reported uncertainty; Bains FILTER==1 never selected",
        "source_rank": SOURCE_RANK,
        "conflict_rule": "|dH| > max(2*hypot(u1,u2), "
        f"{CONFLICT_FLOOR_KCAL} kcal/mol), unknown u = 0",
        "units": {
            "H298_exp": "kcal/mol",
            "uncertainty": "kcal/mol",
            "source_spread": "kcal/mol",
        },
        "kj_per_kcal": KJ_PER_KCAL,
        "rows": len(table),
        "long_rows": len(long),
        "selected_by_source": table["source"].value_counts().to_dict(),
        "conflict_rows": int(table["conflict_flag"].sum()),
        "low_priority_only_rows": int(table["low_priority_only"].sum()),
        "stewart_reference_rows": int(table["stewart_reference"].sum()),
        "isomers_collapsed_rows": int((table["isomers_collapsed"] > 0).sum()),
        "skeleton_duplicates": int(table["inchikey"].str[:14].duplicated().sum()),
        "nheavy": table["nheavy"].value_counts().sort_index().to_dict(),
        "output_sha256": sha256_of(args.output),
        "all_sources_sha256": sha256_of(all_sources),
    }
    if args.cbs_csv is not None:
        manifest["cbs_comparison"] = compare_with_cbs(table, args.cbs_csv)

    manifest_path.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    summary = {k: manifest[k] for k in ("rows", "long_rows", "conflict_rows")}
    if "cbs_comparison" in manifest:
        summary["cbs_overlap"] = manifest["cbs_comparison"]["overlap_inchikey"]
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
