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
4. NIST Chemistry WebBook (SRD 69), compound pages cached by
   ``fetch_nist_webbook.py``. Identity comes from the page InChI; the gas
   value is the WebBook ``AVG`` row when listed. NIST values are
   validation-only, so they rank below every trainable measured value
   (after Bains ``Winget``): a NIST value is selected only when no
   trainable measured source covers the molecule. Gas values from ion
   energetics (method ``Ion``) are kept as ``nist:gas_ion`` with the lowest
   priority.

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
  counts as zero;
* NIST values are validation-only (the site states "All rights reserved"
  and its robots.txt sets ``ai-train=no``): every NIST value, and every
  selected row taken from NIST, carries ``validation_only=True``;
* a NIST InChI must reproduce the page's own formula and InChIKey, so a
  page whose structure and metadata disagree is rejected, not guessed;
* NIST values without a reported uncertainty stay in the long table but are
  never selected: they are mostly single old measurements, and the gross
  errors seen against ATcT and CBS (ethanolamine, hexylamine) are all in
  this group, where no second source can raise a conflict;
* per NIST quantity one value is used: the WebBook ``AVG`` row, else the
  lowest reported uncertainty, else the lower median; ``Ion`` gas rows only
  when the page has no other gas value;
* ``ΔfH°(gas) = ΔfH°(liquid) + ΔvapH°`` (or solid + ``ΔsubH°``) only with
  transition enthalpies at 298.15 K: the ``°`` rows (standard conditions)
  or table rows within 1 K of 298.15 K; other temperatures are ignored,
  not extrapolated. The uncertainty is the quadrature sum, or unknown if
  either part has none.

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
    "nist:gas": 20,
    "nist:liquid+vap": 21,
    "nist:solid+sub": 21,
    "bains:Yaws": 90,
    "bains:Stewart": 91,
    "nist:gas_ion": 92,
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
    "nist": "Linstrom, P. J.; Mallard, W. G. (eds.) NIST Chemistry WebBook, "
    "NIST Standard Reference Database Number 69, National Institute of "
    "Standards and Technology, Gaithersburg MD, doi:10.18434/T4D303 "
    "(retrieved {retrieved})",
}
NIST_TERMS = (
    "NIST SRD 69 data: (c) U.S. Secretary of Commerce, all rights reserved; "
    "webbook.nist.gov robots.txt sets Content-Signal ai-train=no. Every NIST "
    "value carries validation_only=True and must not be used as an ML "
    "training target; NIST rows are not covered by the CC-BY-SA licence."
)


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
    #: NIST values: never to be used as ML training targets (site terms).
    validation_only: bool = False


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


# --------------------------------------------------------------------- NIST


@dataclass(frozen=True)
class NistPick:
    """The value chosen for one quantity of one NIST compound page."""

    value_kj: float
    uncertainty_kj: float | None
    method: str
    n_rows: int


@dataclass(frozen=True)
class NistPage:
    compound_id: str
    name: str
    formula: str
    inchi: str
    inchikey: str
    picks: dict[str, NistPick]
    off_298_rows: int


class _NistPageParser(HTMLParser):
    """Collect identity fields and every data table of a WebBook page."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.name = ""
        self.inchi_texts: list[str] = []
        self.li_texts: list[str] = []
        self.tables: list[dict[str, Any]] = []
        self._section = ""
        self._in_title = False
        self._in_inchi = False
        self._li_stack: list[str] = []
        self._cell: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = {key: value or "" for key, value in attrs}
        if tag == "h1" and attributes.get("id") == "Top":
            self._in_title = True
        elif tag == "h2":
            self._section = attributes.get("id", "")
        elif tag == "span" and attributes.get("class") == "inchi-text":
            self._in_inchi = True
            self.inchi_texts.append("")
        elif tag == "li":
            self._li_stack.append("")
        elif tag == "table":
            self.tables.append(
                {
                    "section": self._section,
                    "label": attributes.get("aria-label", ""),
                    "rows": [],
                }
            )
        elif tag == "tr" and self.tables:
            self.tables[-1]["rows"].append([])
        elif tag in {"td", "th"} and self.tables and self.tables[-1]["rows"]:
            self._cell = []
            self.tables[-1]["rows"][-1].append(self._cell)

    def handle_endtag(self, tag: str) -> None:
        if tag == "h1":
            self._in_title = False
        elif tag == "span":
            self._in_inchi = False
        elif tag == "li" and self._li_stack:
            text = self._li_stack.pop()
            self.li_texts.append(text)
            if self._li_stack:
                self._li_stack[-1] += text
        elif tag in {"td", "th"}:
            self._cell = None

    def handle_data(self, data: str) -> None:
        if self._in_title:
            self.name += data
        if self._in_inchi:
            self.inchi_texts[-1] += data
        if self._li_stack:
            self._li_stack[-1] += data
        if self._cell is not None:
            self._cell.append(data)


#: Quantity label (whitespace removed) of the one-dimensional data tables.
NIST_QUANTITIES = {
    "ΔfH°gas": "gas",
    "ΔfH°liquid": "liquid",
    "ΔfH°solid": "solid",
    "ΔvapH°": "vap",
    "ΔsubH°": "sub",
}
#: Header of the temperature-dependent tables, used only at 298.15 K.
NIST_T_TABLES = {"ΔvapH(kJ/mol)": "vap", "ΔsubH(kJ/mol)": "sub"}
NIST_T_TOLERANCE_K = 1.0
NIST_ION_METHOD = "Ion"
NIST_DERIVED = {
    "nist:liquid+vap": ("liquid", "vap"),
    "nist:solid+sub": ("solid", "sub"),
}


def parse_nist_value(text: str) -> tuple[float, float | None]:
    """``-234. ± 2.`` -> (-234.0, 2.0); ``-323.6`` -> (-323.6, None)."""
    value, _, uncertainty = text.partition("±")
    return float(value), float(uncertainty) if uncertainty.strip() else None


def pick_nist(rows: list[tuple[float, float | None, str]]) -> NistPick:
    """NIST average if listed, else lowest reported uncertainty, else median.

    The WebBook ``AVG`` row is its own average of the individual points.
    Without it, the most precise value wins; with no uncertainty at all, the
    lower median keeps one actual reported value.
    """
    averages = [row for row in rows if row[2] == "AVG"]
    if averages:
        value, uncertainty, method = averages[0]
    else:
        with_u = [row for row in rows if row[1] is not None]
        if with_u:
            value, uncertainty, method = min(with_u, key=lambda row: row[1] or 0.0)
        else:
            ordered = sorted(rows)
            value, uncertainty, method = ordered[(len(ordered) - 1) // 2]
    return NistPick(value, uncertainty, method, len(rows))


def _cell_text(cell: list[str]) -> str:
    return " ".join("".join(cell).split())


def parse_nist_page(page: str, compound_id: str) -> NistPage:
    parser = _NistPageParser()
    parser.feed(page)
    parser.close()

    standard: dict[str, list[tuple[float, float | None, str]]] = {}
    at_298: dict[str, list[tuple[float, float | None, str]]] = {}
    off_298 = 0
    for table in parser.tables:
        rows = [[_cell_text(cell) for cell in row] for row in table["rows"]]
        header = rows[0][0].replace(" ", "") if rows and rows[0] else ""
        if header in NIST_T_TABLES:
            quantity = NIST_T_TABLES[header]
            for cells in rows[1:]:
                try:
                    value, temperature = float(cells[0]), float(cells[1])
                except (IndexError, ValueError):
                    continue
                if abs(temperature - 298.15) <= NIST_T_TOLERANCE_K:
                    at_298.setdefault(quantity, []).append((value, None, "T298"))
                else:
                    off_298 += 1
            continue
        for cells in rows:
            if len(cells) < 4:
                continue
            kind = NIST_QUANTITIES.get(cells[0].replace(" ", ""))
            if kind is None or cells[2] != "kJ/mol":
                continue
            try:
                value, uncertainty = parse_nist_value(cells[1])
            except ValueError:
                continue
            standard.setdefault(kind, []).append((value, uncertainty, cells[3]))

    # Ion-energetics gas values are indirect (appearance energies, typical
    # uncertainty 8-12 kJ/mol): used only when no other gas value exists.
    gas_rows = standard.pop("gas", [])
    direct = [row for row in gas_rows if row[2] != NIST_ION_METHOD]
    if direct:
        standard["gas"] = direct
    elif gas_rows:
        standard["gas_ion"] = gas_rows
    picks = {kind: pick_nist(found) for kind, found in standard.items()}
    for kind, found in at_298.items():
        picks.setdefault(kind, pick_nist(found))

    inchi = next((t.strip() for t in parser.inchi_texts if t.startswith("InChI=")), "")
    inchikey = next(
        (t.strip() for t in parser.inchi_texts if not t.startswith("InChI=")), ""
    )
    formula = next(
        (
            "".join(text.split(":", 1)[1].split())
            for text in parser.li_texts
            if text.strip().startswith("Formula:")
        ),
        "",
    )
    return NistPage(
        compound_id=compound_id,
        name=" ".join(parser.name.split()),
        formula=formula,
        inchi=inchi,
        inchikey=inchikey,
        picks=picks,
        off_298_rows=off_298,
    )


def nist_smiles(page: NistPage) -> tuple[str | None, str]:
    """SMILES of the page InChI after checking it against the page's own data.

    The InChI defines the identity; its formula must equal the listed
    formula and the InChIKey of the SMILES round trip must equal the listed
    InChIKey, so a wrong structure cannot pass silently.
    """
    if not page.inchi:
        return None, "no_inchi"
    mol = Chem.MolFromInchi(page.inchi)  # type: ignore[no-untyped-call]
    if mol is None:
        return None, "unparsable_inchi"
    if CalcMolFormula(mol) != page.formula:
        return None, "formula_mismatch"
    smiles = Chem.MolToSmiles(mol)
    roundtrip = Chem.MolFromSmiles(smiles)
    expected = page.inchikey or str(
        Chem.InchiToInchiKey(page.inchi)  # type: ignore[no-untyped-call]
    )
    if roundtrip is None or inchikey_of(roundtrip) != expected:
        return None, "inchikey_mismatch"
    return smiles, ""


def nist_values(page: NistPage, smiles: str) -> list[SourceValue]:
    """Measured gas value and condensed + vaporization/sublimation routes."""

    def value(source: str, h298_kj: float, uncertainty_kj: float | None) -> SourceValue:
        return SourceValue(
            source=source,
            source_id=page.compound_id,
            name=page.name,
            raw_smiles=smiles,
            h298_kj=h298_kj,
            uncertainty_kj=uncertainty_kj,
            phase="g",
            selectable=uncertainty_kj is not None,
            validation_only=True,
        )

    values: list[SourceValue] = []
    for kind in ("gas", "gas_ion"):
        gas = page.picks.get(kind)
        if gas is not None:
            values.append(value(f"nist:{kind}", gas.value_kj, gas.uncertainty_kj))
    for source, (condensed, transition) in NIST_DERIVED.items():
        first, second = page.picks.get(condensed), page.picks.get(transition)
        if first is None or second is None:
            continue
        uncertainty = (
            math.hypot(first.uncertainty_kj, second.uncertainty_kj)
            if first.uncertainty_kj is not None and second.uncertainty_kj is not None
            else None
        )
        values.append(value(source, first.value_kj + second.value_kj, uncertainty))
    return values


def _fetch_dates(directory: Path) -> list[str]:
    log = directory / "fetch_log.jsonl"
    if not log.is_file():
        return []
    with log.open(encoding="utf-8") as handle:
        times = sorted(str(json.loads(line).get("time", ""))[:10] for line in handle)
    times = [t for t in times if t]
    return [times[0], times[-1]] if times else []


def read_nist(directory: Path) -> tuple[list[SourceValue], dict[str, Any]]:
    """Compound pages cached by ``fetch_nist_webbook.py`` (gzip HTML)."""
    import gzip

    files = sorted((directory / "compound").glob("*.html.gz"))
    values: list[SourceValue] = []
    rejected: dict[str, int] = {}
    off_298 = 0
    for path in files:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            page = parse_nist_page(handle.read(), path.name.removesuffix(".html.gz"))
        off_298 += page.off_298_rows
        smiles, reason = nist_smiles(page)
        if smiles is None:
            rejected[reason] = rejected.get(reason, 0) + 1
            continue
        found = nist_values(page, smiles)
        if not found:
            rejected["no_usable_value"] = rejected.get("no_usable_value", 0) + 1
        values.extend(found)
    info = {
        "compound_pages": len(files),
        "identity_rejects": rejected,
        "transition_rows_off_298_ignored": off_298,
        "retrieved": _fetch_dates(directory),
        "values_by_route": {
            source: sum(v.source == source for v in values)
            for source in ("nist:gas", "nist:gas_ion", *NIST_DERIVED)
        },
    }
    return values, info


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
                "validation_only": value.validation_only,
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
                "validation_only": bool(chosen["validation_only"]),
            }
        )
    table = pd.DataFrame.from_records(rows)
    table["conflict_flag"] = table["n_conflicts"] > 0
    table.insert(0, "exp_id", [f"exp_{i:05d}" for i in range(len(table))])
    return table


# ----------------------------------------------------------- CBS comparison


COMPARISON_COLUMNS = [
    "exp_id",
    "mol_id",
    "inchikey",
    "smiles",
    "nheavy",
    "source",
    "H298_exp",
    "uncertainty",
    "H298_cbs",
    "cbs_minus_exp",
    "conflict_flag",
    "low_priority_only",
    "validation_only",
]


def compare_with_cbs(
    table: pd.DataFrame, cbs_csv: Path
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Overlap and CBS - experimental statistics on shared InChIKeys.

    The CBS file can hold several reference conformers per molecule; the
    lowest H298_cbs per InChIKey is used (most stable conformer). Also
    returns the per-molecule comparison (``COMPARISON_COLUMNS``, sorted by
    ``exp_id``) so plots can use every point, not only the top outliers.
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
    error = merged["H298_cbs"] - merged["H298_exp"]
    merged = merged.assign(cbs_minus_exp=error, abs_error=error.abs())
    comparison = merged.sort_values("exp_id")[COMPARISON_COLUMNS].reset_index(drop=True)
    if merged.empty:
        return report, comparison

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
    return report, comparison


# --------------------------------------------------------------------- main


def build(
    atct_html: Path | None,
    rmg_dir: Path | None,
    bains_xlsx: Path | None,
    nist_dir: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    values: list[SourceValue] = []
    info: dict[str, Any] = {}
    readers = (
        ("atct", atct_html, read_atct),
        ("rmg", rmg_dir, read_rmg),
        ("bains", bains_xlsx, read_bains),
        ("nist", nist_dir, read_nist),
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
    parser.add_argument("--nist-dir", type=Path, default=None, help="fetch cache")
    parser.add_argument("--cbs-csv", type=Path, default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    all_sources = args.output.with_name(args.output.stem + "_all_sources.csv")
    manifest_path = args.output.with_suffix(".manifest.json")
    vs_cbs = args.output.with_name(args.output.stem + "_vs_cbs.csv")
    targets = [args.output, all_sources, manifest_path]
    if args.cbs_csv is not None:
        targets.append(vs_cbs)
    for target in targets:
        if target.exists():
            print(f"refusing to overwrite {target}", file=sys.stderr)
            return 1

    RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]
    table, long, info = build(
        args.atct_html, args.rmg_dir, args.bains_xlsx, args.nist_dir
    )
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
    if "nist" in info:
        retrieved = " to ".join(info["nist"]["retrieved"]) or "unknown date"
        citations.append(CITATIONS["nist"].format(retrieved=retrieved))

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
        "redistributed with citation. RMG-database carries an MIT header."
        + (" " + NIST_TERMS if "nist" in info else ""),
        "filter": "contains C; elements subset of {C,H,N,O}; no isotopes; "
        "formal charge 0; no radical electrons; gas phase",
        "selection": "per InChIKey: within one source keep the lowest H298 "
        "species (stereo/tautomer collapse), then lowest source_rank, then "
        "lowest reported uncertainty; Bains FILTER==1 never selected; "
        "NIST values without reported uncertainty never selected; "
        "validation_only follows the selected value (True for NIST)",
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
        "validation_only_rows": int(table["validation_only"].sum()),
        "isomers_collapsed_rows": int((table["isomers_collapsed"] > 0).sum()),
        "skeleton_duplicates": int(table["inchikey"].str[:14].duplicated().sum()),
        "nheavy": table["nheavy"].value_counts().sort_index().to_dict(),
        "output_sha256": sha256_of(args.output),
        "all_sources_sha256": sha256_of(all_sources),
    }
    if args.cbs_csv is not None:
        report, comparison = compare_with_cbs(table, args.cbs_csv)
        comparison.to_csv(vs_cbs, index=False, float_format="%.4f")
        report["comparison_csv"] = str(vs_cbs)
        report["comparison_sha256"] = sha256_of(vs_cbs)
        manifest["cbs_comparison"] = report

    manifest_path.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    summary = {k: manifest[k] for k in ("rows", "long_rows", "conflict_rows")}
    if "cbs_comparison" in manifest:
        summary["cbs_overlap"] = manifest["cbs_comparison"]["overlap_inchikey"]
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
