"""Reaction enthalpy calculations."""

from __future__ import annotations

from pathlib import Path

from .io import read_csv_rows, read_xlsx_table, write_csv
from .metrics import absolute_deviation, relative_deviation_percent


def reaction_enthalpy(
    reactants: dict[str, float],
    products: dict[str, float],
) -> float:
    """Compute sum(nu products Hf) - sum(nu reactants Hf)."""
    return sum(products.values()) - sum(reactants.values())


def load_reaction_rows(path: Path) -> list[dict[str, str]]:
    if path.suffix.lower() == ".csv":
        return read_csv_rows(path)
    return read_xlsx_table(path, "04_reacoes", {"reaction_id", "coef_a", "coef_d"})


def calculate_reactions(
    reaction_rows: list[dict[str, str]],
    formation_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Calculate AM1/PM3/PM7 reaction enthalpies from formation rows."""
    hfs: dict[str, dict[str, float]] = {}
    for formation_row in formation_rows:
        mol_id = str(formation_row.get("mol_id", ""))
        method = str(formation_row.get("method", ""))
        value = formation_row.get("hof_selected_kJmol")
        if mol_id and method and value not in {None, ""}:
            hfs.setdefault(mol_id, {})[method] = float(value)  # type: ignore[arg-type]

    out = []
    for row in reaction_rows:
        result: dict[str, object] = {
            "reaction_id": _text(row.get("reaction_id")),
            "reaction_type": _text(row.get("reaction_type")),
            "A": _text(row.get("compound_A")) or _text(row.get("A")),
            "B": _text(row.get("compound_B")) or _text(row.get("B")),
            "C": _text(row.get("compound_C")) or _text(row.get("C")),
            "D": _text(row.get("compound_D")) or _text(row.get("D")),
            "a": _coef(_text(row.get("coef_a")) or _text(row.get("a"))),
            "b": _coef(_text(row.get("coef_b")) or _text(row.get("b"))),
            "c": _coef(_text(row.get("coef_c")) or _text(row.get("c"))),
            "d": _coef(_text(row.get("coef_d")) or _text(row.get("d"))),
            "Hr_exp_kJmol": _float(_text(row.get("Hr_exp_kJmol"))),
        }
        for method in ("AM1", "PM3", "PM7"):
            hr = _hr_for_method(result, hfs, method)
            result[f"Hr_{method}_crest_mopac_kJmol"] = hr
            result[f"AD_{method}_kJmol"] = absolute_deviation(
                hr, result["Hr_exp_kJmol"]  # type: ignore[arg-type]
            )
            result[f"RD_{method}_percent"] = relative_deviation_percent(
                hr, result["Hr_exp_kJmol"]  # type: ignore[arg-type]
            )
        out.append(result)
    return out


def write_reactions_csv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = [
        "reaction_id",
        "reaction_type",
        "A",
        "B",
        "C",
        "D",
        "a",
        "b",
        "c",
        "d",
        "Hr_exp_kJmol",
        "Hr_AM1_crest_mopac_kJmol",
        "Hr_PM3_crest_mopac_kJmol",
        "Hr_PM7_crest_mopac_kJmol",
        "AD_AM1_kJmol",
        "AD_PM3_kJmol",
        "AD_PM7_kJmol",
        "RD_AM1_percent",
        "RD_PM3_percent",
        "RD_PM7_percent",
    ]
    write_csv(path, rows, columns)


def _hr_for_method(
    row: dict[str, object], hfs: dict[str, dict[str, float]], method: str
) -> float | None:
    try:
        reactants = {}
        products = {}
        if row["A"]:
            reactants[str(row["A"])] = (
                _required_float(row["a"]) * hfs[str(row["A"])][method]
            )
        if row["B"]:
            reactants[str(row["B"])] = (
                _required_float(row["b"]) * hfs[str(row["B"])][method]
            )
        if row["C"]:
            products[str(row["C"])] = (
                _required_float(row["c"]) * hfs[str(row["C"])][method]
            )
        if row["D"]:
            products[str(row["D"])] = (
                _required_float(row["d"]) * hfs[str(row["D"])][method]
            )
        return reaction_enthalpy(reactants, products)
    except KeyError:
        return None


def _coef(value: str | None) -> float:
    number = _float(value)
    return 0.0 if number is None else number


def _float(value: str | None) -> float | None:
    if value in {None, ""}:
        return None
    try:
        return float(str(value).replace(",", "."))
    except ValueError:
        return None


def _required_float(value: object) -> float:
    return float(str(value))


def _text(value: object | None) -> str:
    return "" if value is None else str(value)
