"""Parsing, screening and selection of the experimental ΔHf builder."""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest


def _load_builder() -> ModuleType:
    repo_root = Path(__file__).resolve().parents[2]
    script_path = repo_root / "scripts" / "build_experimental_hf.py"
    spec = importlib.util.spec_from_file_location("build_experimental_hf", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolve annotations through sys.modules[cls.__module__].
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


builder = _load_builder()

# Synthetic values in the ATcT table layout (not ATcT data).
ATCT_ROW = (
    '<tr id="{rid} i{number} CAS{cas}"><td class="bkgName"><span class="Name"> '
    '<a href="species/?species_number={number}">{name}</a></span></td>'
    '<td class="bkgFormula&gt; &lt;span class=" formula"=""> '
    '<button onclick="copyname(\'{name}\');" type="button">{formula}  (g) </button>'
    ' </td><td class="bkgImage"><span class="Name"><img class="lazy" '
    'data-original="images/{number}.png" alt="{smiles}"></span></td>'
    '<td class="bkgDHf0"><span class="DHf0">0.0</span></td>'
    '<td class="bkgDHf298"><span class="DHf298">{h298}</span></td>'
    '<td class="bkgUncert"><span class="Uncert">{unc}</span></td>'
    '<td class="bkgUnits"><span class="Units">kJ/mol</span></td>'
    '<td class="bkgMass"><span class="Mass">1.0 ±<br>0.1</span></td>'
    '<td class="bkgATcTID"><span class="ATcTID">{cas}*{phase}</span></td></tr>'
)


def _atct_html(rows: list[dict[str, str]]) -> str:
    body = "".join(ATCT_ROW.format(**row) for row in rows)
    return (
        "<html><body><h1>enthalpies of formation based on version 9.999 of the "
        f"Thermochemical Network</h1><table><tbody>{body}</tbody></table>"
        "</body></html>"
    )


ATCT_ROWS = [
    {
        "rid": "s1_6n4_1c0",
        "number": "23",
        "cas": "74-82-8",
        "name": "Methane",
        "formula": "CH4",
        "smiles": "C",
        "h298": "-100.0",
        "unc": "± 0.40",
        "phase": "0",
    },
    {
        "rid": "s1_6n3_1c0",
        "number": "34",
        "cas": "2229-07-4",
        "name": "Methyl",
        "formula": "CH3",
        "smiles": "[CH3]",
        "h298": "150.0",
        "unc": "± 0.05",
        "phase": "0",
    },
    {
        "rid": "s1_6_8n4_1_1c0",
        "number": "40",
        "cas": "67-56-1",
        "name": "Methanol",
        "formula": "CH3OH",
        "smiles": "CO",
        "h298": "-240.0",
        "unc": "± 0.10",
        "phase": "500",
    },
]


def _write_atct(tmp_path: Path, rows: list[dict[str, str]] = ATCT_ROWS) -> Path:
    path = tmp_path / "atct.html"
    path.write_text(_atct_html(rows), encoding="utf-8")
    return path


def _write_rmg(tmp_path: Path) -> Path:
    directory = tmp_path / "rmg"
    directory.mkdir()
    (directory / "Methane.yml").write_text(
        "smiles: C\n"
        "reference_data:\n"
        "  CCCBDB:\n"
        "    class: ReferenceDataEntry\n"
        "    thermo_data:\n"
        "      H298:\n"
        "        units: kJ/mol\n"
        "        uncertainty: 0.5\n"
        "        value: -90.0\n"
    )
    (directory / "Ethanol.yml").write_text(
        "smiles: CCO\n"
        "reference_data:\n"
        "  ATcT:\n"
        "    thermo_data:\n"
        "      H298: {units: kJ/mol, uncertainty: 0.8, value: -230.0}\n"
        "  CCCBDB:\n"
        "    thermo_data:\n"
        "      H298: {units: kJ/mol, uncertainty: 0.2, value: -231.0}\n"
    )
    (directory / "NoThermo.yml").write_text(
        "smiles: CC\nreference_data:\n  CATCH:\n    atomization_energy: 1\n"
    )
    return directory


BAINS_COLUMNS = ["ATcT", "Pedley", "Yaws", "Stewart"]


def _write_bains(tmp_path: Path, rows: list[list[object]]) -> Path:
    openpyxl = pytest.importorskip("openpyxl")
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = builder.BAINS_SHEET
    sheet.append(["Molecules", None, "Enthalpy of formation (kJ/mol)"])
    sheet.append([None])
    sheet.append(["Name", "SMILES", *BAINS_COLUMNS, "FILTER", "PM7"])
    for row in rows:
        sheet.append(row)
    path = tmp_path / "bains.xlsx"
    workbook.save(path)
    return path


BAINS_ROWS: list[list[object]] = [
    # name, smiles, ATcT, Pedley, Yaws, Stewart, FILTER, PM7
    ["Ethanol", "CCO", "", -235.0, -234.0, "", None, -50.0],
    ["Propanol", "CCCO", "", "", -255.0, -256.0, None, -55.0],
    ["Acetone", "CC(C)=O", "", -217.0, "", "", 1, -50.0],
    ["Formic acid", "OC=O", "", -378.0, "", -390.0, None, -90.0],
    ["Sodium chloride", "[Na+].[Cl-]", "", -180.0, "", "", None, -40.0],
    ["Broken", "C1CC", "", -1.0, "", "", None, 0.0],
]


def test_kj_to_kcal_uses_shared_constant() -> None:
    assert builder.kj_to_kcal(41.84) == pytest.approx(41.84 / builder.KJ_PER_KCAL)


@pytest.mark.parametrize(
    ("text", "expected"),
    [("± 0.63", 0.63), ("exact", 0.0), ("", None), ("±0.043", 0.043)],
)
def test_parse_uncertainty(text: str, expected: float | None) -> None:
    assert builder.parse_uncertainty(text) == expected


@pytest.mark.parametrize(
    ("smiles", "reason"),
    [
        ("[CH3]", "open_shell"),
        ("C[NH3+]", "charged"),
        ("[2H]C", "isotope"),
        ("CCl", "elements"),
        ("O", "elements"),
        ("C1CC", "unparsable_smiles"),
        ("", "unparsable_smiles"),
    ],
)
def test_standardize_rejects(smiles: str, reason: str) -> None:
    structure, found = builder.standardize(smiles)
    assert structure is None
    assert found == reason


def test_standardize_accepts_closed_shell_chon() -> None:
    structure, reason = builder.standardize("OCC")
    assert reason == ""
    assert structure.smiles == "CCO"
    assert structure.inchikey == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert (structure.charge, structure.multiplicity, structure.nheavy) == (0, 1, 3)
    assert structure.formula == "C2H6O"


def test_read_atct_parses_rows_phase_and_version(tmp_path: Path) -> None:
    path = _write_atct(tmp_path, ATCT_ROWS + ATCT_ROWS[:1])  # repeated row
    values, info = builder.read_atct(path)

    assert info["version"] == "9.999"
    assert info["species_rows"] == 3
    by_name = {value.name: value for value in values}
    methane = by_name["Methane"]
    assert methane.raw_smiles == "C"
    assert methane.h298_kj == pytest.approx(-100.0)
    assert methane.uncertainty_kj == pytest.approx(0.40)
    assert methane.phase == "g"
    assert methane.source_id == "74-82-8*0"
    assert by_name["Methanol"].phase == "condensed"


def test_read_rmg_reads_every_tagged_value(tmp_path: Path) -> None:
    values, info = builder.read_rmg(_write_rmg(tmp_path))
    assert info["files"] == 3
    found = {(value.source, value.raw_smiles): value for value in values}
    assert set(found) == {
        ("rmg:CCCBDB", "C"),
        ("rmg:ATcT", "CCO"),
        ("rmg:CCCBDB", "CCO"),
    }
    assert found[("rmg:ATcT", "CCO")].uncertainty_kj == pytest.approx(0.8)


def test_read_bains_marks_filtered_rows(tmp_path: Path) -> None:
    values, info = builder.read_bains(_write_bains(tmp_path, BAINS_ROWS))
    assert info == {"rows": 6, "filter_flagged_rows": 1}
    acetone = [value for value in values if value.raw_smiles == "CC(C)=O"]
    assert len(acetone) == 1 and not acetone[0].selectable
    assert all(value.uncertainty_kj is None for value in values)
    sources = {value.source for value in values if value.raw_smiles == "CCO"}
    assert sources == {"bains:Pedley", "bains:Yaws"}


def _built(tmp_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    return builder.build(
        _write_atct(tmp_path), _write_rmg(tmp_path), _write_bains(tmp_path, BAINS_ROWS)
    )


def _row(table: pd.DataFrame, smiles: str) -> pd.Series:
    match = table[table["smiles"] == smiles]
    assert len(match) == 1, smiles
    return match.iloc[0]


def test_build_screens_and_counts_rejects(tmp_path: Path) -> None:
    table, _, info = _built(tmp_path)
    assert set(table["smiles"]) == {"C", "CCO", "CCCO", "O=CO"}
    assert info["rejected_values"] == {
        "open_shell": 1,
        "phase": 1,
        "elements": 1,
        "unparsable_smiles": 1,
    }
    assert table["exp_id"].is_unique
    assert (table["multiplicity"] == 1).all() and (table["charge"] == 0).all()


def test_selection_follows_source_rank_then_uncertainty(tmp_path: Path) -> None:
    table, _, _ = _built(tmp_path)

    methane = _row(table, "C")  # ATcT direct beats RMG CCCBDB
    assert methane["source"] == "atct"
    assert methane["H298_exp"] == pytest.approx(builder.kj_to_kcal(-100.0))
    assert methane["uncertainty"] == pytest.approx(builder.kj_to_kcal(0.40))
    assert methane["uncertainty_kind"] == "reported"

    # rmg:ATcT outranks rmg:CCCBDB even with a larger uncertainty.
    ethanol = _row(table, "CCO")
    assert ethanol["source"] == "rmg:ATcT"
    assert ethanol["n_values"] == 4 and ethanol["n_sources"] == 4


def test_low_priority_and_stewart_flags(tmp_path: Path) -> None:
    table, _, _ = _built(tmp_path)

    propanol = _row(table, "CCCO")  # only Yaws and Stewart
    assert propanol["source"] == "bains:Yaws"
    assert propanol["low_priority_only"]
    assert propanol["stewart_reference"]
    assert propanol["uncertainty_kind"] == "none" and math.isnan(
        propanol["uncertainty"]
    )

    formic = _row(table, "O=CO")  # Pedley beats Stewart; 12 kJ apart
    assert formic["source"] == "bains:Pedley"
    assert not formic["low_priority_only"]
    assert formic["stewart_reference"]
    assert formic["conflict_flag"] and formic["n_conflicts"] == 1
    assert formic["source_spread"] == pytest.approx(builder.kj_to_kcal(12.0))


def test_filtered_bains_row_is_never_selected(tmp_path: Path) -> None:
    table, long, _ = _built(tmp_path)
    assert "CC(C)=O" not in set(table["smiles"])
    assert "CC(C)=O" in set(long["smiles"])


def test_conflict_rule_uses_combined_uncertainty() -> None:
    chosen = pd.Series({"H298_exp": 0.0, "uncertainty": 1.0})
    others = pd.DataFrame(
        {"H298_exp": [2.5, 3.0, -0.9], "uncertainty": [1.0, 1.0, math.nan]}
    )
    # 2.5 and 3.0 against 2*hypot(1, 1) = 2.83; -0.9 against 2*hypot(1, 0) = 2.0.
    assert builder._conflicts(chosen, others) == 1


def test_compare_with_cbs_uses_lowest_conformer(tmp_path: Path) -> None:
    table, _, _ = _built(tmp_path)
    cbs = tmp_path / "cbs.csv"
    pd.DataFrame(
        {
            "mol_id": ["cbs_1", "cbs_2", "cbs_3", "cbs_4"],
            "smiles": ["C", "OCC", "CCO", "CCCC"],
            "multiplicity": [1, 1, 1, 1],
            "charge": [0, 0, 0, 0],
            "H298_cbs": [-20.0, -50.0, -52.0, -30.0],
        }
    ).to_csv(cbs, index=False)

    report = builder.compare_with_cbs(table, cbs)

    assert report["overlap_inchikey"] == 2
    methane = _row(table, "C")["H298_exp"]
    ethanol = _row(table, "CCO")["H298_exp"]
    errors = [-20.0 - methane, -52.0 - ethanol]
    assert report["all"]["bias"] == pytest.approx(sum(errors) / 2)
    assert report["all"]["mae"] == pytest.approx(sum(abs(e) for e in errors) / 2)


def test_main_writes_outputs_and_refuses_overwrite(tmp_path: Path) -> None:
    output = tmp_path / "out" / "exp.csv"
    argv = [
        "--atct-html",
        str(_write_atct(tmp_path)),
        "--rmg-dir",
        str(_write_rmg(tmp_path)),
        "--rmg-commit",
        "abc123",
        "--output",
        str(output),
    ]
    assert builder.main(argv) == 0

    manifest = json.loads(output.with_suffix(".manifest.json").read_text())
    assert manifest["rows"] == len(pd.read_csv(output))
    assert manifest["sources"]["atct"]["version"] == "9.999"
    assert any("abc123" in citation for citation in manifest["citations"])
    assert manifest["kj_per_kcal"] == builder.KJ_PER_KCAL
    assert (output.parent / "exp_all_sources.csv").exists()

    assert builder.main(argv) == 1
