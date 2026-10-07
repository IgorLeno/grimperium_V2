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
    '<button onclick="copyname(\'{name}\');" type="button">{formula}  ({label}) </button>'
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
        "label": "g",
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
        "label": "g",
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
        "label": "cr,l",
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


def test_atct_charge_comes_from_row_id(tmp_path: Path) -> None:
    cation = {
        **ATCT_ROWS[0],
        "rid": "s1_6n4_1c1",
        "number": "99",
        "name": "Methane cation",
        "h298": "1100.0",
    }
    values, _ = builder.read_atct(_write_atct(tmp_path, [ATCT_ROWS[0], cation]))
    assert {value.declared_charge for value in values} == {0, 1}

    long, rejected = builder.long_table(values)
    assert set(long["name"]) == {"Methane"}
    assert rejected == {"charged": 1}


def test_same_source_isomers_keep_lowest_enthalpy() -> None:
    values = [
        builder.SourceValue(
            "atct", "628-92-2*0", "Cycloheptene", "C1=CCCCCC1", -8.7, 0.7, "g"
        ),
        builder.SourceValue(
            "atct", "45509-99-7*0", "trans-Cycloheptene", "C1=CCCCCC1", 103.9, 0.4, "g"
        ),
        builder.SourceValue(
            "bains:Pedley", "row1", "Cycloheptene", "C1=CCCCCC1", -9.0, None, "g"
        ),
    ]
    long, _ = builder.long_table(values)
    row = builder.select(long).iloc[0]

    assert row["source_id"] == "628-92-2*0"
    assert row["isomers_collapsed"] == 1
    assert row["n_values"] == 2
    assert not row["conflict_flag"]
    assert row["source_spread"] == pytest.approx(builder.kj_to_kcal(0.3))


def test_rmg_zero_uncertainty_is_kept(tmp_path: Path) -> None:
    directory = tmp_path / "rmg0"
    directory.mkdir()
    (directory / "Methane.yml").write_text(
        "smiles: C\nreference_data:\n  ATcT:\n    thermo_data:\n"
        "      H298: {units: kJ/mol, uncertainty: 0.0, value: -74.5}\n"
    )
    values, _ = builder.read_rmg(directory)
    assert values[0].uncertainty_kj == 0.0


@pytest.mark.parametrize(
    ("formula", "atct_id", "expected"),
    [
        ("CH4  (g) ", "74-82-8*0", "g"),
        ("CH3COOH  (g, syn) ", "64-19-7*1", "g_variant"),
        ("CH2  (g, singlet) ", "2465-56-7*2", "g_variant"),
        ("C6H6  (cr,l) ", "71-43-2*500", "condensed"),
        ("NH2CH2COOH  (aq) ", "56-40-6*800", "condensed"),
    ],
)
def test_atct_phase_from_formula_label(
    formula: str, atct_id: str, expected: str
) -> None:
    assert builder.atct_phase(formula, atct_id) == expected


# Synthetic values in the NIST WebBook compound-page layout (not NIST data).
NIST_ONE_D_ROW = (
    '<tr class="cal"><td style="text-align: left;">{quantity}</td>'
    '<td class="right-nowrap">{value}</td><td style="text-align: right;">'
    '{units}</td><td style="text-align: center;"><a href="#">{method}</a></td>'
    '<td style="text-align: left;">Ref</td><td style="text-align: left;">x</td></tr>'
)
NIST_ONE_D_HEAD = (
    '<table class="data" aria-label="One dimensional data"><tr>'
    '<th scope="col">Quantity</th><th scope="col">Value</th>'
    '<th scope="col">Units</th><th scope="col">Method</th>'
    '<th scope="col">Reference</th><th scope="col">Comment</th></tr>'
)
NIST_T_TABLE = (
    '<h3>Enthalpy of {what}</h3><table class="data" aria-label="Enthalpy of '
    '{what}"><tr><th scope="col">&#916;<sub>{sub}</sub>H (kJ/mol)</th>'
    '<th scope="col">Temperature (K)</th><th scope="col">Method</th>'
    '<th scope="col">Reference</th><th scope="col">Comment</th></tr>{rows}</table>'
)
ETHANOL_INCHI = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
ETHANOL_KEY = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


def _q(kind: str, phase: str = "") -> str:
    """``ΔfH°gas`` or ``ΔvapH°`` as the WebBook writes them."""
    if kind == "f":
        return f"&#916;<sub>f</sub>H&deg;<sub>{phase}</sub>"
    return f"&#916;<sub>{kind}</sub>H&deg;"


def _nist_section(section: str, rows: list[tuple[str, str, str]]) -> str:
    body = "".join(
        NIST_ONE_D_ROW.format(quantity=q, value=v, units="kJ/mol", method=m)
        for q, v, m in rows
    )
    return f'<h2 id="{section}">x</h2>{NIST_ONE_D_HEAD}{body}</table>'


def _nist_page(
    *,
    formula: str = "C<sub>2</sub>H<sub>6</sub>O",
    inchi: str = ETHANOL_INCHI,
    inchikey: str = ETHANOL_KEY,
    gas: list[tuple[str, str, str]] | None = None,
    condensed: list[tuple[str, str, str]] | None = None,
    phase: list[tuple[str, str, str]] | None = None,
    vap_table: list[tuple[str, str]] | None = None,
) -> str:
    parts = [
        "<html><body><h1>NIST Chemistry WebBook, SRD 69</h1>"
        '<main id="main"><h1 id="Top">Ethanol</h1><ul>'
        '<li><strong><a href="#">Formula</a>:</strong> ' + formula + "</li>"
        '<li><strong><a href="#">Molecular weight</a>:</strong> 46.0684</li>'
        '<li><div><strong>IUPAC Standard InChI:</strong> <span class="inchi-text">'
        + inchi
        + "</span></div></li>"
        '<li><div><strong>IUPAC Standard InChIKey:</strong> <span class="inchi-text">'
        + inchikey
        + "</span></div></li></ul>"
    ]
    if gas is not None:
        parts.append(_nist_section("Thermo-Gas", gas))
    if condensed is not None:
        parts.append(_nist_section("Thermo-Condensed", condensed))
    if phase is not None or vap_table is not None:
        parts.append(_nist_section("Thermo-Phase", phase or []))
        if vap_table is not None:
            rows = "".join(
                f'<tr class="exp"><td class="right-nowrap">{v}</td>'
                f'<td class="right-nowrap">{t}</td><td>N/A</td><td>Ref</td>'
                "<td>&nbsp;</td></tr>"
                for v, t in vap_table
            )
            parts.append(NIST_T_TABLE.format(what="vaporization", sub="vap", rows=rows))
    parts.append(
        '<h2 id="Notes">Notes</h2><table class="data"><tr><td>'
        + _q("f", "gas")
        + "</td><td>Enthalpy of formation of gas at standard conditions</td></tr>"
        "</table></main></body></html>"
    )
    return "".join(parts)


def _write_nist(tmp_path: Path, pages: dict[str, str]) -> Path:
    import gzip

    directory = tmp_path / "nist"
    (directory / "compound").mkdir(parents=True)
    for compound_id, page in pages.items():
        with gzip.open(directory / "compound" / f"{compound_id}.html.gz", "wt") as fh:
            fh.write(page)
    (directory / "fetch_log.jsonl").write_text(
        '{"url": "u", "status": 200, "time": "2026-10-07T01:00:00+00:00"}\n'
        '{"url": "u", "status": 200, "time": "2026-10-08T02:00:00+00:00"}\n'
    )
    return directory


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("-234. ± 2.", (-234.0, 2.0)),
        ("-323.6", (-323.6, None)),
        ("95.5 ± 0.3", (95.5, 0.3)),
    ],
)
def test_parse_nist_value(text: str, expected: tuple[float, float | None]) -> None:
    assert builder.parse_nist_value(text) == expected


def test_pick_nist_prefers_average_then_uncertainty_then_median() -> None:
    average = builder.pick_nist([(-1.0, 0.1, "Ccb"), (-2.0, 2.0, "AVG")])
    assert (average.value_kj, average.uncertainty_kj, average.n_rows) == (-2.0, 2.0, 2)
    precise = builder.pick_nist(
        [(-1.0, 0.5, "Ccb"), (-3.0, 0.2, "Ccb"), (-9.0, None, "N/A")]
    )
    assert precise.value_kj == -3.0
    median = builder.pick_nist(
        [(-5.0, None, "a"), (-1.0, None, "b"), (-3.0, None, "c")]
    )
    assert (median.value_kj, median.uncertainty_kj) == (-3.0, None)


def test_parse_nist_page_reads_identity_and_average() -> None:
    page = builder.parse_nist_page(
        _nist_page(gas=[(_q("f", "gas"), "-200. &plusmn; 2.", "AVG")]), "C1"
    )
    assert page.name == "Ethanol"
    assert page.formula == "C2H6O"
    assert (page.inchi, page.inchikey) == (ETHANOL_INCHI, ETHANOL_KEY)
    assert page.picks["gas"] == builder.NistPick(-200.0, 2.0, "AVG", 1)
    assert set(page.picks) == {"gas"}  # the Notes symbol table is not data


def test_nist_condensed_plus_vaporization_at_298() -> None:
    html = _nist_page(
        gas=[(_q("f", "gas"), "-200. &plusmn; 2.", "AVG")],
        condensed=[(_q("f", "liquid"), "-240. &plusmn; 3.", "AVG")],
        phase=[(_q("vap"), "41. &plusmn; 4.", "AVG")],
    )
    page = builder.parse_nist_page(html, "C1")
    smiles, reason = builder.nist_smiles(page)
    assert (smiles, reason) == ("CCO", "")

    values = {v.source: v for v in builder.nist_values(page, smiles)}
    assert set(values) == {"nist:gas", "nist:liquid+vap"}
    derived = values["nist:liquid+vap"]
    assert derived.h298_kj == pytest.approx(-199.0)
    assert derived.uncertainty_kj == pytest.approx(5.0)
    assert all(v.validation_only and v.phase == "g" for v in values.values())


def test_nist_vaporization_off_298_is_ignored() -> None:
    html = _nist_page(
        condensed=[(_q("f", "liquid"), "-240. &plusmn; 3.", "AVG")],
        vap_table=[("38.6", "351.5"), ("40.0", "326.")],
    )
    page = builder.parse_nist_page(html, "C1")
    assert "vap" not in page.picks
    assert page.off_298_rows == 2
    assert builder.nist_values(page, "CCO") == []


def test_nist_vaporization_table_row_at_298_is_used() -> None:
    html = _nist_page(
        condensed=[(_q("f", "liquid"), "-240.", "Ccb")],
        vap_table=[("42.0", "298."), ("38.6", "351.5")],
    )
    page = builder.parse_nist_page(html, "C1")
    (value,) = builder.nist_values(page, "CCO")
    assert value.source == "nist:liquid+vap"
    assert value.h298_kj == pytest.approx(-198.0)
    assert value.uncertainty_kj is None


def test_nist_compound_without_gas_or_route_gives_nothing(tmp_path: Path) -> None:
    pages = {
        "C1": _nist_page(condensed=[(_q("f", "liquid"), "-240.", "Ccb")]),
        "C2": _nist_page(inchi="", inchikey=""),
    }
    values, info = builder.read_nist(_write_nist(tmp_path, pages))
    assert values == []
    assert info["identity_rejects"] == {"no_usable_value": 1, "no_inchi": 1}


@pytest.mark.parametrize(
    ("kwargs", "reason"),
    [
        ({"formula": "C<sub>2</sub>H<sub>4</sub>O"}, "formula_mismatch"),
        ({"inchikey": "IKHGUXGNUITLKF-UHFFFAOYSA-N"}, "inchikey_mismatch"),
        ({"inchi": "InChI=1S/garbage"}, "unparsable_inchi"),
    ],
)
def test_nist_identity_must_match_page_metadata(
    kwargs: dict[str, str], reason: str
) -> None:
    page = builder.parse_nist_page(_nist_page(**kwargs), "C1")
    assert builder.nist_smiles(page) == (None, reason)


def test_nist_ranks_below_trainable_sources_and_is_validation_only(
    tmp_path: Path,
) -> None:
    ethanol = _nist_page(gas=[(_q("f", "gas"), "-229.5 &plusmn; 0.1", "AVG")])
    propanol = _nist_page(
        formula="C<sub>3</sub>H<sub>8</sub>O",
        inchi="InChI=1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3",
        inchikey="BDERNNFJNOPAEC-UHFFFAOYSA-N",
        gas=[(_q("f", "gas"), "-250. &plusmn; 1.", "AVG")],
    )
    nist_dir = _write_nist(tmp_path, {"C64175": ethanol, "C71238": propanol})
    bains = _write_bains(tmp_path, BAINS_ROWS)

    table, long, info = builder.build(None, None, bains, nist_dir)
    by_smiles = table.set_index("smiles")
    # Bains Pedley (trainable) beats NIST; NIST beats the low-priority Yaws.
    assert by_smiles.loc["CCO", "source"] == "bains:Pedley"
    assert not by_smiles.loc["CCO", "validation_only"]
    assert by_smiles.loc["CCCO", "source"] == "nist:gas"
    assert bool(by_smiles.loc["CCCO", "validation_only"])
    assert not by_smiles.loc["CCCO", "low_priority_only"]
    assert long.loc[long["source"] == "nist:gas", "validation_only"].all()
    assert info["nist"]["retrieved"] == ["2026-10-07", "2026-10-08"]
    rank = builder.SOURCE_RANK
    assert rank["bains:Winget"] < rank["nist:gas"] < rank["nist:liquid+vap"]


def test_nist_ion_gas_value_is_lowest_priority_fallback() -> None:
    mixed = _nist_page(
        gas=[
            (_q("f", "gas"), "-230. &plusmn; 9.", "Ion"),
            (_q("f", "gas"), "-234. &plusmn; 2.", "Ccb"),
        ]
    )
    page = builder.parse_nist_page(mixed, "C1")
    assert set(page.picks) == {"gas"}
    assert page.picks["gas"].value_kj == -234.0

    ion_only = _nist_page(gas=[(_q("f", "gas"), "-230. &plusmn; 9.", "Ion")])
    page = builder.parse_nist_page(ion_only, "C1")
    (value,) = builder.nist_values(page, "CCO")
    assert value.source == "nist:gas_ion"
    assert builder.SOURCE_RANK["nist:gas_ion"] >= builder.LOW_PRIORITY_RANK


def test_main_cites_nist_terms(tmp_path: Path) -> None:
    gas_only = _nist_page(gas=[(_q("f", "gas"), "-229.5 &plusmn; 0.1", "AVG")])
    output = tmp_path / "exp.csv"
    argv = ["--nist-dir", str(_write_nist(tmp_path, {"C1": gas_only}))]
    assert builder.main([*argv, "--output", str(output)]) == 0

    manifest = json.loads(output.with_suffix(".manifest.json").read_text())
    assert manifest["validation_only_rows"] == 1
    assert any(
        "Number 69" in c and "2026-10-07 to 2026-10-08" in c
        for c in manifest["citations"]
    )
    assert "validation_only=True" in manifest["license_note"]


def test_nist_value_without_uncertainty_is_never_selected(tmp_path: Path) -> None:
    bare = _nist_page(gas=[(_q("f", "gas"), "-229.5", "Ccb")])
    values, _ = builder.read_nist(_write_nist(tmp_path, {"C1": bare}))
    assert [(v.source, v.selectable) for v in values] == [("nist:gas", False)]

    methane = builder.SourceValue("atct", "74-82-8*0", "Methane", "C", -74.5, 0.1, "g")
    long, _ = builder.long_table([*values, methane])
    assert len(long) == 2
    assert list(builder.select(long)["smiles"]) == ["C"]
