"""CSV and minimal XLSX IO for the TCC workbook."""

from __future__ import annotations

import csv
import posixpath
import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from .models import PipelineTask

NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}

METHODS = {"AM1", "PM3", "PM7"}


def _cell_col(ref: str) -> int:
    letters = "".join(ch for ch in ref if ch.isalpha())
    value = 0
    for ch in letters:
        value = value * 26 + (ord(ch.upper()) - ord("A") + 1)
    return value - 1


def _coerce(value: str | None) -> str:
    return "" if value is None else str(value).strip()


def _float(value: str | None) -> float | None:
    text = _coerce(value)
    if not text:
        return None
    try:
        return float(text.replace(",", "."))
    except ValueError:
        return None


def _int(value: str | None, default: int = 0) -> int:
    number = _float(value)
    return default if number is None else int(number)


def _bool(value: str | None, default: bool = True) -> bool:
    text = _coerce(value).lower()
    if text in {"false", "falso", "0", "no", "nao", "não"}:
        return False
    if text in {"true", "verdadeiro", "1", "yes", "sim"}:
        return True
    return default


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


class MiniWorkbook:
    """Minimal read-only XLSX workbook reader for plain tables."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self._zip = ZipFile(path)
        self._shared_strings = self._read_shared_strings()
        self._sheets = self._read_sheet_paths()

    @property
    def sheet_names(self) -> list[str]:
        return list(self._sheets)

    def close(self) -> None:
        self._zip.close()

    def _read_shared_strings(self) -> list[str]:
        if "xl/sharedStrings.xml" not in self._zip.namelist():
            return []
        root = ET.fromstring(self._zip.read("xl/sharedStrings.xml"))
        values = []
        for si in root.findall("a:si", NS):
            values.append("".join(t.text or "" for t in si.findall(".//a:t", NS)))
        return values

    def _read_sheet_paths(self) -> dict[str, str]:
        workbook = ET.fromstring(self._zip.read("xl/workbook.xml"))
        rels = ET.fromstring(self._zip.read("xl/_rels/workbook.xml.rels"))
        relmap = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
        result: dict[str, str] = {}
        sheets = workbook.find("a:sheets", NS)
        if sheets is None:
            return result
        for sheet in sheets:
            name = sheet.attrib["name"]
            rid = sheet.attrib[
                "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
            ]
            target = relmap[rid]
            path = (
                target.lstrip("/")
                if target.startswith("/xl/")
                else posixpath.normpath(posixpath.join("xl", target))
            )
            result[name] = path
        return result

    def read_rows(self, sheet_name: str) -> list[list[str]]:
        if sheet_name not in self._sheets:
            raise KeyError(f"Sheet not found: {sheet_name}")
        root = ET.fromstring(self._zip.read(self._sheets[sheet_name]))
        rows = []
        for row in root.findall(".//a:sheetData/a:row", NS):
            values: list[str] = []
            for cell in row.findall("a:c", NS):
                idx = _cell_col(cell.attrib.get("r", "A1"))
                while len(values) <= idx:
                    values.append("")
                values[idx] = self._cell_value(cell)
            rows.append(values)
        return rows

    def table(self, sheet_name: str, required: set[str]) -> list[dict[str, str]]:
        rows = self.read_rows(sheet_name)
        header_idx = None
        headers: list[str] = []
        for idx, row in enumerate(rows):
            normalized = [_coerce(value) for value in row]
            if required.issubset(set(normalized)):
                header_idx = idx
                headers = normalized
                break
        if header_idx is None:
            raise ValueError(f"Required columns not found in sheet {sheet_name}")
        records = []
        for row in rows[header_idx + 1 :]:
            if not any(_coerce(value) for value in row):
                continue
            row = [*row, *([""] * (len(headers) - len(row)))]
            records.append({headers[i]: _coerce(row[i]) for i in range(len(headers))})
        return records

    def _cell_value(self, cell: ET.Element) -> str:
        if cell.attrib.get("t") == "inlineStr":
            return "".join(t.text or "" for t in cell.findall(".//a:t", NS))
        value = cell.find("a:v", NS)
        text = "" if value is None or value.text is None else value.text
        if cell.attrib.get("t") == "s" and text:
            return self._shared_strings[int(text)]
        return text


def read_xlsx_table(
    path: Path, sheet_name: str, required: set[str]
) -> list[dict[str, str]]:
    workbook = MiniWorkbook(path)
    try:
        return workbook.table(sheet_name, required)
    finally:
        workbook.close()


def load_completed_mol_ids(summary_path: Path) -> set[str]:
    """Read the summary xlsx and return mol_ids whose status is 'success'."""
    if not summary_path.exists():
        return set()
    try:
        workbook = MiniWorkbook(summary_path)
        try:
            rows = workbook.table(_SUMMARY_SHEET, {"mol_id", "status"})
        except (KeyError, ValueError):
            if not workbook.sheet_names:
                workbook.close()
                return set()
            rows = workbook.table(workbook.sheet_names[0], {"mol_id", "status"})
        finally:
            workbook.close()
        return {
            row["mol_id"]
            for row in rows
            if row.get("status", "").strip().lower() == "success"
        }
    except Exception:
        return set()


def load_molecule_reference(path: Path) -> dict[str, dict[str, str]]:
    rows = read_xlsx_table(path, "01_moleculas_cbthermo", {"mol_id", "smiles"})
    return {row["mol_id"]: row for row in rows}


def load_molecules(
    path: Path,
    limit: int | None = None,
) -> list[dict[str, object]]:
    """Load unique molecules from the workbook (one entry per molecule, no method duplication)."""
    rows = read_xlsx_table(path, "01_moleculas_cbthermo", {"mol_id", "smiles"})
    molecules: list[dict[str, object]] = []
    for row in rows:
        molecules.append(
            {
                "mol_id": row.get("mol_id", ""),
                "molecule": row.get("molecule", ""),
                "smiles": row.get("smiles", ""),
                "hf_exp_kJmol": _float(row.get("Hf_exp_kJmol")),
                "formula": row.get("formula", ""),
                "group": row.get("group", ""),
                "charge": _int(row.get("charge"), 0),
                "multiplicity": _int(row.get("multiplicity"), 1),
                "run_crest": _bool(row.get("run_crest"), True),
            }
        )
        if limit is not None and len(molecules) >= limit:
            break
    return molecules


_SUMMARY_SHEET = "grimperium_mini_summary"


def append_summary_row(
    path: Path, row: dict[str, object], columns: list[str]
) -> None:
    """Append one result row to the summary xlsx, creating the file if needed.

    If the file exists but has a different schema (e.g., from a prior pipeline
    version), it is silently overwritten starting from this row.
    """
    rows_so_far: list[dict[str, object]] = []
    if path.exists():
        try:
            prior = read_xlsx_table(path, _SUMMARY_SHEET, set())
            rows_so_far = [dict(r) for r in prior]
        except KeyError:
            pass  # existing file has a different sheet layout; start fresh
    rows_so_far.append(row)
    write_simple_xlsx(path, {_SUMMARY_SHEET: (columns, rows_so_far)})


def load_pipeline_tasks(
    path: Path,
    methods: list[str] | None = None,
    limit: int | None = None,
) -> list[PipelineTask]:
    selected_methods = {m.upper() for m in methods} if methods else METHODS
    refs = load_molecule_reference(path) if path.suffix.lower() == ".xlsx" else {}
    if path.suffix.lower() == ".csv":
        rows = read_csv_rows(path)
    else:
        rows = read_xlsx_table(
            path, "02_pipeline_grimperium_mini", {"mol_id", "smiles", "method"}
        )

    tasks: list[PipelineTask] = []
    for row in rows:
        method = row.get("method", "").upper()
        if method not in selected_methods or method not in METHODS:
            continue
        ref = refs.get(row.get("mol_id", ""), {})
        max_conf = _int(row.get("max_conformers_to_process"), 0) or None
        tasks.append(
            PipelineTask(
                mol_id=row.get("mol_id", ""),
                molecule=row.get("molecule", "") or ref.get("molecule", ""),
                smiles=row.get("smiles", "") or ref.get("smiles", ""),
                method=method,  # type: ignore[arg-type]
                run_crest=_bool(row.get("run_crest"), True),
                charge=_int(row.get("charge"), 0),
                multiplicity=_int(row.get("multiplicity"), 1),
                formula=ref.get("formula", ""),
                group=ref.get("group", ""),
                max_conformers_to_process=max_conf,
                hf_exp_kJmol=_float(
                    row.get("ref_Hf_exp_kJmol") or ref.get("Hf_exp_kJmol")
                ),
                hf_tcc_am1_kJmol=_float(ref.get("Hf_AM1_TCC_VEGA_MOPAC_kJmol")),
                hf_tcc_pm3_kJmol=_float(ref.get("Hf_PM3_TCC_VEGA_MOPAC_kJmol")),
                hf_tcc_pm7_kJmol=_float(
                    row.get("ref_Hf_PM7_TCC_VEGA_MOPAC_kJmol")
                    or ref.get("Hf_PM7_TCC_VEGA_MOPAC_kJmol")
                ),
                hf_joback_kJmol=_float(
                    row.get("ref_Hf_JR_kJmol") or ref.get("Hf_JR_kJmol")
                ),
                hf_constantinou_gani_kJmol=_float(
                    row.get("ref_Hf_CG_kJmol") or ref.get("Hf_CG_kJmol")
                ),
                hf_marrero_gani_kJmol=_float(
                    row.get("ref_Hf_MG_kJmol") or ref.get("Hf_MG_kJmol")
                ),
            )
        )
        if limit is not None and len(tasks) >= limit:
            break
    return tasks


def validate_workbook(path: Path) -> list[str]:
    workbook = MiniWorkbook(path)
    try:
        missing = []
        required_sheets = {
            "01_moleculas_cbthermo",
            "02_pipeline_grimperium_mini",
            "03_resultados_formacao",
            "04_reacoes",
            "05_resumo_metricas",
            "06_grimperium_original",
        }
        for sheet in sorted(required_sheets - set(workbook.sheet_names)):
            missing.append(f"missing sheet: {sheet}")
        molecules = workbook.table("01_moleculas_cbthermo", {"mol_id", "smiles"})
        tasks = workbook.table(
            "02_pipeline_grimperium_mini", {"mol_id", "smiles", "method"}
        )
        reactions = workbook.table(
            "04_reacoes",
            {"reaction_id", "coef_a", "coef_b", "coef_c", "coef_d"},
        )
        if not molecules:
            missing.append("01_moleculas_cbthermo has no data rows")
        if not tasks:
            missing.append("02_pipeline_grimperium_mini has no data rows")
        bad_methods = sorted(
            {
                r.get("method", "")
                for r in tasks
                if r.get("method", "").upper() not in METHODS
            }
        )
        if bad_methods:
            missing.append(
                f"unsupported methods in pipeline sheet: {', '.join(bad_methods)}"
            )
        if not reactions:
            missing.append("04_reacoes has no data rows")
        return missing
    finally:
        workbook.close()


def _xml_escape(value: object) -> str:
    text = "" if value is None else str(value)
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_simple_xlsx(
    path: Path, sheets: dict[str, tuple[list[str], Sequence[Mapping[str, object]]]]
) -> None:
    """Write a small XLSX workbook with inline strings."""
    path.parent.mkdir(parents=True, exist_ok=True)
    content_types = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
    ]
    for idx in range(1, len(sheets) + 1):
        content_types.append(
            f'<Override PartName="/xl/worksheets/sheet{idx}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )
    content_types.append("</Types>")
    workbook_sheets = []
    rels = []
    worksheet_files = {}
    for idx, (name, (columns, rows)) in enumerate(sheets.items(), start=1):
        workbook_sheets.append(
            f'<sheet name="{_xml_escape(name)}" sheetId="{idx}" r:id="rId{idx}"/>'
        )
        rels.append(
            f'<Relationship Id="rId{idx}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{idx}.xml"/>'
        )
        worksheet_files[f"xl/worksheets/sheet{idx}.xml"] = _worksheet_xml(columns, rows)
    with ZipFile(path, "w", ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", "\n".join(content_types))
        z.writestr(
            "_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>',
        )
        z.writestr(
            "xl/workbook.xml",
            '<?xml version="1.0" encoding="UTF-8"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
            + "".join(workbook_sheets)
            + "</sheets></workbook>",
        )
        z.writestr(
            "xl/_rels/workbook.xml.rels",
            '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            + "".join(rels)
            + "</Relationships>",
        )
        for filename, xml in worksheet_files.items():
            z.writestr(filename, xml)


def _worksheet_xml(columns: list[str], rows: Sequence[Mapping[str, object]]) -> str:
    all_rows: list[Mapping[str, object]] = [dict(zip(columns, columns)), *rows]
    xml_rows = []
    for r_idx, row in enumerate(all_rows, start=1):
        cells = []
        for c_idx, column in enumerate(columns):
            ref = f"{_column_name(c_idx)}{r_idx}"
            value = row.get(column, "")
            if isinstance(value, int | float) and value != "":
                cells.append(f'<c r="{ref}"><v>{value}</v></c>')
            else:
                cells.append(
                    f'<c r="{ref}" t="inlineStr"><is><t>{_xml_escape(value)}</t></is></c>'
                )
        xml_rows.append(f'<row r="{r_idx}">' + "".join(cells) + "</row>")
    return (
        '<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'
        + "".join(xml_rows)
        + "</sheetData></worksheet>"
    )


def _column_name(idx: int) -> str:
    name = ""
    idx += 1
    while idx:
        idx, rem = divmod(idx - 1, 26)
        name = chr(ord("A") + rem) + name
    return name
