"""Unit repair of the ThermoCBS enthalpy columns in the CHON builder."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest


def _load_builder() -> ModuleType:
    repo_root = Path(__file__).resolve().parents[2]
    script_path = repo_root / "scripts" / "build_thermo_cbs_chon.py"
    spec = importlib.util.spec_from_file_location("build_thermo_cbs_chon", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


builder = _load_builder()

WATER_XYZ = "3\ncomment\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0"
SULFUR_XYZ = "3\ncomment\nS 0 0 0\nH 1.3 0 0\nH -0.3 1.3 0"

# Values as written in the source: integers are cal/mol.
SOURCE = f""",smiles,xyz,multiplicity,charge,nheavy,H298_cbs,H298_b3,cbs_b3
0,CCO,"{WATER_XYZ}",1,0,3,-54.11786,-40.0,-14.11786
1,C=CCC(=O)O,"{WATER_XYZ}",1,0,6,-79523,-64.74815,-14.77485
2,CC(=O)NC(C(N)=O)C(N)=O,"{WATER_XYZ}",1,0,9,-145781,-124351,-21.43
3,CCN,"{WATER_XYZ}",1,0,3,-11.0,-5.0,-6000
4,S,"{SULFUR_XYZ}",1,0,1,-4.9,-1.0,-3.9
5,CO,"{WATER_XYZ}",1,0,2,-48.0,-30.0,-18.0
6,OCO,"{WATER_XYZ}",1,0,3,-90.5,-60.25,-30.25
"""


def _source(tmp_path: Path, text: str = SOURCE) -> Path:
    path = tmp_path / "thermo_cbs.csv"
    path.write_text(text, encoding="utf-8")
    return path


def test_integer_formatted_enthalpies_are_converted_from_cal(tmp_path: Path) -> None:
    table = builder.build(_source(tmp_path), None).set_index("mol_id")

    assert table.loc["cbs_00000", "H298_cbs"] == pytest.approx(-54.11786)
    assert table.loc["cbs_00001", "H298_cbs"] == pytest.approx(-79.523)
    assert table.loc["cbs_00002", "H298_b3"] == pytest.approx(-124.351)
    assert table.loc["cbs_00003", "cbs_b3"] == pytest.approx(-6.0)
    assert table.loc["cbs_00000", "h298_unit_repair"] == ""
    assert table.loc["cbs_00001", "h298_unit_repair"] == "H298_cbs"
    assert table.loc["cbs_00002", "h298_unit_repair"] == "H298_cbs;H298_b3"
    assert table.loc["cbs_00003", "h298_unit_repair"] == "cbs_b3"
    residual = table["H298_cbs"] - table["H298_b3"] - table["cbs_b3"]
    assert residual.abs().max() < builder.IDENTITY_TOLERANCE
    assert "cbs_00004" not in table.index  # sulfur still filtered out


def test_identity_failure_stops_the_build(tmp_path: Path) -> None:
    broken = SOURCE.replace("-48.0,-30.0,-18.0", "-48.0,-30.0,-1.0")

    with pytest.raises(ValueError, match="H298_cbs = H298_b3 \\+ cbs_b3"):
        builder.build(_source(tmp_path, broken), None)


def test_text_format_decides_the_unit() -> None:
    assert builder.is_integer_text("-79523")
    assert not builder.is_integer_text("-79523.0")
    assert not builder.is_integer_text("1e-3")
    frame = pd.DataFrame(
        {
            "H298_cbs": ["-17307", "2.5"],
            "H298_b3": ["3.22464", "1.0"],
            "cbs_b3": ["-20.53164", "1.5"],
        }
    )
    repaired = builder.repair_enthalpy_units(frame)
    assert repaired["H298_cbs"].tolist() == pytest.approx([-17.307, 2.5])
    assert repaired["h298_unit_repair"].tolist() == ["H298_cbs", ""]
