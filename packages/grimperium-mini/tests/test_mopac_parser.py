from types import SimpleNamespace

from grimperium_mini.config import MiniConfig
from grimperium_mini.mopac import parse_hof_kcalmol, run_mopac


def test_parse_final_heat_of_formation():
    text = "FINAL HEAT OF FORMATION =      -17.93603 KCAL/MOL = -75.04 KJ/MOL"
    assert parse_hof_kcalmol(text) == -17.93603


def test_parse_heat_of_formation_without_final():
    text = "HEAT OF FORMATION       =        12.500 KCAL/MOL"
    assert parse_hof_kcalmol(text) == 12.5


def test_nonzero_mopac_with_extractable_hof_is_warning(tmp_path, monkeypatch):
    xyz = tmp_path / "c.xyz"
    xyz.write_text("1\nC\nC 0 0 0\n", encoding="utf-8")

    def fake_run(cmd, cwd, capture_output, text, timeout):
        out = tmp_path / "mopac" / "conf_PM7.out"
        out.write_text("FINAL HEAT OF FORMATION = -10.0 KCAL/MOL\n", encoding="utf-8")
        return SimpleNamespace(returncode=7, stdout="", stderr="warning")

    monkeypatch.setattr("grimperium_mini.mopac.subprocess.run", fake_run)
    result = run_mopac(xyz, tmp_path / "mopac", "PM7", MiniConfig(), stem="conf")
    assert result.status == "warning"
    assert result.hof_kJmol == -41.84


def test_mopac_failure_without_hof_is_failed(tmp_path, monkeypatch):
    xyz = tmp_path / "c.xyz"
    xyz.write_text("1\nC\nC 0 0 0\n", encoding="utf-8")

    def fake_run(cmd, cwd, capture_output, text, timeout):
        out = tmp_path / "mopac" / "conf_AM1.out"
        out.write_text("NO HOF HERE\n", encoding="utf-8")
        return SimpleNamespace(returncode=1, stdout="", stderr="bad")

    monkeypatch.setattr("grimperium_mini.mopac.subprocess.run", fake_run)
    result = run_mopac(xyz, tmp_path / "mopac", "AM1", MiniConfig(), stem="conf")
    assert result.status == "failed"
    assert "Heat of formation not found" in result.error_message
