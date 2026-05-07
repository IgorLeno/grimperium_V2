from pathlib import Path

from grimperium_mini.mopac import create_mopac_input


def _xyz(path: Path) -> None:
    path.write_text(
        "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 1.0 0.0 0.0\n",
        encoding="utf-8",
    )


def test_mopac_input_uses_requested_method(tmp_path):
    xyz = tmp_path / "water.xyz"
    _xyz(xyz)
    for method in ("AM1", "PM3", "PM7"):
        mop = tmp_path / f"{method}.mop"
        create_mopac_input(xyz, mop, method)  # type: ignore[arg-type]
        text = mop.read_text(encoding="utf-8")
        assert text.startswith(f"{method} PRECISE EF AUX\n")
        assert "PM7 PRECISE" not in text if method != "PM7" else True


def test_mopac_input_includes_charge_and_multiplicity(tmp_path):
    xyz = tmp_path / "water.xyz"
    _xyz(xyz)
    mop = tmp_path / "charged.mop"
    create_mopac_input(xyz, mop, "AM1", charge=-1, multiplicity=2)
    assert mop.read_text(encoding="utf-8").splitlines()[0] == (
        "AM1 PRECISE EF AUX CHARGE=-1 DOUBLET"
    )
