from __future__ import annotations

import hashlib
import stat
from pathlib import Path

import pytest

from semi_imperium.training.program_info import probe_programs


def _script(folder: Path, name: str, output: str) -> Path:
    path = folder / name
    path.write_text(f"#!/bin/sh\ncat <<'EOF'\n{output}\nEOF\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_version_hash_and_path_of_each_program(tmp_path: Path) -> None:
    crest = _script(tmp_path, "crest", "  banner\n       Version 3.0.2, Tue, 05 August")
    xtb = _script(tmp_path, "xtb", "   * xtb version 6.7.1 (edcfbbe) compiled by")
    mopac = _script(tmp_path, "mopac", "MOPAC version 23.2.2")

    info = probe_programs({"crest": str(crest), "xtb": str(xtb), "mopac": str(mopac)})

    assert info == {
        "crest": "3.0.2",
        "crest_sha256": _sha256(crest),
        "crest_path": str(crest),
        "xtb": "6.7.1",
        "xtb_sha256": _sha256(xtb),
        "xtb_path": str(xtb),
        "mopac": "23.2.2",
        "mopac_sha256": _sha256(mopac),
        "mopac_path": str(mopac),
    }


def test_a_bare_name_is_resolved_on_the_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    xtb = _script(tmp_path, "xtb", "xtb version 6.7.1")
    monkeypatch.setenv("PATH", str(tmp_path))

    info = probe_programs({"xtb": "xtb"})

    assert info["xtb_path"] == str(xtb) and info["xtb_sha256"] == _sha256(xtb)


def test_a_missing_program_is_recorded_as_unavailable(tmp_path: Path) -> None:
    info = probe_programs({"mopac": str(tmp_path / "nowhere")})

    assert info == {"mopac": "unavailable", "mopac_path": str(tmp_path / "nowhere")}


def test_unrecognised_version_output_keeps_the_hash(tmp_path: Path) -> None:
    crest = _script(tmp_path, "crest", "no version line here")

    info = probe_programs({"crest": str(crest)})

    assert info["crest"] == "unknown" and info["crest_sha256"] == _sha256(crest)
