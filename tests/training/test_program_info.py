from __future__ import annotations

import hashlib
import stat
from pathlib import Path

import pytest

from semi_imperium.training.program_info import probe_programs, program_mismatches


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
        "environment_sha256": "unavailable",
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

    assert info == {
        "mopac": "unavailable",
        "mopac_path": str(tmp_path / "nowhere"),
        "environment_sha256": "unavailable",
    }


def test_unrecognised_version_output_keeps_the_hash(tmp_path: Path) -> None:
    crest = _script(tmp_path, "crest", "no version line here")

    info = probe_programs({"crest": str(crest)})

    assert info["crest"] == "unknown" and info["crest_sha256"] == _sha256(crest)


def _env(root: Path, packages: list[str]) -> Path:
    """A conda-style environment: ``bin/`` with the programs, ``conda-meta/``."""
    (root / "bin").mkdir(parents=True)
    meta = root / "conda-meta"
    meta.mkdir()
    (meta / "history").write_text("==> 2026-10-20 <==\n")
    for package in packages:
        (meta / f"{package}.json").write_text("{}")
    return root / "bin"


LOCKED = ["crest-3.0.2-h8876d29_1", "libopenblas-0.3.34-pthreads_hf13c14d_2"]


def test_the_environment_is_hashed_by_its_package_list(tmp_path: Path) -> None:
    bin_a = _env(tmp_path / "a", LOCKED)
    bin_b = _env(tmp_path / "b", list(reversed(LOCKED)))
    for folder in (bin_a, bin_b):
        _script(folder, "crest", "Version 3.0.2")
        _script(folder, "xtb", "xtb version 6.7.1")

    a = probe_programs({"crest": str(bin_a / "crest"), "xtb": str(bin_a / "xtb")})
    b = probe_programs({"crest": str(bin_b / "crest"), "xtb": str(bin_b / "xtb")})

    # Same lockfile in another folder: same hash, so the cell matches.
    assert a["environment_sha256"] == b["environment_sha256"]
    assert (
        a["environment_sha256"]
        == hashlib.sha256(
            "\n".join(sorted(f"{p}.json" for p in LOCKED)).encode()
        ).hexdigest()
    )
    assert program_mismatches(a, b) == []


def test_another_library_build_changes_the_environment_hash(tmp_path: Path) -> None:
    ours = _env(tmp_path / "a", LOCKED)
    theirs = _env(
        tmp_path / "b", ["crest-3.0.2-h8876d29_1", "libopenblas-0.3.30-pthreads_0"]
    )
    for folder in (ours, theirs):
        _script(folder, "crest", "Version 3.0.2")

    central = probe_programs({"crest": str(ours / "crest")})
    cell = probe_programs({"crest": str(theirs / "crest")})

    # The binary is the same; only a library the hash of crest misses differs.
    assert central["crest_sha256"] == cell["crest_sha256"]
    assert [m.split(":")[0] for m in program_mismatches(central, cell)] == [
        "environment_sha256"
    ]


CENTRAL = {
    "crest": "3.0.2",
    "crest_sha256": "aaa",
    "crest_path": "/home/central/env/bin/crest",
    "mopac": "23.2.2",
    "mopac_sha256": "bbb",
    "mopac_path": "/x/mopac",
}


def test_same_programs_in_another_folder_match() -> None:
    cell = {**CENTRAL, "crest_path": "/home/aluno/env/bin/crest", "mopac_path": "m"}
    assert program_mismatches(CENTRAL, cell) == []


def test_version_hash_and_missing_programs_are_each_reported() -> None:
    cell = {**CENTRAL, "crest": "3.0.1", "mopac_sha256": "ccc"}
    del cell["crest_sha256"]
    assert program_mismatches(CENTRAL, cell) == [
        "crest: 3.0.1 (central 3.0.2)",
        "crest_sha256: ausente (central aaa)",
        "mopac_sha256: ccc (central bbb)",
    ]


def test_unreadable_programs_never_match_even_on_both_sides() -> None:
    central = {**CENTRAL, "crest": "unknown"}
    cell = {**CENTRAL, "crest": "unknown", "mopac": "unavailable"}
    assert program_mismatches(central, cell) == [
        "crest: unknown (central unknown)",
        "mopac: unavailable (central 23.2.2)",
    ]


def test_a_program_the_central_does_not_use_is_a_mismatch() -> None:
    assert program_mismatches(CENTRAL, {**CENTRAL, "orca": "6.0"}) == [
        "orca: 6.0 (central ausente)"
    ]
