"""Tests for the rdkit-stubs patch helper."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType


def _load_patch_module() -> ModuleType:
    repo_root = Path(__file__).resolve().parents[2]
    script_path = repo_root / "scripts" / "patch_rdkit_stubs.py"
    spec = importlib.util.spec_from_file_location("patch_rdkit_stubs", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_patch_rdchem_stubs_fixes_invalid_getprop_overloads(tmp_path: Path) -> None:
    patch_mod = _load_patch_module()
    broken = patch_mod._BROKEN_SIGNATURE
    fixed = patch_mod._FIXED_SIGNATURE

    stub_root = tmp_path / "rdkit-stubs"
    target = stub_root / "Chem" / "rdchem.pyi"
    target.parent.mkdir(parents=True)
    target.write_text(
        "\n".join(
            [
                "@typing.overload",
                f"def GetProp(self, key: str, {broken} -> typing.Any:",
                "    ...",
            ]
        ),
        encoding="utf-8",
    )

    replacements, remaining = patch_mod.patch_rdchem_stubs(stub_root)

    assert replacements == 1
    assert remaining == 0
    assert fixed in target.read_text(encoding="utf-8")


def test_patch_rdchem_stubs_is_idempotent(tmp_path: Path) -> None:
    patch_mod = _load_patch_module()
    fixed = patch_mod._FIXED_SIGNATURE

    stub_root = tmp_path / "rdkit-stubs"
    target = stub_root / "Chem" / "rdchem.pyi"
    target.parent.mkdir(parents=True)
    target.write_text(
        f"def GetProp(self, key: str, {fixed} -> typing.Any:\n    ...\n",
        encoding="utf-8",
    )

    replacements, remaining = patch_mod.patch_rdchem_stubs(stub_root)

    assert replacements == 0
    assert remaining == 0
