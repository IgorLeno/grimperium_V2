#!/usr/bin/env python3
"""Patch invalid GetProp overload signatures in bundled rdkit-stubs.

RDKit 2026.3.3 ships ``rdkit-stubs/Chem/rdchem.pyi`` overloads where
``default`` lacks a default value after ``autoConvert=False``. Mypy treats
that as a syntax error and aborts before checking project code.

Upstream: https://github.com/rdkit/rdkit/issues/9335
"""

from __future__ import annotations

import sys
from pathlib import Path

try:
    import rdkit
except ImportError:
    rdkit = None  # type: ignore[assignment,misc]

_BROKEN_SIGNATURE = "autoConvert: bool = False, default: typing.Any)"
_FIXED_SIGNATURE = "autoConvert: bool = False, default: typing.Any = ...)"
_TARGET_REL = Path("Chem") / "rdchem.pyi"


def find_rdkit_stubs_root() -> Path | None:
    """Return the installed ``rdkit-stubs`` directory, if present."""
    if rdkit is None:
        return None

    candidate = Path(rdkit.__file__).resolve().parent.parent / "rdkit-stubs"
    if candidate.is_dir():
        return candidate
    return None


def patch_rdchem_stubs(stub_root: Path) -> tuple[int, int]:
    """Patch ``rdchem.pyi`` in-place. Returns (replacements, remaining_broken)."""
    target = stub_root / _TARGET_REL
    if not target.is_file():
        return 0, 0

    original = target.read_text(encoding="utf-8")
    patched = original.replace(_BROKEN_SIGNATURE, _FIXED_SIGNATURE)
    if patched == original:
        remaining = patched.count(_BROKEN_SIGNATURE)
        return 0, remaining

    target.write_text(patched, encoding="utf-8")
    remaining = patched.count(_BROKEN_SIGNATURE)
    replacements = original.count(_BROKEN_SIGNATURE) - remaining
    return replacements, remaining


def main() -> int:
    stub_root = find_rdkit_stubs_root()
    if stub_root is None:
        print(
            "patch_rdkit_stubs: rdkit not installed; nothing to patch", file=sys.stderr
        )
        return 0

    replacements, remaining = patch_rdchem_stubs(stub_root)
    target = stub_root / _TARGET_REL
    if remaining:
        print(
            f"patch_rdkit_stubs: {remaining} invalid signature(s) remain in {target}",
            file=sys.stderr,
        )
        return 1

    if replacements:
        print(
            f"patch_rdkit_stubs: patched {replacements} GetProp overload(s) in {target}"
        )
    else:
        print(f"patch_rdkit_stubs: {target} already patched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
