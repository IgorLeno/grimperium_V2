"""Shared helpers for the ``scripts/conformer_validation`` tests.

The experiment lives under ``scripts/`` (not an installed package), so the
tests load it by path, the same way ``test_patch_rdkit_stubs.py`` loads
its script.
"""

from __future__ import annotations

import importlib
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

from semi_imperium.conformers import ConformerGeometry
from semi_imperium.mopac import DisplacementLineage, ForceRun, OptimizationRun

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = REPO_ROOT / "scripts" / "conformer_validation"
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "crest_atom_order"


def load(name: str) -> ModuleType:
    """Import ``conformer_validation.<name>`` from ``scripts/`` by path."""
    if "conformer_validation" not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            "conformer_validation",
            SCRIPT_DIR / "__init__.py",
            submodule_search_locations=[str(SCRIPT_DIR)],
        )
        assert spec is not None and spec.loader is not None
        package = importlib.util.module_from_spec(spec)
        sys.modules["conformer_validation"] = package
        spec.loader.exec_module(package)
    return importlib.import_module(f"conformer_validation.{name}")


def force_output(*frequencies: float) -> str:
    """Minimal MOPAC FORCE output the classifier accepts."""
    vibrations = "\n".join(
        f"VIBRATION {mode} A1 FREQ. {frequency:.3f}"
        for mode, frequency in enumerate(frequencies, start=1)
    )
    return (
        "START OF FORCE CALCULATION OUTPUT\n"
        "GRADIENTS WERE INITIALLY ACCEPTABLY SMALL\n"
        "GRADIENT NORM = 0.041\n"
        "DESCRIPTION OF VIBRATIONS\n"
        f"{vibrations}\n"
        "== MOPAC DONE ==\n"
    )


class FakeMopac:
    """MOPAC double: geometry unchanged, energy from the conformer index."""

    def __init__(self, *, fail: bool = False) -> None:
        self.fail = fail
        self.optimized: list[int] = []

    def optimize(
        self,
        *,
        hamiltonian: str,
        geometry: ConformerGeometry,
        source_conformer_index: int,
        attempt_id: str,
        displacement: DisplacementLineage | None,
    ) -> OptimizationRun:
        del hamiltonian, attempt_id, displacement
        if self.fail:
            raise RuntimeError("fake MOPAC exploded")
        self.optimized.append(source_conformer_index)
        return OptimizationRun(
            converged=True,
            geometry=geometry,
            heat_of_formation_kcal_mol=-100.0 + 0.1 * source_conformer_index,
        )

    def verify_force(
        self,
        *,
        hamiltonian: str,
        optimization: OptimizationRun,
        attempt_id: str,
    ) -> ForceRun:
        del hamiltonian, optimization, attempt_id
        return ForceRun(output=force_output(-2.5, 56.0, 1240.0))
