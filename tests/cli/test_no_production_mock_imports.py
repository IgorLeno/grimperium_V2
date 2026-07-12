"""Boundary: production code under src/grimperium must not import cli.mock_data."""

from __future__ import annotations

import ast
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[2] / "src" / "grimperium"
ALLOWED = {
    # mock_data itself may exist; nothing else in production may import it.
}


def _imports_mock_data(path: Path) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if mod == "grimperium.cli.mock_data" or mod.endswith(".mock_data"):
                return True
            if mod == "grimperium.cli" and any(
                alias.name == "mock_data" for alias in node.names
            ):
                return True
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "grimperium.cli.mock_data":
                    return True
    return False


def test_no_production_imports_of_mock_data() -> None:
    offenders: list[str] = []
    for path in SRC_ROOT.rglob("*.py"):
        rel = path.relative_to(SRC_ROOT)
        if rel.as_posix() == "cli/mock_data.py":
            continue
        if path in ALLOWED:
            continue
        if _imports_mock_data(path):
            offenders.append(str(rel))
    assert (
        offenders == []
    ), "Production modules must not import grimperium.cli.mock_data: " + ", ".join(
        offenders
    )


def test_database_registry_does_not_embed_legacy_default_aliases() -> None:
    registry_path = SRC_ROOT / "cli" / "database_registry.py"
    source = registry_path.read_text(encoding="utf-8")

    assert "_create_defaults" not in source
    assert '"CBS"' not in source
    assert '"PM7"' not in source
    assert '"NIST"' not in source
