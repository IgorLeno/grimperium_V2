from __future__ import annotations

import ast
from pathlib import Path


def test_delta_ensemble_does_not_import_type_aliases_from_package_root() -> None:
    """Avoid package-root re-entry while grimperium.__init__ is initializing."""
    source_path = (
        Path(__file__).parents[2] / "src" / "grimperium" / "ml" / "delta_ensemble.py"
    )
    tree = ast.parse(source_path.read_text(encoding="utf-8"))

    root_alias_imports = [
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module == "grimperium"
        for alias in node.names
        if alias.name in {"DictStrAny", "MatrixFloat"}
    ]

    assert root_alias_imports == []
