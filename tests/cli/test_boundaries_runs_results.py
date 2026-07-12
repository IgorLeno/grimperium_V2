from __future__ import annotations

import ast
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[2] / "src" / "grimperium"


def _imports_for(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imports.add(node.module)
    return imports


def test_runs_and_results_do_not_import_cli_or_views() -> None:
    offenders: list[str] = []
    for package in ("runs", "results"):
        for path in (SRC_ROOT / package).rglob("*.py"):
            imports = _imports_for(path)
            if any(
                name == "grimperium.cli"
                or name.startswith("grimperium.cli.")
                or ".views" in name
                for name in imports
            ):
                offenders.append(str(path.relative_to(SRC_ROOT)))

    assert offenders == []


def test_results_view_delegates_scientific_analysis_to_service() -> None:
    path = SRC_ROOT / "cli" / "views" / "results_view.py"
    source = path.read_text(encoding="utf-8")

    assert "_compute_divergence_stats" not in source
    assert "relative_error_pct" not in source
    assert "np." not in source
    assert "analyze_dataset" in source or "ResultsService" in source
