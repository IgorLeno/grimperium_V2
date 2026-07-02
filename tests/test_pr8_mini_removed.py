from __future__ import annotations

from pathlib import Path

import tomllib

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_pr8_removes_mini_trees() -> None:
    mini_paths = [
        REPO_ROOT / "packages" / "grimperium-mini",
        REPO_ROOT / "src" / "grimperium_mini",
        REPO_ROOT / "tests" / "grimperium_mini",
    ]

    assert [
        str(path.relative_to(REPO_ROOT)) for path in mini_paths if path.exists()
    ] == []


def test_pr8_removes_mini_poetry_dependency() -> None:
    pyproject = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
    dev_dependencies = pyproject["tool"]["poetry"]["group"]["dev"]["dependencies"]

    assert "grimperium-mini" not in dev_dependencies


def test_pr8_removes_mini_ci_wiring() -> None:
    active_files = [
        REPO_ROOT / ".github" / "workflows" / "ci.yml",
        REPO_ROOT / ".gitignore",
        REPO_ROOT / ".pre-commit-config.yaml",
        REPO_ROOT / "poetry.lock",
    ]
    forbidden_fragments = [
        "mini-test",
        "mini-lint",
        "mini-typecheck",
        "grimperium-mini",
        "grimperium_mini",
        "packages/grimperium-mini",
    ]

    offenders = {
        str(path.relative_to(REPO_ROOT)): [
            fragment for fragment in forbidden_fragments if fragment in path.read_text()
        ]
        for path in active_files
    }

    assert {path: matches for path, matches in offenders.items() if matches} == {}
