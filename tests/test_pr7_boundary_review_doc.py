"""Guard test: the PR7 Package Boundary Review doc claimed by CHANGELOG exists.

CHANGELOG.md records the creation of ``docs/plans/pr7-package-boundary-review.md``
with a dependency matrix and a recorded location decision. This guard keeps that
claim honest — if the doc is removed, the claim must be removed too.
"""

from __future__ import annotations

from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_DOC = _REPO_ROOT / "docs" / "plans" / "pr7-package-boundary-review.md"


def test_boundary_review_doc_exists() -> None:
    assert _DOC.is_file(), f"Missing boundary review doc: {_DOC}"


def test_boundary_review_doc_records_decision() -> None:
    text = _DOC.read_text(encoding="utf-8")
    # The recorded decision: adopt src/grimperium/results/, not an external
    # packages/grimperium-results/ package.
    assert "src/grimperium/results/" in text
    assert "packages/grimperium-results/" in text
    assert "Adopt Option B" in text


def test_changelog_reference_is_backed_by_the_doc() -> None:
    changelog = (_REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    if "pr7-package-boundary-review.md" in changelog:
        assert _DOC.is_file()
