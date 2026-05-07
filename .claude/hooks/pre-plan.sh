#!/usr/bin/env bash
# pre-plan.sh — Informational context check at session start.
# Always exits 0 (non-blocking). Prints environment summary.

set -euo pipefail

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Grimperium — Session Context"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# 1. Confirm we're at the project root
if [[ -f "pyproject.toml" ]]; then
    echo "  ✓ Project root detected (pyproject.toml found)"
else
    echo "  ⚠ pyproject.toml not found — check working directory"
fi

# 2. Virtual environment status
if [[ -n "${VIRTUAL_ENV:-}" ]]; then
    echo "  ✓ venv active: ${VIRTUAL_ENV}"
else
    echo "  ⚠ No virtual environment active"
fi

# 3. Git context
if git rev-parse --git-dir > /dev/null 2>&1; then
    BRANCH=$(git branch --show-current 2>/dev/null || echo "unknown")
    MODIFIED=$(git diff --name-only 2>/dev/null | wc -l | tr -d ' ')
    STAGED=$(git diff --cached --name-only 2>/dev/null | wc -l | tr -d ' ')
    echo "  ✓ Branch: ${BRANCH}"
    echo "    Modified files: ${MODIFIED}  |  Staged: ${STAGED}"
else
    echo "  ⚠ Not inside a git repository"
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

exit 0
