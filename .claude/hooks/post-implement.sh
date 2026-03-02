#!/usr/bin/env bash
# post-implement.sh — Quality gate after implementation.
# Run manually: bash .claude/hooks/post-implement.sh
#
# Exit codes:
#   0 = all checks pass (or only warnings)
#   1 = ruff violations found (blocking)
#
# Ruff is the only blocker — black/mypy/pytest emit warnings but don't block.

set -uo pipefail

PASS="✓"
WARN="⚠"
FAIL="✗"
EXIT_CODE=0

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Grimperium — Quality Gate"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# 1. Ruff — BLOCKING
echo ""
echo "  [1/4] Ruff lint..."
if ruff check src/ tests/ --quiet 2>&1; then
    echo "  ${PASS} Ruff: no violations"
else
    echo "  ${FAIL} Ruff: violations found (blocking)"
    EXIT_CODE=1
fi

# 2. Black — WARNING only
echo ""
echo "  [2/4] Black format check..."
if black --check src/ tests/ --quiet 2>&1; then
    echo "  ${PASS} Black: formatting OK"
else
    echo "  ${WARN} Black: formatting issues (run 'black src/ tests/' to fix)"
fi

# 3. Mypy — WARNING only
echo ""
echo "  [3/4] Mypy type check..."
if mypy src/grimperium --ignore-missing-imports --quiet 2>&1; then
    echo "  ${PASS} Mypy: no type errors"
else
    echo "  ${WARN} Mypy: type issues found (review above)"
fi

# 4. Pytest unit tests — WARNING only
echo ""
echo "  [4/4] Pytest unit tests..."
if pytest tests/unit/ -q --tb=no 2>&1; then
    echo "  ${PASS} Pytest: all unit tests pass"
else
    echo "  ${WARN} Pytest: test failures found (review above)"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "  Result: PASS"
else
    echo "  Result: FAIL (fix ruff errors before committing)"
fi
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

exit ${EXIT_CODE}
