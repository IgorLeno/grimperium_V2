#!/usr/bin/env bash
# pre-commit-reminder.sh — Remind about CHANGELOG.md before committing.
# Install as a git hook: cp .claude/hooks/pre-commit-reminder.sh .git/hooks/pre-commit
#
# Not connected to Claude Code hooks (would trigger on every Bash tool call).
# Use manually or as a git pre-commit hook.

set -euo pipefail

STAGED_FILES=$(git diff --cached --name-only 2>/dev/null || true)

if echo "${STAGED_FILES}" | grep -q "CHANGELOG.md"; then
    # Changelog included — all good
    exit 0
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  ⚠ CHANGELOG.md not staged"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  Staged files:"
echo "${STAGED_FILES}" | sed 's/^/    /'
echo ""
echo "  Did you forget to update CHANGELOG.md?"
echo "  Options:"
echo "    1. Stage it:  git add CHANGELOG.md && git commit"
echo "    2. Skip:      SKIP_CHANGELOG=1 git commit"
echo ""

# Allow bypass via env variable
if [[ "${SKIP_CHANGELOG:-0}" == "1" ]]; then
    echo "  Skipping changelog check (SKIP_CHANGELOG=1)"
    exit 0
fi

# Interactive: wait for user decision
if [[ -t 1 ]]; then
    read -rp "  Continue without CHANGELOG.md? [y/N] " answer
    case "${answer}" in
        [yY]) exit 0 ;;
        *) echo "  Aborted. Stage CHANGELOG.md and retry."; exit 1 ;;
    esac
fi

# Non-interactive (e.g., CI): block by default
echo "  Non-interactive mode: blocking commit."
exit 1
