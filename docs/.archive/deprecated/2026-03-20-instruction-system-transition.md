# Instruction System Transition Archive

Archived on 2026-03-20 during the instruction-system cleanup.

## What Changed

The repository previously had multiple competing sources of agent guidance:

- `AGENTS.md`
- `docs/CLAUDE.md`
- `CONTEXT.md`
- `docs/WORKFLOW.md`
- `.claude/README.md`
- `.claude/DEVELOPMENT_PROCEDURES.md`
- several specialized skill files with their own workflow assumptions

Those documents had drifted into conflicting authority and contained stale
references such as:

- fixed "Phase C / BATCH 12" status
- `tasks/*` files that do not exist in the repository
- slash commands and Serena/awesome-claude-code workflows not present here
- outdated dataset names and MOPAC `.aux` parsing guidance
- duplicated or inconsistent quality-gate commands

## New Design

The current instruction hierarchy is:

1. `/CLAUDE.md` as the single authoritative local behavior contract
2. `/AGENTS.md` as a short compatibility entry point
3. `.claude/*` supporting workflow references and specialized skills
4. `docs/.archive/` for retired historical context

## Notes

The old material remains recoverable through git history. This archive note
exists so the cleanup decision is visible from inside the repository without
leaving the current documentation tree.
