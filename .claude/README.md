# `.claude/` Directory Guide

This directory contains supporting material for local agent workflows.
[`/CLAUDE.md`](/home/igor/Projetos/grimperium/CLAUDE.md) is the only
authoritative local behavior contract. Nothing inside `.claude/` should act as
an independent constitution.

## What Lives Here

- `DEVELOPMENT_PROCEDURES.md`: supporting execution and verification checklist.
- `CHANGELOG_TEMPLATE.md`: template for `[Unreleased]` changelog entries.
- `agents/`: specialized review personas.
- `skills/`: bounded task-specific playbooks.
- `hooks/`: optional helper scripts. Useful, but not a substitute for explicit
  verification.
- `settings*.json`: local tool configuration.

## How To Use This Directory

1. Start with `CLAUDE.md`.
2. Use `DEVELOPMENT_PROCEDURES.md` for practical workflow details.
3. Use a skill or subagent spec only when the task clearly matches its scope.
4. If a support file conflicts with `CLAUDE.md`, fix the support file instead of
   treating the conflict as two valid authorities.

## Current Skills And Agents

- `skills/grimperium-ci.md`: reproduce and triage CI failures with repo-native
  commands.
- `skills/grimperium-tests.md`: test execution strategy and command shortcuts.
- `skills/grimperium-format.md`: formatting and linting workflow.
- `skills/grimperium-docs.md`: documentation maintenance and changelog updates.
- `skills/grimperium-csv.md`: current CSV schema reference and migration guardrails.
- `skills/grimperium-mopac.md`: current MOPAC `.out` parsing reference.
- `agents/scientific-reviewer.md`: chemistry- and numerics-focused review lens.

## Maintenance Rule

When a command, path, schema, or workflow changes:

1. Update `CLAUDE.md` if the main contract changed.
2. Update the affected support file or skill.
3. Archive or redirect stale guidance instead of leaving competing instructions
   in place.
