# Grimperium Development Procedures

This document is a supporting workflow reference. The authoritative local
contract lives in
[`/CLAUDE.md`](/home/igor/Projetos/grimperium/CLAUDE.md).

## Standard Work Loop

### 1. Scope And Plan

- Clarify the requested outcome and affected modules.
- For non-trivial work, write down a short plan before editing.
- Use existing code, tests, and docs as source of truth. Do not carry forward
  stale batch or phase assumptions from historical notes.

### 2. Inspect Before Editing

- Read the relevant source files first.
- Read the closest tests before changing behavior.
- For pipeline/schema work, check the corresponding integration points, not just
  the file being edited.

### 3. Implement In Small Verified Steps

- Prefer the smallest coherent change that resolves the issue.
- Add or update tests when behavior changes.
- Keep edits local unless the task explicitly requires broader refactoring.

### 4. Validate

Repo-native commands:

```bash
black src/ tests/
ruff check src/ tests/
mypy src/ --strict
pytest tests/ -v --cov=src/grimperium
```

Verification policy:

- Run the narrowest relevant check while iterating.
- Before handoff, run the broadest relevant gate you can.
- State exactly what you ran and what you did not run.

### 5. Update Docs And Changelog

- Update `CHANGELOG.md` under `[Unreleased]` for completed batches or meaningful
  workflow-visible changes.
- Use `.claude/CHANGELOG_TEMPLATE.md` to keep entries consistent.
- Update nearby docs when commands, file paths, schemas, or workflows changed.
- Do not leave historical notes as current authority. Redirect or archive them.

### 6. Commit Hygiene

- Keep commit messages descriptive and scoped to the actual change.
- Reference a batch label only when one was actually provided or used for the
  work.
- Do not claim a full verification sweep if you only ran targeted checks.

## Change Types That Need Extra Care

- CLI navigation and settings flows: verify the affected view path and
  cancellation behavior.
- CSV schema changes: confirm `BatchCSVManager.get_schema()` and the schema
  tests still agree.
- MOPAC parsing changes: validate against the current `.out` parser tests and
  integration tests.
- Data and ML changes: validate loaders, fusion, and downstream training or
  prediction paths as needed.

## Changelog Checklist

- Entry is under `[Unreleased]`
- Section matches the change type (`Added`, `Changed`, `Fixed`, `Deprecated`,
  `Removed`, `Security`)
- Date uses absolute format: `(YYYY-MM-DD)`
- Files or user-visible areas are named explicitly
- Verification is mentioned when it adds value
- No invented version number or future release note is added
