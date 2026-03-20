# AGENTS.md - Grimperium

This file is a short compatibility entry point for agents that look for
`AGENTS.md` first. The authoritative project-specific guidance now lives in
[`CLAUDE.md`](/home/igor/Projetos/grimperium/CLAUDE.md). If this file and
another local guide differ, follow `CLAUDE.md`.

## Quick Orientation

- Python 3.10+ scientific/ML repository with CLI, CREST+MOPAC pipeline,
  data-fusion, and ML training/prediction code.
- Active package layout:
  `src/grimperium/{cli,crest_pm7,core,data,ml,models,utils}`.
- Active tests:
  `tests/{unit,integration,experiments,cli,ml}` plus shared fixtures.

## Standard Commands

- `black src/ tests/`
- `ruff check src/ tests/`
- `mypy src/ --strict`
- `pytest tests/ -v --cov=src/grimperium`
- `pre-commit run --all-files`

## Local Doc Map

- Main authority: [`CLAUDE.md`](/home/igor/Projetos/grimperium/CLAUDE.md)
- Workflow reference:
  [`.claude/DEVELOPMENT_PROCEDURES.md`](/home/igor/Projetos/grimperium/.claude/DEVELOPMENT_PROCEDURES.md)
- Changelog template:
  [`.claude/CHANGELOG_TEMPLATE.md`](/home/igor/Projetos/grimperium/.claude/CHANGELOG_TEMPLATE.md)
- Historical materials:
  [`docs/.archive/`](/home/igor/Projetos/grimperium/docs/.archive/README.md)
