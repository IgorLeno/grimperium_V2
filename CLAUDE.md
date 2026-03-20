# Grimperium Agent Guide

This file is the authoritative local instruction file for Grimperium. It adds
project-specific rules on top of the global agent behavior model. If another
local document disagrees with this file, follow this file.

## Operating Principles

- Read the relevant code, tests, and docs before making claims about behavior.
- For non-trivial work, make the scope explicit and plan before editing.
- Verify before declaring work done. Report exact commands run and anything not
  verified.
- Do not widen scope silently. Ask before changing architecture, dependencies,
  or unrelated modules.

## Repo Snapshot

Current top-level package layout:

- `src/grimperium/cli/`: Rich/questionary CLI, controllers, views, settings.
- `src/grimperium/crest_pm7/`: CREST + MOPAC PM7 pipeline, CSV backend,
  monitoring, parsing, validation.
- `src/grimperium/core/`: batch orchestration, delta-learning, metrics,
  molecule/value utilities.
- `src/grimperium/data/`: dataset loading, fusion, semiempirical adapters.
- `src/grimperium/ml/`: feature pipelines, training, evaluation, prediction,
  persistence, charts.
- `src/grimperium/models/`: model wrappers and ensembles.
- `src/grimperium/utils/`: validation, logging, feature helpers.

Current test layout:

- `tests/unit/`
- `tests/integration/`
- `tests/experiments/`
- `tests/cli/`
- `tests/ml/`
- shared fixtures in `tests/conftest.py` and package-specific `conftest.py`
  files.

Do not describe the project as locked to an old fixed phase or numbered batch
unless the current task explicitly provides that context.

## Quality Gates

Use repo-native commands and keep them aligned with `pyproject.toml`.

Primary commands:

- `black src/ tests/`
- `ruff check src/ tests/`
- `mypy src/ --strict`
- `pytest tests/ -v --cov=src/grimperium`
- `pre-commit run --all-files`

Execution rule:

- Run targeted checks while iterating.
- Before handing off substantial code changes, run the broadest relevant gate
  you can afford.
- If a gate is skipped, blocked, or only partially run, say so explicitly.

## Project-Specific Rules

- CSV schema authority lives in
  `src/grimperium/crest_pm7/batch/csv_manager.py` and `tests/test_csv_schema.py`.
  Do not trust old docs, exported CSVs, or historical notes over the current
  schema code and tests.
- MOPAC descriptor extraction currently reads `.out` files via
  `src/grimperium/crest_pm7/mopac_descriptors.py`. Do not reintroduce `.aux`
  parsing assumptions without checking the current code and tests.
- Active dataset names are `data/thermo_cbs_chon.csv` and `data/thermo_pm7.csv`.
  Removed dataset names can appear in historical docs and changelog entries, but
  they are not current workflow authority.
- CLI changes must preserve `BaseView` navigation, visible Rich feedback for
  state changes, and graceful cancellation paths.
- Scientific or data-pipeline changes need grounded verification across adjacent
  loaders, fusion, ML, and pipeline tests. Do not validate only the edited
  module if the behavior crosses package boundaries.

## Batch And Changelog Habit

- Use "batch" as a scoped unit of work when a task is non-trivial, but do not
  invent fixed historical batch numbers.
- Keep `CHANGELOG.md` updated under `[Unreleased]` for completed batches or
  materially user-visible or workflow-relevant changes.
- Use `.claude/CHANGELOG_TEMPLATE.md` for entry shape.
- When commands, paths, or workflow expectations change, update the nearby docs
  in the same batch.

## Local Authority Order

1. `CLAUDE.md`
2. `AGENTS.md`
3. `.claude/DEVELOPMENT_PROCEDURES.md`
4. `.claude/README.md`, `CONTEXT.md`, `docs/WORKFLOW.md`, `docs/TESTING_GUIDE.md`
5. `.claude/skills/*.md` and `.claude/agents/*.md`
6. `docs/.archive/`
