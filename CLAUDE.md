# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

This is the authoritative local instruction file for Grimperium. It adds
project-specific rules on top of the global agent behavior model. If another
local document disagrees with this file, follow this file.

## Operating Principles

- Read the relevant code, tests, and docs before making claims about behavior.
- For non-trivial work, make the scope explicit and plan before editing.
- Verify before declaring work done. Report exact commands run and anything not
  verified.
- Do not widen scope silently. Ask before changing architecture, dependencies,
  or unrelated modules.

## Development Commands

### Quality Gates (aligned with pyproject.toml)

```bash
poetry run black src/ tests/                         # Format
poetry run ruff check src/ tests/                    # Lint
poetry run mypy src/ --strict                        # Type check
poetry run pytest tests/ -v --cov=src/grimperium     # Full test suite with coverage
poetry run pre-commit run --all-files                # All gates at once
```

### Running Tests

```bash
poetry run pytest tests/ -v                          # All tests
poetry run pytest tests/unit/ -v                     # Unit tests only
poetry run pytest tests/ml/ -v                       # ML pipeline tests only
poetry run pytest tests/cli/ -v                      # CLI tests only
poetry run pytest tests/integration/ -v              # Integration tests
poetry run pytest tests/experiments/ -v              # Hypothesis/stress tests
poetry run pytest tests/path/test_file.py -v         # Single file
poetry run pytest tests/path/test_file.py::test_fn -v  # Single test
poetry run pytest -m "not slow" -v                   # Skip slow-marked tests
poetry run pytest -m integration -v                  # Only integration-marked tests
```

### Execution Rule

- Run targeted checks while iterating.
- Before handing off substantial code changes, run the broadest relevant gate
  you can afford.
- If a gate is skipped, blocked, or only partially run, say so explicitly.

## Architecture Overview

### Delta-Learning Pipeline (the core idea)

Grimperium corrects fast-but-inaccurate quantum chemistry (PM7) using ML:

```
delta = H298_CBS (expensive, accurate) − H298_PM7 (cheap, ~80% correct)
prediction = H298_PM7 + ensemble.predict(features)   ← the correction
```

The model learns the *error* (delta), not the absolute value. This converges
with less training data.

### Package Interaction

```
CLI (cli/)  →  API (api.py)  →  Core (core/)
                                   ├── DeltaLearner (delta_learning.py)
                                   ├── BatchOrchestrator
                                   └── Molecule model
                                        ↓
                               ML Pipeline (ml/)
                               data_loader → features → trainer → evaluator
                                             → persistence → predictor → charts
                                        ↓
                            CREST+MOPAC Pipeline (crest_pm7/)
                            pipeline.py → conformer_generator → conformer_selector
                                        → mopac_optimizer → energy_extractor
                                        → mopac_descriptors
                                        ↓
                               Data Layer (data/)
                               loader.py, fusion.py, semiempirical adapters
```

### Data Leakage Prevention (critical design constraint)

The ML pipeline enforces strict ordering to prevent leakage:
1. CSV filtering (status/quality) **before** train/test split
2. `FeaturePipeline` imputer fitted on **train only** (`features.py`)
3. Test set transformed using train statistics
4. `evaluator.py` requires a **pre-fitted** pipeline (never re-fits on eval data)
5. `persistence.py` bundles the fitted pipeline with the model

### Ensemble Architecture

`DeltaLearningEnsemble` (models/delta_ensemble.py):
- KRR (Kernel Ridge Regression) — captures smooth molecular similarity
- XGBoost — captures non-linear interactions
- Weighted average prediction (default 50/50, configurable)

### ML Quality Gate (ml/gate.py)

Auto-rejects models that fail: MAE ≤ 3.5 kcal/mol, R² ≥ 0.97, RMSE ≤ 5.0 kcal/mol.

### Status Tracking (CSV columns)

| Column             | Values                                    |
|--------------------|-------------------------------------------|
| `status`           | OK, PENDING, RERUN, SKIP, FAILED          |
| `cbs_quality_flag` | OK, SUSPECT                               |
| `crest_status`     | SUCCESS, TIMEOUT, ERROR, NOT_ATTEMPTED    |
| `mopac_status`     | SUCCESS, TIMEOUT, ERROR, NOT_ATTEMPTED    |

ML training filters to `status == 'OK' AND cbs_quality_flag == 'OK'`.

### Configuration (config.py)

Three nested dataclasses: `GrimperiumConfig` → `ModelConfig` + `FeatureConfig` +
`DataConfig`. All defaults are in Python dataclass fields. No YAML/JSON loader yet.

### Feature Vector (270 dimensions)

- Tabular (3): nheavy, charge, multiplicity
- Morgan Fingerprints (256): ECFP-like circular substructure patterns
- RDKit descriptors (11): MolWt, TPSA, LogP, rotatable bonds, etc.

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

Test layout:

- `tests/unit/`, `tests/integration/`, `tests/experiments/`, `tests/cli/`,
  `tests/ml/`
- Shared fixtures in `tests/conftest.py` and package-specific `conftest.py` files.
- Markers: `@pytest.mark.slow`, `@pytest.mark.integration`, `@pytest.mark.unit`

Do not describe the project as locked to an old fixed phase or numbered batch
unless the current task explicitly provides that context.

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
