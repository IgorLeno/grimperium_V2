# Phase D — ML Pipeline (TDD)

## Context

Create `src/grimperium/ml/` module with data loading, feature engineering, training, and evaluation for the delta-learning pipeline. TDD approach: tests first, then implementation.

## Key Dependencies (existing code)

- `DeltaLearner`: `src/grimperium/core/delta_learning.py` — fit(X, y_cbs, y_pm7), predict(X, y_pm7), evaluate(X, y_cbs, y_pm7)
- `CSVDataLoader`: `src/grimperium/core/csv_data_loader.py` — load_dataframe() with strict/permissive modes
- `compute_all_metrics`: returns dict with keys rmse, mae, r2, mape, max_error
- `MoleculeStatus`: `src/grimperium/crest_pm7/batch/enums.py` — valid status values (OK, Pending, etc.)
- `MatrixFloat`: `NDArray[np.floating[Any]]` from `src/grimperium/__init__.py`

## Implementation Steps

### Step 1: Create test files (RED phase)

- `tests/ml/__init__.py`
- `tests/ml/conftest.py` — fixtures with 10-row synthetic CSV
- `tests/ml/test_data_loader.py` — 5 tests
- `tests/ml/test_features.py` — 4 tests
- `tests/ml/test_trainer.py` — 3 tests
- `tests/ml/test_evaluator.py` — 2 tests

### Step 2: Confirm all tests fail (ModuleNotFoundError)

### Step 3: Implement modules (GREEN phase)

- `src/grimperium/ml/__init__.py`
- `src/grimperium/ml/data_loader.py`
- `src/grimperium/ml/features.py`
- `src/grimperium/ml/trainer.py`
- `src/grimperium/ml/evaluator.py`

### Step 4: All tests pass

### Step 5: Quality gates (mypy, ruff, black)

## Model selection & experimentation

- **Cross-validation strategy:** Use K-fold (default k=5) when samples are i.i.d.; switch to grouped CV when multiple rows belong to the same molecule family/source split to avoid fold leakage. Integrate CV orchestration around `trainer.train(...)` inputs and call `DeltaLearner.fit(...)` inside each fold.
- **Hyperparameter tuning:** Start with constrained grid/random search for fast baselines, then optionally Optuna for broader exploration. Evaluate each trial with validation folds (not test split), aggregating fold metrics before selecting candidate parameters.
- **Experiment tracking/versioning:** Prefer MLflow when available (params, fold metrics, artifacts). If MLflow is not available, enforce deterministic artifact naming (`model_<date>_<hash>_<cvscore>.pkl`) plus JSON metadata containing `feature_cols`, split seed, fold metrics, and selected params.
- **Model selection criteria beyond a single metric:** Define metric hierarchy (`rmse` primary, `mae` secondary, `r2` tertiary), then apply tie-breakers on calibration/stability across folds, robustness under distribution shift slices, and inference cost/latency.
- **Pipeline integration points:** Keep train/test holdout logic in `trainer.train(...)`; use `DeltaLearner.evaluate(...)` for fold/holdout reporting; run final full-dataset reporting through evaluator only with a pre-fitted `FeaturePipeline` from training to preserve leakage guarantees.

## Important Notes

- CSVDataLoader validates rows against MoleculeStatus enum — test CSV must use valid enum values
- CSVDataLoader requires: mol_id, smiles, nheavy, status columns
- 3 features (mopac_homo/lumo/gap_ev) have 12.4% NaN — SimpleImputer median
- H298_cbs already in thermo_pm7.csv — no external merge needed
- trainer.py spec in user prompt has truncated FeaturePipeline/train code — will reconstruct from context

## CRITICAL: Data Leakage Prevention in trainer.py

The split MUST happen on the DataFrame BEFORE fit_transform:

1. `train_test_split(df, y_cbs, y_pm7)` → df_train, df_test, y_cbs_train, y_cbs_test, y_pm7_train, y_pm7_test
2. `pipeline.fit_transform(df_train)` → X_train (imputer fitted on train only)
3. `pipeline.transform(df_test)` → X_test (imputer transforms test, no leakage)
