# Grimperium v0.2.0 - Data Layer Implementation Plan

## Overview

Implement the real data layer for Grimperium, evolving existing stubs in `ChemperiumLoader` and `DataFusion` into functional implementations, while creating comprehensive test fixtures.

**Design Decisions (confirmed with user):**
- **Loader**: Evolve existing interface (keep `load(path)`) but add 3-way split and filters
- **Fusion**: Extend existing class - add `select_task_view()` alongside merge methods
- **Columns**: Keep `H298_b3` and `S298` as OPTIONAL (only core columns required)

## Critical Files to Modify

| File | Action | Purpose |
|------|--------|---------|
| `src/grimperium/data/loader.py` | Evolve | Implement methods, keep interface |
| `src/grimperium/data/fusion.py` | Extend | Add task views to existing class |
| `tests/fixtures/mock_data.py` | Enhance | Add realistic fixtures |
| `tests/unit/test_loader.py` | Enhance | Implement skipped tests |
| `tests/unit/test_fusion.py` | Enhance | Implement skipped tests |
| `tests/integration/test_pipeline.py` | Enhance | Add data pipeline tests |
| `docs/architecture.md` | Update | Document data layer |

---

## Phase 1: ChemperiumLoader Implementation

### File: `src/grimperium/data/loader.py`

**Design:** Evolve existing interface - keep constructor and `load(path)` signature.

**Keep existing:**
- `__init__(self, validate: bool = True)`
- `REQUIRED_COLUMNS` = `["smiles", "charge", "multiplicity", "nheavy", "H298_cbs"]`
- `OPTIONAL_COLUMNS` = `["xyz", "H298_b3", "S298", "A", "B"]`
- `CP_COLUMNS` = `[f"cp_{i}" for i in range(1, 46)]`

**Methods to implement:**

1. **`load(path, columns=None, max_nheavy=None, allowed_elements=None)`**
   - Read CSV/Parquet by suffix detection
   - Validate required columns if `self.validate=True`
   - Apply optional filters (max_nheavy, allowed_elements)
   - Store in `self.data`

2. **`split(df=None, test_size=0.2, random_state=42, stratify_by=None)`**
   - Returns: `(train_df, test_df)` - existing 2-way split

3. **`train_val_test_split(df=None, test_size=0.2, val_size=0.1, random_state=42)` (NEW)**
   - Returns: `(train_df, val_df, test_df)` - 3-way split
   - Proportions: train=70%, val=10%, test=20% (default)

4. **`get_features(df=None, include_cp=False)`**
   - Returns tabular feature columns (nheavy, charge, multiplicity)
   - Optionally includes cp_1...cp_45

5. **`get_targets(df=None, target="H298_cbs")`**
   - Returns target column as Series

6. **`_validate_dataframe(df)` (private)**
   - Check REQUIRED_COLUMNS exist
   - Raise ValueError with clear message if missing

7. **`_apply_filters(df, max_nheavy, allowed_elements)` (private)**
   - Filter by `nheavy <= max_nheavy` if provided
   - Element filter is a stub for now (future RDKit integration)

---

## Phase 2: DataFusion Implementation

### File: `src/grimperium/data/fusion.py`

**Design:** Extend existing class - add task views alongside merge methods.

**Keep existing:**
- `__init__(target_column, semiempirical_column, delta_column)`
- `merged_data` attribute

**Methods to implement:**

1. **`merge(chemperium_df, semiempirical_df, on="smiles", how="inner", validate_merge=True)`**
   - Merge DataFrames on specified column
   - Validate merge if requested (check row counts)
   - Store in `self.merged_data`

2. **`compute_deltas(df=None)`**
   - Add delta column: `delta = target - semiempirical`
   - Works on merged_data if df not provided

3. **`select_task_view(df, task="enthalpy")` (NEW)**
   - `"enthalpy"`: X=features, y=H298_cbs
   - `"entropy"`: X=features, y=S298
   - `"heat_capacity"`: X=features, Y=cp_1...cp_45 (multi-output)
   - Returns: `(X: DataFrame, y: Series|DataFrame)`

4. **`get_training_data(df=None)`**
   - Returns `(features_df, delta_array)` for ML training

5. **`analyze_deltas(df=None)`**
   - Returns dict with stats: mean, std, min, max, median

6. **`_validate_merge(merged, expected_rows)` (private)**

7. **`_default_feature_columns(df, exclude)` (private helper)**
   - Returns `["nheavy", "charge", "multiplicity"]` minus excluded columns

---

## Phase 3: Mock Data Fixtures

### File: `tests/fixtures/mock_data.py`

**Enhance existing generators + add new functions:**

1. **Keep existing:** `generate_chemperium_mock()`, `generate_pm7_mock()`, `generate_features_mock()`, `generate_training_data()`

2. **Add `make_small_chemperium_df(n=100, random_state=42) -> DataFrame`**
   - Wrapper with sensible defaults for testing
   - Uses `generate_chemperium_mock()` internally
   - Ensures all optional columns present (H298_b3, S298, cp_1...cp_45)

3. **Add `make_chemperium_with_pm7_df(n=100, random_state=42) -> DataFrame`**
   - Extends `make_small_chemperium_df`
   - Adds `H298_pm7` column with systematic error (~3±1.5 kcal/mol from CBS)
   - Used to test delta computation

---

## Phase 4: Unit Tests

### File: `tests/unit/test_loader.py`

**Unskip and implement existing tests:**

1. `test_load_csv()` - Load from mock CSV (using tmp_path fixture)
2. `test_load_missing_file()` - Expect FileNotFoundError
3. `test_split()` - Verify 2-way split proportions
4. `test_get_features()` - Returns feature columns
5. `test_get_targets()` - Returns target Series
6. `test_validate_missing_columns()` - Expect ValueError
7. `test_validate_valid_dataframe()` - Should not raise

**Add new tests:**

8. `test_load_parquet()` - Load from mock Parquet file
9. `test_load_with_max_nheavy_filter()` - Verify filtering works
10. `test_train_val_test_split()` - Verify 3-way split proportions
11. `test_split_reproducibility()` - Same seed = same split

### File: `tests/unit/test_fusion.py`

**Unskip and implement existing tests:**

1. `test_merge_inner()` - Merge two DataFrames
2. `test_compute_deltas()` - Adds delta column
3. `test_analyze_deltas()` - Returns stats dict

**Add new tests:**

4. `test_select_task_view_enthalpy()` - X, y for H298_cbs
5. `test_select_task_view_entropy()` - X, y for S298
6. `test_select_task_view_heat_capacity()` - X, Y with 45 columns
7. `test_select_task_view_invalid()` - Expect ValueError

---

## Phase 5: Integration Tests

### File: `tests/integration/test_pipeline.py`

**Unskip `TestDataPipeline` class and implement:**

1. `test_load_and_fuse()`:
   - Load mock CSV → Merge with PM7 → Compute deltas → Verify

**Add new tests:**

2. `test_loader_to_task_view()`:
   - Loader.load() → DataFusion.select_task_view() → Verify X, y shapes

3. `test_full_data_pipeline()`:
   - Load → Filter → 3-way Split → Task View → Verify all steps

---

## Phase 6: Documentation

### File: `docs/architecture.md`

**Add section after existing content:**

```markdown
## Data Layer (v0.2.0)

### ChemperiumLoader

Loads and validates the Chemperium thermochemistry dataset (~52k molecules).

**Supported formats:** CSV, Parquet
**Required columns:** smiles, charge, multiplicity, nheavy, H298_cbs
**Optional columns:** xyz, H298_b3, S298, cp_1...cp_45

**Key methods:**
- `load(path)` - Load and validate dataset
- `split()` - 2-way train/test split
- `train_val_test_split()` - 3-way train/val/test split
- `get_features()` / `get_targets()` - Extract X, y

### DataFusion

Combines data sources and creates task-specific views.

**Key methods:**
- `merge()` - Combine Chemperium + PM7 data
- `compute_deltas()` - Calculate CBS - PM7 corrections
- `select_task_view(task)` - Get X, y for specific task:
  - `"enthalpy"` → targets H298_cbs
  - `"entropy"` → targets S298
  - `"heat_capacity"` → targets cp_1...cp_45 (multioutput)

### Extension Points (Batch 6)
- PM7/MOPAC integration via SemiempiricalHandler
- Automatic delta computation from SMILES
```

---

## Implementation Order

1. **Mock Data Fixtures** (Phase 3) - Required by tests
2. **ChemperiumLoader** (Phase 1) - Core data loading
3. **DataFusion** (Phase 2) - Task views + merge
4. **Unit Tests** (Phase 4) - Validate implementations
5. **Integration Tests** (Phase 5) - End-to-end validation
6. **Documentation** (Phase 6) - Update architecture docs

---

## Validation Criteria

After implementation, all must pass:

```bash
# All tests pass (previously 42 passed, 30 skipped)
pytest tests/ -v

# Code style
ruff check src/

# Imports work
python -c "from grimperium.data.loader import ChemperiumLoader"
python -c "from grimperium.data.fusion import DataFusion"

# Quick smoke test
python -c "
from grimperium.data import ChemperiumLoader, DataFusion
loader = ChemperiumLoader()
print('ChemperiumLoader:', loader)
fusion = DataFusion()
print('DataFusion:', fusion)
"
```

---

## Notes

- **No RDKit dependency** - Element filtering is a stub returning df unchanged
- **No PM7/MOPAC** - SemiempiricalHandler stays as stub (Batch 6)
- **Preserve existing interface** - All existing tests should pass unchanged
- **Add, don't replace** - New methods added alongside existing ones
