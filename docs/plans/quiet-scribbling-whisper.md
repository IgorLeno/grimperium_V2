# Plan: CSV Schema Restructuring + RDKit Descriptors + Bug Fixes

## Context

The CSV tracking system for batch molecule processing has several issues:
1. RDKit descriptors are computed but never written to CSV (stored only in `result.rdkit_descriptors` dict)
2. Column order is disorganized - RDKit columns (independent of CREST/MOPAC) should come before computation results
3. `abs_diff` is redundant with `target_delta_kcalmol` (same info, different sign)
4. Interrupted batches leave molecules stuck in RUNNING status with no recovery
5. HOMO/LUMO/GAP parsing may silently fail for some molecules
6. `conformer_details/` directory files are being tracked by git

The user provided the exact 58-column target schema (see below).

---

## Target Schema (58 columns, exact order)

```
mol_id, status, smiles, multiplicity, charge, nheavy, H298_cbs, H298_pm7,
target_delta_kcalmol, abs_diff_%, batch_id, timestamp, total_time, reruns,
rdkit_nrotbonds, rdkit_tpsa, rdkit_num_rings, rdkit_fsp3, rdkit_mol_weight,
rdkit_hbond_donors, rdkit_hbond_acceptors, rdkit_nC, rdkit_nH, rdkit_nO, rdkit_nN,
rdkit_bonds_single, rdkit_bonds_double, rdkit_bonds_triple, rdkit_bonds_aromatic,
crest_status, xtb, v3, qm, nci, c_method, energy_window, rmsd_threshold, opt_lvl,
crest_conformers_generated, num_conformers_selected, crest_time_s,
mopac_status, k_selected_pm7, mopac_dipole_debye, mopac_ionization_potential_ev,
mopac_homo_ev, mopac_lumo_ev, mopac_gap_ev, mopac_cosmo_area_a2, mopac_cosmo_volume_a3,
mopac_gradient_norm, mopac_num_scf_cycles, mopac_point_group, mopac_time_s,
batch_order, batch_failure_policy, assigned_crest_timeout, assigned_mopac_timeout
```

### Changes vs Current Schema
- **NEW columns:** `multiplicity`, `charge`, 15x `rdkit_*` columns (atom/bond counts)
- **RENAMED:** `reference_hof` -> `H298_cbs`, `crest_time` -> `crest_time_s`
- **REMOVED:** `abs_diff`, `has_heteroatoms`, `quality_grade`, `success`, `error_message`, `total_execution_time`, `crest_error`, `precise_scf`, `scf_threshold`, `crest_optlev`, `threads`, `retry_count`, `last_error_message`, `max_retries`, all `reserved_*` columns

---

## Migration Strategy for Existing CSVs

### Schema Version Field
- Write a `schema_version` column (integer) to every CSV when first writing with the new schema. Set `schema_version = 2` for new files; absence of this field or value `1` indicates the old format.
- The loader (`CSVDataLoader`) must check for `schema_version` and emit a `DeprecationWarning` if a v1 CSV is loaded without migration.

### Idempotent Migration Script
Create `scripts/migrate_csv_schema.py` with the following behaviour:
1. **Column renames**: copy values from old names to new (`reference_hof` → `H298_cbs`, `crest_time` → `crest_time_s`, `mopac_time` → `mopac_time_s`, `nrotbonds` → `rdkit_nrotbonds`).
2. **New columns**: fill `multiplicity` with `1`, `charge` with `0`, all 15 `rdkit_*` columns with `None`/empty, `schema_version` with `2`.
3. **Removed columns**: drop `abs_diff`, `has_heteroatoms`, `quality_grade`, `success`, `error_message`, `total_execution_time`, `crest_error`, `precise_scf`, `scf_threshold`, `crest_optlev`, `threads`, `retry_count`, `last_error_message`, `max_retries`, all `reserved_*` columns.
4. **Column reorder**: rewrite CSV with columns in the exact 58-column order from `BatchCSVManager.get_schema()`.
5. The script is idempotent: running it twice on a v2 CSV produces no change.

### Partially-Completed Batches
- Rows with `status == RUNNING` at migration time: **preferred approach is option 2 (no schema change)**. Under an advisory file lock, reset `status` from `RUNNING` to `PENDING` and increment `reruns`. Do not introduce `migration_state` in the target 58-column schema unless a separate migration-tracking design is explicitly approved.
- Rows with `status == OK` or `status == SKIP` are migrated as-is (rename/fill columns, preserve data).
- Any row whose migration produces inconsistent data (e.g., non-numeric where numeric expected) should be written to a `migration_errors.csv` sidecar file for manual review.

### Migration Mode: Batch (one-time) vs Automatic On-Load
- **Recommended**: one-time batch migration run before restarting the pipeline. This avoids coupling migration logic into the hot path.
- **On-load migration**: `CSVDataLoader._add_missing_columns()` already fills missing columns with defaults (see `OPTIONAL_COLUMNS` dict) and maps old column names to new ones (see alias mapping in `_add_missing_columns`). This provides transparent backwards compatibility for small/incremental differences.
- **Rollback**: before migration, target backup `<filename>.pre_migration_backup.csv` and first check whether it already exists. If it exists, log and abort with a clear error that references the backup filename and instructs removing or versioning the existing backup before continuing. In this error path, do **not** call `BatchCSVManager.get_schema()` and do not run migration. If backup creation succeeds, continue; if post-migration validation fails (run `BatchCSVManager.get_schema()` and verify all columns present; check row count matches), restore from backup and abort.

---

## Implementation Steps

### Step 1: `.gitignore` fix
**File:** `.gitignore`
- Change `data/molecules_pm7/conformer_details/*.json` to `data/molecules_pm7/conformer_details/`
- Run `git rm --cached data/molecules_pm7/conformer_details/` to untrack existing files

### Step 2: Update `extract_all_rdkit_descriptors()`
**File:** `src/grimperium/core/descriptors.py`

Replace the 13-key dict with 15 keys matching the CSV schema:
- Keep (renamed): `nrotbonds` -> `rdkit_nrotbonds`, `tpsa` -> `rdkit_tpsa`, `fraction_csp3` -> `rdkit_fsp3`, `mol_wt` -> `rdkit_mol_weight`, `hbd` -> `rdkit_hbond_donors`, `hba` -> `rdkit_hbond_acceptors`, `arom_ring_count` -> `rdkit_num_rings`
- **NEW** atom counts: `rdkit_nC`, `rdkit_nH`, `rdkit_nO`, `rdkit_nN` (iterate `Chem.AddHs(mol).GetAtoms()`)
  - Only count atoms with symbol in `{"C", "H", "O", "N"}`; increment the corresponding counter.
  - Do **not** count S, P, halogens, or other elements in any of the four counters.
  - If any atom with a symbol outside `{C, H, O, N}` is encountered, emit `LOG.warning("Molecule contains non-CHON atom: %s", symbol)` so callers are aware; do not crash.
- **NEW** bond counts: `rdkit_bonds_single`, `rdkit_bonds_double`, `rdkit_bonds_triple`, `rdkit_bonds_aromatic` (iterate `mol.GetBonds()`, check `GetBondType()`)
- **REMOVE:** `logp`, `ring_count`, `sat_ring_count`, `num_heteroatoms`, `labute_asa`, `bertz_ct`

### Step 3: Rewrite CSV schema in `BatchCSVManager`
**File:** `src/grimperium/crest_pm7/batch/csv_manager.py`

Replace all class-level column lists to produce the exact 58-column schema:

```python
IDENTITY_COLUMNS = ["mol_id", "status", "smiles"]

MOLECULAR_PROPERTIES_COLUMNS = [
    "multiplicity", "charge", "nheavy", "H298_cbs",
    "H298_pm7", "target_delta_kcalmol", "abs_diff_%",
]

BATCH_INFO_COLUMNS = ["batch_id", "timestamp", "total_time", "reruns"]

RDKIT_COLUMNS = [  # NEW list
    "rdkit_nrotbonds", "rdkit_tpsa", "rdkit_num_rings", "rdkit_fsp3",
    "rdkit_mol_weight", "rdkit_hbond_donors", "rdkit_hbond_acceptors",
    "rdkit_nC", "rdkit_nH", "rdkit_nO", "rdkit_nN",
    "rdkit_bonds_single", "rdkit_bonds_double", "rdkit_bonds_triple",
    "rdkit_bonds_aromatic",
]

CREST_COLUMNS = [
    "crest_status", "xtb", "v3", "qm", "nci", "c_method",
    "energy_window", "rmsd_threshold", "opt_lvl",
    "crest_conformers_generated", "num_conformers_selected", "crest_time_s",
]

MOPAC_COLUMNS = [
    "mopac_status", "k_selected_pm7",
    "mopac_dipole_debye", "mopac_ionization_potential_ev",
    "mopac_homo_ev", "mopac_lumo_ev", "mopac_gap_ev",
    "mopac_cosmo_area_a2", "mopac_cosmo_volume_a3",
    "mopac_gradient_norm", "mopac_num_scf_cycles", "mopac_point_group",
    "mopac_time_s",
]

BATCH_CONFIG_COLUMNS = [
    "batch_order", "batch_failure_policy",
    "assigned_crest_timeout", "assigned_mopac_timeout",
]
```

**Remove:** `RETRY_TRACKING_COLUMNS`, `PHASE_B_RESERVED_COLUMNS`, `MOPAC_CONFIG_COLUMNS`, `CREST_CONFIG_COLUMNS` (replaced by new names)

**Update `get_schema()`:** Concatenate new lists in order.

**Update `load_csv()`:** Remove dtype refs for deleted columns (`quality_grade`, `crest_error`, `error_message`, `last_error_message`). Add `"mopac_status": str`.

**Update `get_reference_hof()`:** Change `reference_hof` to `H298_cbs`.

**Update `select_batch()`:** Change `row.get("nrotbonds")` to `row.get("rdkit_nrotbonds")`.

**Update `_apply_sorting_strategy()`:** Change `"nrotbonds"` to `"rdkit_nrotbonds"`.

**Update `mark_rerun()`/`mark_skip()`:** Remove `retry_count`/`max_retries`/`last_error_message` writes. Use `reruns` column instead.

### Step 4: Rewrite `pm7result_to_csv_update()`
**File:** `src/grimperium/crest_pm7/batch/csv_manager.py`

- Map all 15 `rdkit_*` columns from `result.rdkit_descriptors` dict
- Remove `abs_diff` calculation (keep `abs_diff_%`)
- Rename: `crest_time` -> `crest_time_s`
- Remove: `quality_grade`, `success`, `error_message`, `total_execution_time`, `crest_error`, `mopac_time`
- MOPAC descriptor columns still set to `None` (filled by csv_enhancements)

### Step 5: Add `reset_stuck_running()`
**File:** `src/grimperium/crest_pm7/batch/csv_manager.py`

```python
def reset_stuck_running(self, stale_after_seconds: float = 3600.0) -> int:
    """Reset stale RUNNING molecules to PENDING on startup recovery.

    Only resets rows that have been in RUNNING status longer than
    `stale_after_seconds` (default: 1 hour), using the `timestamp` column
    as the start-of-processing reference.  Rows with a `pid` column set to a
    currently-running OS process are treated as owned and skipped.

    Multi-worker requirement: callers MUST acquire an advisory file lock
    (e.g. `filelock.FileLock`) before calling this method and MUST keep the
    same lock held before any subsequent `select_batch()` call so exactly one
    process performs startup recovery + selection.

    Args:
        stale_after_seconds: Rows in RUNNING status for longer than this many
            seconds are considered stuck and reset (default 3600 = 1 hour).

    Returns:
        Number of rows reset to PENDING.
    """
```

Implementation details:
- Parse and normalize `timestamp` to timezone-aware UTC before comparing to `datetime.now(UTC)`: parse with `datetime.fromisoformat(...)`, treat naive timestamps as UTC by attaching `timezone.utc`, then convert via `.astimezone(timezone.utc)`. Compute `elapsed = (now_utc - ts_utc).total_seconds()` and skip when `elapsed < stale_after_seconds` (no naive-vs-aware comparisons).
- If a `pid` column exists, call `psutil.pid_exists(pid)` (or `os.kill(pid, 0)`) and skip rows whose process is still alive.
- Rows that are reset: set `status = PENDING`, increment `reruns`, clear `timestamp`.
- Wrap CSV read/write with the same advisory file lock (`filelock.FileLock`) used by the caller, so recovery and immediate selection preserve single-process semantics.
- Save CSV and return count of reset rows.

**Call sites** (add before `select_batch()`, under the same file lock):
- `src/grimperium/cli/views/batch_view.py`
- `src/grimperium/cli/views/databases_view.py`

### Step 6: Update `csv_enhancements.py`
**File:** `src/grimperium/crest_pm7/csv_enhancements.py`

- Remove `abs_diff` from `updates` dict in `update_molecule_with_mopac_results()`
- Remove `precise_scf`, `scf_threshold`, `crest_optlev`, `threads` from updates
- Update `BatchSettingsCapture.capture_batch_settings()` to remove deleted settings

### Step 7: Harden HOMO/LUMO/GAP parser
**File:** `src/grimperium/crest_pm7/mopac_descriptors.py`

- Add negative lookbehind `(?<![A-Z_])` to the EIGENVALUES regex so that tokens like `ALPHA_EIGENVALUES` are not matched (only bare `EIGENVALUES[N]` should trigger the parser).
- After extracting eigenvalues, compare the parsed count to the declared count `N` from `EIGENVALUES[N]`; if they differ emit `LOG.warning("EIGENVALUES count mismatch: declared %d, parsed %d", declared, parsed)` and treat the result as unreliable (return `None` for HOMO/LUMO/GAP).
- Before indexing the eigenvalues array for HOMO index `homo_idx` and LUMO index `lumo_idx = homo_idx + 1`, assert `0 <= homo_idx < len(eigenvalues)` and `lumo_idx < len(eigenvalues)`; if bounds fail emit a `LOG.warning` with the molecule identifier and the bad index value, and return `None` for the affected descriptor(s).
- Promote any existing `LOG.debug` calls that indicate parsing failures (e.g., regex no-match, empty eigenvalue list) to `LOG.warning` so issues are visible in normal log output.

### Step 8: Update `models.py`
**File:** `src/grimperium/crest_pm7/batch/models.py`

Update `BatchRowCSV` Pydantic model to match new 58-column schema.

### Step 9: Update tests
**Files:**
- `tests/unit/test_csv_schema.py` - Update column count (→ 58), fix docstring breakdown, update column lists, remove RETRY/PHASE_B tests, add RDKIT tests
- `tests/test_csv_schema.py` - Update `EXPECTED_COLUMNS` list (58 entries) and `reruns` semantics (not `retry_count`)
- `tests/test_csv_manager_retry.py` - Adapt to use `reruns` instead of `retry_count`; ensure rows written in helper `_write_csv` include `crest_status` and `mopac_status` columns
- `tests/unit/test_pm7_pipeline_refactor.py` - Remove `abs_diff` from test DataFrames
- `tests/test_csv_enhancements_batch_settings.py` - Remove deleted settings

**NEW tests to add:**
- **Full pipeline integration test**: exercise CSV ingestion through `BatchCSVManager.load_csv()` → `select_batch()` → `mark_ok()` / `mark_rerun()` with a 58-column CSV (matching `get_schema()`); assert rerun semantics use `reruns` column.
- **RDKit edge-case test**: call `extract_all_rdkit_descriptors()` with SMILES containing S/P/halogen atoms (e.g., `"CCCS"` for propanethiol); assert `rdkit_nC`, `rdkit_nH`, `rdkit_nO`, `rdkit_nN` contain only C/H/O/N counts, a warning was emitted, and no KeyError is raised.
- **Concurrency test** (`tests/test_csv_manager_concurrency.py`): spawn two threads both calling `reset_stuck_running()` simultaneously on the same CSV; assert final `reruns` count equals exactly 1 (not 2), explicitly validating that `filelock.FileLock` prevents concurrent recovery and enforces single-process semantics.
- **Backward-compatibility test** (`tests/test_csv_compatibility.py`): write a minimal old-format CSV (with columns `reference_hof`, `crest_time`, `nrotbonds`) and load it through `CSVDataLoader`; assert the loader produces clear, actionable error or deprecation messages and that data is preserved under `H298_cbs`, `crest_time_s`, `rdkit_nrotbonds`.

---

## Verification

1. `pytest tests/ -x` - All tests pass
2. `python -c "from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager; m = BatchCSVManager(None); print(len(m.get_schema()), m.get_schema())"` - Verify 58 columns in exact order
3. `python -c "from grimperium.core.descriptors import extract_all_rdkit_descriptors; print(extract_all_rdkit_descriptors('CCO'))"` - Verify 15 keys with `rdkit_` prefix and expected value ranges (e.g., `rdkit_nC == 2`, `rdkit_nO == 1` for ethanol)
4. `python -c "from grimperium.core.descriptors import extract_all_rdkit_descriptors; print(extract_all_rdkit_descriptors('c1ccccc1'))"` - Spot-check benzene: `rdkit_nC == 6`, `rdkit_nH == 6`, `rdkit_bonds_aromatic == 6`
5. Verify `.gitignore` excludes `conformer_details/` directory
6. `ruff check src/` - No lint errors
7. **CSV I/O integration**: instantiate `BatchCSVManager`, load a real CSV, modify a row, save, reload; assert row count and column set match the 58-column schema from `get_schema()`.
8. **HOMO/LUMO regression**: run the MOPAC descriptor parser on previously failing molecules (molecules where `ALPHA_EIGENVALUES` appears in output) and confirm no `None` is returned due to the regex fix.
9. **Performance micro-benchmark**: measure `BatchCSVManager.load_csv()` + `save_csv()` round-trip for a 1000-row CSV with the new 58-column schema; ensure it completes under 2 seconds (to confirm no regression from added columns).
10. **`reset_stuck_running()` idempotency**: call `reset_stuck_running()` twice on the same CSV with one stale RUNNING row; assert count is `1` on first call and `0` on second call (no double-increment of `reruns`).
