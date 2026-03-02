# Plan: CSV Schema Restructuring + RDKit Descriptors + Bug Fixes

## Context

The CSV tracking system for batch molecule processing has several issues:
1. RDKit descriptors are computed but never written to CSV (stored only in `result.rdkit_descriptors` dict)
2. Column order is disorganized - RDKit columns (independent of CREST/MOPAC) should come before computation results
3. `abs_diff` is redundant with `target_delta_kcalmol` (same info, different sign)
4. Interrupted batches leave molecules stuck in RUNNING status with no recovery
5. HOMO/LUMO/GAP parsing may silently fail for some molecules
6. `conformer_details/` directory files are being tracked by git

The user provided the exact 57-column target schema (see below).

---

## Target Schema (57 columns, exact order)

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
- **REMOVED:** `abs_diff`, `has_heteroatoms`, `quality_grade`, `success`, `error_message`, `total_execution_time`, `crest_error`, `mopac_time`, `precise_scf`, `scf_threshold`, `crest_optlev`, `threads`, `retry_count`, `last_error_message`, `max_retries`, all `reserved_*` columns

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
- **NEW** bond counts: `rdkit_bonds_single`, `rdkit_bonds_double`, `rdkit_bonds_triple`, `rdkit_bonds_aromatic` (iterate `mol.GetBonds()`, check `GetBondType()`)
- **REMOVE:** `logp`, `ring_count`, `sat_ring_count`, `num_heteroatoms`, `labute_asa`, `bertz_ct`

### Step 3: Rewrite CSV schema in `BatchCSVManager`
**File:** `src/grimperium/crest_pm7/batch/csv_manager.py`

Replace all class-level column lists to produce the exact 57-column schema:

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
def reset_stuck_running(self) -> int:
    """Reset RUNNING molecules to PENDING on startup recovery."""
```

Resets status to PENDING, increments `reruns`, saves CSV, returns count.

**Call sites** (add before `select_batch()`):
- `src/grimperium/cli/views/batch_view.py`
- `src/grimperium/cli/views/databases_view.py`

### Step 6: Update `csv_enhancements.py`
**File:** `src/grimperium/crest_pm7/csv_enhancements.py`

- Remove `abs_diff` from `updates` dict in `update_molecule_with_mopac_results()`
- Remove `precise_scf`, `scf_threshold`, `crest_optlev`, `threads` from updates
- Update `BatchSettingsCapture.capture_batch_settings()` to remove deleted settings

### Step 7: Harden HOMO/LUMO/GAP parser
**File:** `src/grimperium/crest_pm7/mopac_descriptors.py`

- Add negative lookbehind `(?<![A-Z_])` to EIGENVALUES regex to avoid matching `ALPHA_EIGENVALUES`
- Validate parsed eigenvalue count against declared count `EIGENVALUES[N]`
- Add bounds checking for HOMO/LUMO indices with diagnostic logging
- Upgrade silent failures from `LOG.debug` to `LOG.warning`

### Step 8: Update `models.py`
**File:** `src/grimperium/crest_pm7/batch/models.py`

Update `BatchRowCSV` Pydantic model to match new 57-column schema.

### Step 9: Update tests
**Files:**
- `tests/unit/test_csv_schema.py` - Update column count (-> 57), column lists, remove RETRY/PHASE_B tests, add RDKIT tests
- `tests/test_csv_schema.py` - Update `EXPECTED_COLUMNS` list
- `tests/test_csv_manager_retry.py` - Adapt to use `reruns` instead of `retry_count`
- `tests/unit/test_pm7_pipeline_refactor.py` - Remove `abs_diff` from test DataFrames
- `tests/test_csv_enhancements_batch_settings.py` - Remove deleted settings
- **NEW tests:** `reset_stuck_running()`, RDKit descriptor atom/bond counts

---

## Verification

1. `pytest tests/ -x` - All tests pass
2. `python -c "from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager; m = BatchCSVManager(None); print(len(m.get_schema()), m.get_schema())"` - Verify 57 columns in exact order
3. `python -c "from grimperium.core.descriptors import extract_all_rdkit_descriptors; print(extract_all_rdkit_descriptors('CCO'))"` - Verify 15 keys with `rdkit_` prefix
4. Verify `.gitignore` excludes `conformer_details/` directory
5. `ruff check src/` - No lint errors
