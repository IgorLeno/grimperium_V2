# Dataset Migration Guide

## Overview

On 2026-01-10, Grimperium migrated from `thermo_cbs_opt.csv` (original) to `thermo_cbs_clean.csv` (filtered).

## Comparison

| Aspect | thermo_cbs_opt.csv (OLD) | thermo_cbs_clean.csv (NEW) |
|--------|------------------------|--------------------------|
| **Source** | Original CBS dataset | Cleaned/filtered version |
| **Molecules** | All (original count) | ~30,026 (filtered) |
| **Halogenated** | ✅ Included | ❌ Removed |
| **Sulfur** | ✅ Included | ❌ Removed |
| **Non-essential cols** | ✅ Included | ❌ Removed |
| **Phase A Use** | ⚠️ Not recommended | ✅ Recommended |
| **Status** | ⚠️ Deprecated | ✅ ACTIVE |

## Filtering Applied

### Halogenated Molecules Removed
Reason: Not in scope for Phase A study

Examples removed:
- Chlorinated: `Cl-CH2-CH3` (chloroethane)
- Brominated: `Br-CH2-CH3` (bromoethane)
- Fluorinated: `F-CH2-CH3` (fluoroethane)
- Iodinated: `I-CH2-CH3` (iodoethane)

### Sulfur-Containing Molecules Removed
Reason: Not in scope for Phase A study

Examples removed:
- Thiols: `CH3-SH` (methanethiol)
- Sulfides: `CH3-S-CH3` (dimethyl sulfide)
- Disulfides: `CH3-S-S-CH3` (dimethyl disulfide)

### Non-Essential Columns Removed
Kept: molecule_id, smiles, thermodynamic_properties

Removed: metadata columns not needed for Phase A

## Migration Steps

### Step 1: Update Imports
```python
from grimperium.data import ChemperiumLoader

loader = ChemperiumLoader()
```

### Step 2: Use New Method
```python
# ❌ OLD (shows deprecation warning)
df = loader.load_thermo_cbs_opt()

# ✅ NEW (recommended)
df = loader.load_thermo_cbs_clean()
```

### Step 3: If You Need Original Data
```python
# Explicit original file
df_original = loader.load_thermo_cbs_opt(
    path="data/thermo_cbs_opt.csv"
)
```

## FAQ

**Q: Can I still use the original dataset?**  
A: Yes, call `load_thermo_cbs_opt(path="data/thermo_cbs_opt.csv")`

**Q: When will load_thermo_cbs_opt() be removed?**  
A: Timeline: 2026-01-10 (deprecated) → 2026-06-10 (formal deprecation) → 2026-12-10 (possible removal)

**Q: Why remove halogenated molecules?**  
A: Phase A scope focuses on organic molecules without halogens. Simplifies initial validation.

**Q: How many molecules were filtered out?**  
A: See dataset statistics in Phase A docs (exact count TBD)

**Q: Does this affect my Phase B code?**  
A: Phase B models should train on clean dataset. If you need original data, use explicit path.

## Deprecation Timeline

| Date | Event |
|------|-------|
| 2026-01-10 | `load_thermo_cbs_clean()` introduced as primary method |
| 2026-01-10 | `load_thermo_cbs_opt()` default changed to clean dataset |
| 2026-01-10 | DeprecationWarning added to `load_thermo_cbs_opt()` |
| 2026-06-10 | `load_thermo_cbs_opt()` officially marked deprecated |
| 2026-12-10 | `load_thermo_cbs_opt()` may be removed (TBD) |

## Impact Assessment

### Code That Needs Updates

1. **Test Fixtures** (`tests/fixtures/real_data.py`)
   - Update to use `load_thermo_cbs_clean()`
   - Document dataset version in comments

2. **Integration Tests** (`tests/integration/`)
   - Replace `load_thermo_cbs_opt()` with `load_thermo_cbs_clean()`
   - Update assertions for filtered molecule count

3. **Scripts** (`scripts/`)
   - Phase A scripts: use clean dataset
   - Benchmarking scripts: may need original dataset

### Code That Doesn't Need Updates

1. **Core API** (`grimperium/api.py`)
   - No direct dataset loading

2. **Models** (`grimperium/models/`)
   - Dataset-agnostic

3. **Utils** (`grimperium/utils/`)
   - No dataset dependencies

## Release Notes Entry

### Breaking Changes (v0.2.0 - 2026-01-10)

**Dataset Migration: thermo_cbs_opt.csv → thermo_cbs_clean.csv**

The default dataset used by `ChemperiumLoader.load_thermo_cbs_opt()` has changed:

**OLD BEHAVIOR:**
- Default: `data/thermo_cbs_opt.csv` (original, unfiltered)
- Contains: All molecules, non-essential columns

**NEW BEHAVIOR:**
- Default: `data/thermo_cbs_clean.csv` (filtered for Phase A)
- Filtering applied:
  - ❌ Halogenated molecules removed (not in Phase A scope)
  - ❌ Sulfur-containing molecules removed (not in Phase A scope)
  - ❌ Non-essential columns removed
- ✅ ~30,026 molecules (filtered subset)

**MIGRATION GUIDE:**

Option 1: Use new explicit method (RECOMMENDED)
```python
# OLD (shows DeprecationWarning)
df = loader.load_thermo_cbs_opt()

# NEW (recommended)
df = loader.load_thermo_cbs_clean()
```

Option 2: Load original dataset explicitly
```python
# If you need the original unfiltered data:
df = loader.load_thermo_cbs_opt(path="data/thermo_cbs_opt.csv")
```

**Deprecation Timeline:**
- 2026-01-10: New method introduced, deprecation warning added
- 2026-06-10: `load_thermo_cbs_opt()` officially marked deprecated
- 2026-12-10: May be removed in future version

**Why This Change?**
Phase A validation requires molecules in study scope:
- Focus: Organic molecules without halogens or sulfur
- Rationale: Simplify initial validation, reduce dataset complexity
- Impact: Production-ready dataset (30K molecules) for Phase A onwards

## Technical Details

### File Locations
- Original: `data/thermo_cbs_opt.csv` (legacy)
- Cleaned: `data/thermo_cbs_clean.csv` (active)

### Column Differences

#### thermo_cbs_opt.csv (Original)
- All original columns from CBS dataset
- Non-essential metadata included

#### thermo_cbs_clean.csv (Cleaned)
Essential columns only:
- `molecule_id`: Unique identifier
- `smiles`: SMILES representation
- `charge`: Total molecular charge
- `multiplicity`: Spin multiplicity
- `nheavy`: Number of heavy atoms
- `H298_cbs`: Enthalpy at 298K (CBS reference)
- Other thermodynamic properties

### Filtering Logic

Halogenated molecule detection:
```python
# Pseudocode for filtering
smiles_contains_any(['Cl', 'Br', 'F', 'I'])
```

Sulfur molecule detection:
```python
# Pseudocode for filtering
smiles_contains('S')
```

## Support

If you encounter issues during migration:

1. Check the FAQ above
2. Review the [Phase A documentation](PHASE-A-START-HERE.md)
3. Open an issue on GitHub with the `dataset-migration` label
4. Contact the development team

## Acknowledgments

This migration was implemented to align the codebase with Phase A validation requirements and improve dataset quality for production use.
