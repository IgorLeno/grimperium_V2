# Dataset Reference Guide

## Overview

Grimperium uses two primary datasets for delta-learning and semiempirical optimization:

1. **thermo_cbs_chon.csv** - Primary thermochemistry reference (CHON molecules only)
2. **thermo_pm7.csv** - Semiempirical optimization results (PM7 method)

---

## thermo_cbs_chon.csv (Primary)

### What is it?

A curated thermochemistry database of **29,568** small-to-medium organic molecules composed exclusively of **C, H, O, N** atoms (CHON).

### Size & Composition

- **Total molecules:** 29,568
- **Elemental composition:** C, H, O, N only
- **Removed elements:** Halogens (F, Cl, Br, I), sulfur (S), rare heteroatoms (B, P, As, Ge)

### Creation Process

Starting from the original CBS dataset of 52,837 molecules:

1. **Step 1:** Remove halogen-containing molecules (F, Cl, Br, I) → 30,026 molecules
2. **Step 2:** Remove sulfur-containing molecules (S) → 30,026 molecules (overlap with step 1)
3. **Step 3:** Remove rare heteroatoms (B, P, As, Ge) → **29,568 molecules (final)**

**Total removed:** 23,269 molecules (44% of original)
**Final dataset:** 29,568 molecules (56% of original)

### Properties Included

| Column | Description |
|--------|-------------|
| **smiles** | Molecular structure representation |
| **charge** | Total molecular charge |
| **multiplicity** | Spin multiplicity |
| **nheavy** | Number of heavy atoms |
| **H298_cbs** | Enthalpy at 298K (CBS-level, high-accuracy reference) |
| **H298_b3** | Enthalpy at 298K (B3LYP-level, cheap DFT baseline) |
| **S298** | Entropy at 298K (optional) |
| **xyz** | Optimized molecular geometry (optional) |

### Why CHON-Only?

#### Physical Reasoning
- **Homogeneous electronic physics:** No exotic valence states or hypervalence extremes
- **No relativistic effects:** Avoids heavy atoms requiring relativistic corrections
- **Well-described by semiempirical methods:** GFN-xTB, PM7, and common DFTs work well for CHON
- **No highly polarizable atoms:** Avoids challenges with dispersion corrections

#### Statistical Reasoning
- **Rare heteroatoms act as structural outliers:** Few examples prevent proper generalization
- **Imbalanced representation:** 23,269 removed molecules create a more balanced dataset
- **Remove without information loss:** Rare cases don't improve model performance

#### Methodological Reasoning (Delta-Learning)
- **Delta = CBS - B3LYP correction:** The learning target is the energy correction
- **Homogeneous error patterns are learnable:** CHON molecules share similar systematic DFT errors
- **Mixing chemical regimes violates delta-learning assumptions:** Halogens/sulfur introduce different error physics

### Use Cases

✅ **Recommended:**
- Delta-learning model training (primary use case)
- Benchmark for semiempirical methods
- Validation of thermochemistry predictions
- Test set for model evaluation

❌ **Not suitable for:**
- Modeling halogenated compounds
- Drug-like molecules with sulfur or halogens
- General-purpose thermochemistry (use full CBS dataset instead)

---

## thermo_pm7.csv (Secondary)

### What is it?

Results from **CREST conformer generation + PM7 semiempirical optimization** pipeline.

### Purpose

- Experimental validation target for delta-learning models
- Alternative baseline for comparison with CBS/B3LYP
- Production semiempirical predictions

### Properties Included

| Column | Description |
|--------|-------------|
| **smiles** | Molecular structure |
| **pm7_enthalpy** or **H298_pm7** | PM7-optimized enthalpy |
| **conformer_count** | Number of conformers generated |
| **quality_grade** | Quality assessment (A/B/C/D/F) |
| **batch_metadata** | Processing batch information |

### Use Cases

✅ **Recommended:**
- Validation of PM7-based predictions
- Alternative cheap baseline for delta-learning
- Experimental pipeline testing

❌ **Not suitable for:**
- High-accuracy thermochemistry (use CBS dataset)
- Training machine learning models (use thermo_cbs_chon.csv)

---

## Deprecated Datasets

### thermo_cbs_clean.csv ❌

**Status:** REMOVED (January 18, 2026)
**Replaced by:** `thermo_cbs_chon.csv`
**Reason:** CHON-only is more descriptive of actual content (30,026 → 29,568 molecules after rare heteroatom removal)

**Migration:**
```python
# OLD (deprecated)
df = loader.load_thermo_cbs_clean()

# NEW (use this)
df = loader.load_thermo_cbs_chon()
```

---

### thermo_cbs_opt.csv ❌

**Status:** REMOVED (January 18, 2026)
**Reason:** Not used in current workflow; `thermo_cbs_chon.csv` is sufficient

**Note:** Original unfiltered dataset (52,837 molecules) is no longer maintained in the active codebase.

---

### test_batch_final.csv ❌

**Status:** REMOVED (January 18, 2026)
**Reason:** Test data now managed via fixtures (`tests/fixtures/real_data.py`), not CSV files

---

## Loading Data in Code

### Primary dataset (CHON)

```python
from grimperium.data.loader import ChemperiumLoader

loader = ChemperiumLoader()
df = loader.load_thermo_cbs_chon(max_nheavy=50)

print(f"Loaded {len(df)} CHON molecules")
# Output: Loaded 29568 CHON molecules
```

### Secondary dataset (PM7)

```python
from grimperium.data.loader import ChemperiumLoader

loader = ChemperiumLoader()
df = loader.load_thermo_pm7(max_nheavy=50)

print(f"Loaded {len(df)} PM7 results")
```

### Using path constants directly

```python
from grimperium.data import THERMO_CBS_CHON_PATH, THERMO_PM7_PATH

print(f"Primary dataset: {THERMO_CBS_CHON_PATH}")
# Output: Primary dataset: data/thermo_cbs_chon.csv

print(f"Secondary dataset: {THERMO_PM7_PATH}")
# Output: Secondary dataset: data/thermo_pm7.csv
```

---

## Migration Guide

### For existing code using old methods:

| Old Code | New Code | Notes |
|----------|----------|-------|
| `loader.load_thermo_cbs_clean()` | `loader.load_thermo_cbs_chon()` | Method renamed |
| `loader.load_thermo_cbs_opt()` | `loader.load_thermo_cbs_chon()` | Method replaced |
| `loader.load_thermo_batch_final()` | `loader.load_thermo_pm7()` | Method renamed |
| `THERMO_CBS_CLEAN_PATH` | `THERMO_CBS_CHON_PATH` | Constant renamed |
| `THERMO_CBS_OPT_PATH` | `THERMO_CBS_CHON_PATH` | Constant replaced |

### Error handling:

If you see `FileNotFoundError` mentioning old dataset names:

```
FileNotFoundError: Primary dataset not found: data/thermo_cbs_chon.csv
Expected: 29,568 CHON molecules with CBS and B3LYP enthalpies
Note: thermo_cbs_clean.csv and thermo_cbs_opt.csv are no longer used
```

**Solution:** Use `load_thermo_cbs_chon()` instead of deprecated methods.

---

## Dataset Statistics

### thermo_cbs_chon.csv

```
Total molecules:        29,568
Average nheavy:         ~12-15 heavy atoms
N_heavy_atoms range:    1-50
Elemental composition:  C, H, O, N only
Properties:             CBS & B3LYP enthalpies, descriptors
Quality:                High (curated, homogeneous chemical space)
File size:              ~1.5 MB (CSV)
```

### thermo_pm7.csv

```
Total molecules:        Varies (CREST pipeline output)
Properties:             PM7 enthalpy, conformer data, grades
Quality:                Semiempirical (experimental)
File size:              ~2 MB (CSV, with metadata)
```

---

## File Locations

```
data/
├── thermo_cbs_chon.csv   ← Primary dataset (29,568 CHON molecules)
└── thermo_pm7.csv         ← Secondary dataset (PM7 optimization results)
```

**Note:** Only these two CSV files should exist in `data/`. Any additional CSVs indicate orphaned files from previous versions.

---

## Validation & Quality Assurance

### Expected file structure:

```bash
# Verify dataset files
ls -lh data/*.csv

# Expected output:
# -rw-rw-r-- 1 user user 1.5M data/thermo_cbs_chon.csv
# -rw-rw-r-- 1 user user 2.0M data/thermo_pm7.csv
```

### Programmatic validation:

```python
from pathlib import Path
from grimperium.data import THERMO_CBS_CHON_PATH, THERMO_PM7_PATH

# Check files exist
assert Path(THERMO_CBS_CHON_PATH).exists(), "Primary dataset missing"
assert Path(THERMO_PM7_PATH).exists(), "Secondary dataset missing"

# Load and validate
from grimperium.data.loader import ChemperiumLoader

loader = ChemperiumLoader(validate=True)
df_chon = loader.load_thermo_cbs_chon()
df_pm7 = loader.load_thermo_pm7()

print(f"✓ Loaded {len(df_chon)} CHON molecules")
print(f"✓ Loaded {len(df_pm7)} PM7 results")
```

---

## References

- **CHON dataset curation:** Internal process (January 2026)
- **Delta-learning framework:** See `docs/DELTALEARNING.md` (if available)
- **CREST + PM7 pipeline:** See `docs/CREST_INTEGRATION.md`
- **Dataset migration:** See `CHANGELOG.md` (January 18, 2026 entry)

---

## FAQ

### Q: Why not include all 52,837 molecules?

**A:** Halogens, sulfur, and rare heteroatoms introduce:
- Different DFT error physics (polarization, dispersion)
- Outlier cases that harm delta-learning generalization
- Heterogeneous chemical space that violates delta-learning assumptions

### Q: Can I use this dataset for halogenated molecules?

**A:** No. Models trained on CHON-only data will perform poorly on molecules with F, Cl, Br, I, or S. Use a different dataset for those cases.

### Q: What happened to `thermo_cbs_clean.csv`?

**A:** Renamed to `thermo_cbs_chon.csv` to be more descriptive. After removing rare heteroatoms, the dataset decreased from 30,026 → 29,568 molecules.

### Q: How do I test with real data in CI?

**A:** Use `tests/fixtures/real_data.py`:

```python
from tests.fixtures.real_data import load_real_subset

# Load 1000 molecules for fast CI tests
df = load_real_subset(n=1000, stratified=True, random_state=42)
```

---

**Last Updated:** January 18, 2026
**Maintained By:** Grimperium Team
**Review Cycle:** Quarterly (or on major dataset changes)
