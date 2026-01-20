# 📋 DOCUMENTO 1: GRIMPERIUM SPEC

# GRIMPERIUM V2 - PHASE A COMPLETION SPECIFICATION

**Status:** ✅ Final Specification  
**Date:** 2026-01-20  
**Version:** 1.0  
**Target:** Phase A Validation & Fixes  

---

## EXECUTIVE SUMMARY

Grimperium V2 has completed Phase A validation with 3 test molecules. This spec captures all required refinements to enable Phase B (29k molecule dataset preparation).

**Key Changes:**
- CSV schema: 35 → 49 columns
- Settings UI: Toggles → Dropdowns
- Data validation: Proper delta calculations
- Documentation: Add RDKit integration guide
- Tracking: Settings persistence + rerun counter

---

## PART 1: CSV SCHEMA REFINEMENT

### 1.1 Current State (After Phase A)

```
mol_id | status | smiles | ... | most_stable_hof | batch_order | batch_failure_policy
mol_00001 | OK | C=Cc1ccc(C)o1 | ... | -3.68 | 1.0 | partial_ok
mol_00002 | OK | Cc1ccoc1 | ... | -14.66 | 2.0 | partial_ok
mol_00003 | OK | Cc1cccc2occc12 | ... | -0.14 | 3.0 | partial_ok
```

**Issues:**
- ❌ Column name `most_stable_hof` doesn't match pattern (should be `H298_pm7` like `H298_cbs`)
- ❌ Columns `crest_timeout` and `mopac_timeout` duplicate `assigned_crest_timeout` and `assigned_mopac_timeout`
- ❌ Missing metrics for comparing CBS vs PM7
- ❌ Missing settings tracking (CREST/MOPAC configuration used)
- ❌ Missing rerun counter for failure handling
- ❌ Delta names unclear (`delta_e_12` vs `delta_e_13` vs `delta_e_15`)
- ❌ No conformer selection tracking

### 1.2 Final Schema (49 Columns)

**NEW COLUMN ORDER:**

```
1.   mol_id
2.   status
3.   smiles
4.   multiplicity
5.   charge
6.   nheavy
7.   H298_cbs
8.   H298_pm7              ← Renamed from most_stable_hof
9.   abs_diff              ← NEW: |H298_cbs - H298_pm7|
10.  abs_diff_%            ← NEW: Percentage difference
11.  batch_id
12.  timestamp
13.  reruns                ← NEW: Retry counter (0-3)
14.  nrotbonds             ← RDKit descriptor
15.  tpsa                  ← RDKit descriptor
16.  aromatic_rings        ← RDKit descriptor
17.  crest_status
18.  xtb                   ← NEW: CREST setting
19.  v3                    ← NEW: CREST setting
20.  qm                    ← NEW: CREST quick mode (renamed from quick)
21.  nci                   ← NEW: CREST setting
22.  c_method              ← NEW: CREST method (from crest_gfnff toggle)
23.  energy_window         ← NEW: CREST setting
24.  rmsd_threshold        ← NEW: CREST setting
25.  threads               ← NEW: CREST setting
26.  crest_conformers_generated
27.  crest_time
28.  num_conformers_selected
29.  mopac_status
30.  precise_scf           ← NEW: MOPAC setting
31.  scf_threshold         ← NEW: MOPAC setting
32.  mopac_time            ← FIX: Currently empty, must be aggregated
33.  delta_1               ← NEW: |H298_cbs - hof_1| (was delta_e_12)
34.  delta_2               ← NEW: |H298_cbs - hof_2| (was delta_e_13)
35.  delta_3               ← NEW: |H298_cbs - hof_3| (was delta_e_15)
36.  conformer_selected    ← NEW: Which conformer used (1/2/3)
37.  error_message
38.  batch_order
39.  batch_failure_policy
40.  assigned_crest_timeout
41.  assigned_mopac_timeout
42-49: [Reserved for future metrics - to be defined in Phase B]
```

**COLUMNS TO REMOVE:**
- `crest_timeout` (duplicate of `assigned_crest_timeout`)
- `mopac_timeout` (duplicate of `assigned_mopac_timeout`)

**COLUMNS TO RENAME:**
- `most_stable_hof` → `H298_pm7` (matches `H298_cbs` naming pattern)
- `delta_e_12` → `delta_1`
- `delta_e_13` → `delta_2`
- `delta_e_15` → `delta_3`
- `quick` → `qm` (in settings, to match CREST terminology)
- `gfnff` → `c_method` (in settings, clearer naming)

### 1.3 Column Definitions

#### Core Molecular Properties
| Column | Type | Source | Definition |
|--------|------|--------|-----------|
| mol_id | str | Input | Unique molecule identifier (e.g., mol_00001) |
| smiles | str | Input | SMILES string representation |
| multiplicity | int | Input | Spin multiplicity (1=singlet, 2=doublet, etc) |
| charge | int | Input | Molecular charge |
| nheavy | int | Input | Number of heavy atoms (non-hydrogen) |

#### Reference Data (from database)
| Column | Type | Source | Definition |
|--------|------|--------|-----------|
| H298_cbs | float | External DB | Reference enthalpy at 298K from CBS-QB3 (kcal/mol) |

#### RDKit Descriptors (NEW)
| Column | Type | Source | Definition |
|--------|------|--------|-----------|
| nrotbonds | int | RDKit | Number of rotatable bonds |
| tpsa | float | RDKit | Topological Polar Surface Area (Ų) |
| aromatic_rings | int | RDKit | Count of aromatic rings |

#### PM7 Calculation Results
| Column | Type | Source | Definition |
|--------|------|--------|-----------|
| H298_pm7 | float | MOPAC | Enthalpy at 298K from selected PM7 conformer (kcal/mol) |
| crest_conformers_generated | int | CREST | Total conformers created by CREST |
| num_conformers_selected | int | CREST→MOPAC | Top-N conformers selected for MOPAC optimization |
| mopac_time | float | MOPAC | Total execution time for MOPAC (aggregated across all selected conformers, seconds) |
| crest_time | float | CREST | CREST execution time (seconds) |

#### Metrics (NEW)
| Column | Type | Calculation | Definition |
|--------|------|-------------|-----------|
| abs_diff | float | \|H298_cbs - H298_pm7\| | Absolute difference in enthalpies (kcal/mol) |
| abs_diff_% | float | (abs_diff / \|H298_cbs\|) × 100 | Percentage error relative to CBS value |

#### Delta Calculations (NEW)
| Column | Type | Calculation | Definition |
|--------|------|-------------|-----------|
| delta_1 | float | \|H298_cbs - hof_conformer_1\| | Energy difference vs best conformer (lowest energy) |
| delta_2 | float | \|H298_cbs - hof_conformer_2\| | Energy difference vs 2nd best conformer |
| delta_3 | float | \|H298_cbs - hof_conformer_3\| | Energy difference vs 3rd best conformer |
| conformer_selected | int | argmin([delta_1, delta_2, delta_3]) | Index of selected conformer (1=best, 2=2nd, 3=3rd) |

#### Conformer Selection Logic
The best conformer is selected by:
1. Collect successful HOF values from all conformers
2. Sort by HOF energy (lowest first)
3. Take top 3 conformers
4. Calculate delta for each: δ = |H298_cbs - hof|
5. **Select conformer with MINIMUM |delta|** (best agreement with CBS)
6. Use that conformer's HOF as H298_pm7
7. Store all deltas and conformer index

**Example:**
```
CREST Generated: 5 conformers
MOPAC Optimized: 5 conformers successful

HOF values (sorted): [-3.68, -3.52, -3.45, -2.90, -2.80]
Top 3: [-3.68, -3.52, -3.45]

H298_cbs = 0.15234

Deltas:
  delta_1 = |0.15234 - (-3.68)| = 3.83234
  delta_2 = |0.15234 - (-3.52)| = 3.67234
  delta_3 = |0.15234 - (-3.45)| = 3.60234  ← MINIMUM

conformer_selected = 3  (third best conformer has best agreement)
H298_pm7 = -3.45       (use that conformer's HOF)
abs_diff = 3.60234
abs_diff_% = (3.60234 / 0.15234) × 100 = 2364.32%
```

#### CREST Configuration (NEW)
| Column | Type | Source | Definition |
|--------|------|--------|-----------|
| xtb | bool | Settings | Use xTB method (true/false) |
| v3 | bool | Settings | Use v3 algorithm (true/false) |
| qm | str | Settings | Quick mode: "off", "quick", "squick", or "mquick" |
| nci | bool | Settings | Use NCI mode (true/false) |
| c_method | str | Settings | CREST QM method: "gfn2", "gfnff", or "gfn2//gfnff" |
| energy_window | float | Settings | Energy window in kcal/mol for conformer selection |
| rmsd_threshold | float | Settings | RMSD threshold for duplicate removal |
| threads | int | Settings | CPU threads allocated to CREST |

#### MOPAC Configuration (NEW)
| Column | Type | Source | Definition |
|--------|------|--------|-----------|
| precise_scf | bool | Settings | Use precise SCF convergence (true/false) |
| scf_threshold | float | Settings | SCF convergence threshold (e.g., 1e-4) |

#### Status & Tracking
| Column | Type | Definition | Valid Values |
|--------|------|-----------|--------------|
| status | str | Molecule processing state | Pending, Current Batch, Running, OK, Rerun, Skip |
| reruns | int | Rerun attempt counter | 0-3 (max 3 retries) |
| error_message | str | Error description if failed | Any string or empty |
| batch_id | str | Batch identifier | e.g., batch_0001 |
| timestamp | str | Processing timestamp | Format: DD/MM-HH:MM |
| batch_order | int | Processing order within batch | 1.0, 2.0, 3.0, ... |
| batch_failure_policy | str | How to handle partial failures | "strict", "partial_ok", "lenient" |

#### CREST/MOPAC Timeouts (OLD - kept for reference)
| Column | Type | Definition |
|--------|------|-----------|
| assigned_crest_timeout | float | Timeout allocated to CREST (minutes) |
| assigned_mopac_timeout | float | Timeout allocated to MOPAC (minutes) |

#### Status Information
| Column | Type | Definition | Example |
|--------|------|-----------|---------|
| crest_status | str | CREST result status | "SUCCESS", "TIMEOUT", "FAILED" |
| mopac_status | str | MOPAC result status | "SUCCESS", "TIMEOUT", "FAILED" |

### 1.4 Test Data Validation

**Expected values for 3 test molecules:**

**mol_00001: C=Cc1ccc(C)o1**
```
H298_cbs = 0.15234
H298_pm7 = -3.68
abs_diff = |0.15234 - (-3.68)| = 3.83234
abs_diff_% = 3.83234 / 0.15234 * 100 = 2516.84%
nrotbonds = 1
tpsa ≈ 12.03
aromatic_rings = 1
crest_conformers_generated = 5
num_conformers_selected = 1
crest_status = "SUCCESS"
mopac_status = "SUCCESS"
```

**mol_00002: Cc1ccoc1**
```
H298_cbs = -12.64087
H298_pm7 = -14.66
abs_diff = |-12.64087 - (-14.66)| = 2.01913
abs_diff_% = 2.01913 / 12.64087 * 100 = 15.97%
nrotbonds = 0
aromatic_rings = 1
```

**mol_00003: Cc1cccc2occc12**
```
H298_cbs = -0.57999
H298_pm7 = -0.14
abs_diff = |-0.57999 - (-0.14)| = 0.43999
abs_diff_% = 0.43999 / 0.57999 * 100 = 75.86%
nrotbonds = 0
aromatic_rings = 2
```

---

## PART 2: RDKIT INTEGRATION

### 2.1 Current State

RDKit is already being used for descriptor calculation but lacks documentation.

### 2.2 Requirements

Create new documentation file: `docs/RDKIT_INTEGRATION.md`

Must follow same pattern as:
- `docs/CREST_INTEGRATION.md`
- `docs/MOPAC_INTEGRATION.md`

### 2.3 Documentation Scope

The document should include:

**Sections Required:**
1. **Overview** - What is RDKit, why use it
2. **Key Capabilities** - What descriptors it computes
3. **Phase A Integration** - How it's used in Phase A
4. **Descriptors Computed**
   - nrotbonds (Number of rotatable bonds)
   - tpsa (Topological Polar Surface Area)
   - aromatic_rings (Count of aromatic rings)
5. **Descriptor Details** - Definition, use case, typical range
6. **Workflow Diagram** - Where RDKit fits in pipeline
7. **Implementation** - Code location, key function
8. **Performance** - Execution time, scalability
9. **Limitations & Future** - Known issues, Phase B improvements
10. **References** - RDKit docs, academic papers

### 2.4 RDKit Capability Verification

**Must verify RDKit can compute:**
- ✓ nrotbonds - Yes, using `rdkit.Chem.Descriptors.NumRotatableBonds()`
- ✓ tpsa - Yes, using `rdkit.Chem.Descriptors.TPSA()`
- ✓ aromatic_rings - Yes, using ring analysis + aromaticity checks

All three descriptors are standard RDKit functionality.

---

## PART 3: CLI ENHANCEMENTS

### 3.1 Settings UI Improvements

#### Current State
```
CREST PM7 Batch Calculation Configuration

How many molecules to calculate?  3
CREST timeout per molecule (minutes)?  20
MOPAC/PM7 timeout per molecule (minutes)?  20

Configuration Summary:
  -  Molecules to calculate: 3
  -  CREST timeout: 20 min
  -  MOPAC timeout: 20 min
```

**Problems:**
- ❌ No visibility into individual tool execution
- ❌ No logs from RDKit, CREST, MOPAC
- ❌ Settings menu has boolean toggles for method selection (limited UX)

#### 3.1.1 Settings Menu Changes

**Current (Settings Menu):**
```
CREST Configuration
  Toggle GFN-FF Force Field (gfnff): [ON/OFF]
  Toggle Quick Mode (quick): [ON/OFF]
```

**New (Settings Menu):**
```
CREST Configuration
  Set CREST Method: gfn2 | gfnff | gfn2//gfnff  [→ Select]
  Set Quick Mode: off | quick | squick | mquick [→ Select]
```

**Rationale:**
- Dropdowns allow more options without multiple toggles
- Clearer semantics (method vs speed)
- Extensible for future CREST variants

**Settings Changes:**
- Remove: `gfnff: bool` toggle
- Remove: `quick: bool` toggle
- Add: `crest_method: str` ∈ {"gfn2", "gfnff", "gfn2//gfnff"}
- Add: `quick_mode: str` ∈ {"off", "quick", "squick", "mquick"}

#### 3.1.2 Logging Enhancements

**Current Batch Output:**
```
Starting batch batch_0001: 3 molecules
  Processing mol_00003 ━━━━━━━━━━━━━━━━━━ 3/3 0:24:46
✓ Batch completed!
```

**New Batch Output (With Logs):**
```
Starting batch batch_0001: 3 molecules

[mol_00001] C=Cc1ccc(C)o1
  [RDKit] Computing descriptors... ✓ nrotbonds=1, tpsa=12.03, aromatic_rings=1
  [CREST] Generating conformers (gfn2)... ✓ 5 conformers in 32s
  [MOPAC] Optimizing conformers (PM7)... ✓ 5/5 successful
    - Conformer 1: HOF = -3.68 kcal/mol
    - Conformer 2: HOF = -3.52 kcal/mol
    - Conformer 3: HOF = -3.45 kcal/mol (SELECTED)
  [Result] H298_pm7 = -3.45, abs_diff_% = 2364%
  Status: ✓ OK

[mol_00002] Cc1ccoc1
  [RDKit] Computing descriptors... ✓ nrotbonds=0, tpsa=..., aromatic_rings=1
  [CREST] Generating conformers... ✓ 1 conformer in 12s
  [MOPAC] Optimizing... ✓ 1/1 successful
  [Result] H298_pm7 = -14.66, abs_diff_% = 16%
  Status: ✓ OK

[mol_00003] Cc1cccc2occc12
  ...

Results Summary:
  -  Total: 3
  -  OK: 3
  -  Rerun: 0
  -  Skip: 0
  -  Success rate: 100.0%
```

**Logs Should Include:**
- ✓ RDKit: Descriptor values as computed
- ✓ CREST: Number of conformers, execution time
- ✓ MOPAC: HOF for each conformer, selection logic
- ✓ Result: Final H298_pm7, metrics, status

### 3.2 Help Documentation Updates

Update CLI help text to reflect new settings:

```
CREST Method (c_method):
  - gfn2: GFN2-xTB method (default, balanced speed/accuracy)
  - gfnff: GFN-FF method (faster, less accurate)
  - gfn2//gfnff: Two-step: GFN2 → GFN-FF refinement

Quick Mode (qm):
  - off: Full optimization
  - quick: Fast optimization (minimal iterations)
  - squick: Super-quick (very few iterations)
  - mquick: Meta-quick (fastest, least accurate)
```

---

## PART 4: RERUN & STATUS SYSTEM

### 4.1 Status Enum

**Valid states:**

```
Pending        → Initial state, not yet selected
Current Batch  → Selected for current batch, awaiting execution
Running        → Currently being processed
OK             → Successfully completed
Rerun          → Failed, will retry in next batch
Skip           → Failed 3 times, skipping permanently
```

### 4.2 Rerun Logic

**Rules:**
- Maximum 3 rerun attempts per molecule
- If fails: status → "Rerun", reruns += 1
- If fails with reruns ≥ 3: status → "Skip"
- In next batch: pick up "Rerun" status molecules

**CSV Tracking:**
```
mol_id  | status  | reruns | error_message
--------|---------|--------|----------------------------------
mol_X   | OK      | 0      | ""
mol_Y   | Rerun   | 1      | "MOPAC timeout after 20 min"
mol_Z   | Skip    | 3      | "Failed after 3 retries"
```

### 4.3 Column: reruns

- Type: `int`
- Range: 0-3
- Incremented on each failure
- If already 3 and fails again: status → "Skip", reruns stays at 3
- Never exceeds 3

---

## PART 5: MOPAC TIME FIX

### 5.1 Current Issue

`mopac_time` column appears empty in CSV.

### 5.2 Root Cause

Not aggregating execution time across multiple conformers.

### 5.3 Solution

In result aggregation:

```python
total_mopac_time = 0
for conformer in successful_conformers:
    total_mopac_time += conformer.mopac_execution_time

result["mopac_time"] = total_mopac_time
```

**Should be:**
- Sum of all conformer times (not average)
- In seconds
- Non-zero for successful calculations
- Present in CSV output

---

## PART 6: FILES TO MODIFY

### Priority 1 (Critical Path)

| File | Change Type | Impact | Lines | 
|------|------------|--------|-------|
| `src/grimperium/cli/settings_manager.py` | Modify | Settings UI, backward compat | ~150 |
| `[CSV Writer Location]` | Modify | Column schema, metrics | ~200 |
| `src/grimperium/crest_pm7/molecule_processor.py` | Modify | Data aggregation | ~15 |

### Priority 2 (Supporting)

| File | Change Type | Impact | Lines |
|------|------------|--------|-------|
| `src/grimperium/crest_pm7/config.py` | Verify/Modify | Status enum | ~30 |
| `docs/RDKIT_INTEGRATION.md` | Create | Documentation | ~400 |

### Priority 3 (Polish)

| File | Change Type | Impact | Lines |
|------|------------|--------|-------|
| Various | Update | Help text | ~50 |

---

## PART 7: SUCCESS CRITERIA

### Code Changes Complete When:

✓ Settings menu shows method and quick-mode as dropdowns  
✓ CSV has exactly 49 columns in specified order  
✓ All column names match specification  
✓ Backward compatibility: old CSV loads correctly  

### Data Validation When:

✓ mol_00001: nrotbonds=1, tpsa≈12, aromatic_rings=1  
✓ mol_00001: abs_diff≈3.83, abs_diff_%≈2517%  
✓ mol_00001: delta_1, delta_2, delta_3 calculated  
✓ mol_00001: conformer_selected ∈ {1,2,3}  
✓ All molecules: H298_cbs and H298_pm7 populated  
✓ All molecules: Settings columns (xtb, v3, qm, etc) populated  
✓ All molecules: mopac_time > 0  

### Testing Complete When:

✓ Settings roundtrip: dict → CSV → dict preserves values  
✓ 3-molecule batch runs successfully  
✓ Output CSV passes schema validation  
✓ All calculations verified (deltas, metrics)  
✓ No regressions in existing functionality  

### Documentation When:

✓ RDKIT_INTEGRATION.md exists  
✓ Follows same structure as CREST/MOPAC docs  
✓ Includes all required sections  
✓ Help text updated for new settings  

---

## PART 8: DATA FLOW DIAGRAM

```
┌──────────────────────────────────────────┐
│ Input: SMILES + H298_cbs (from database) │
└────────────────┬─────────────────────────┘
                 │
                 ▼
    ┌────────────────────────────────┐
    │ RDKit: Compute Descriptors     │
    │ - nrotbonds                    │
    │ - tpsa                         │
    │ - aromatic_rings               │
    └────────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ CREST: Generate Conformers │
        │ Settings: method, qm, etc  │
        │ Output: ~5 conformers      │
        └────────────────┬───────────┘
                         │
                         ▼
        ┌────────────────────────────────┐
        │ MOPAC: Optimize Each Conformer │
        │ Settings: precise_scf, etc     │
        │ Output: HOF for each           │
        └────────────────┬───────────────┘
                         │
                         ▼
        ┌────────────────────────────────┐
        │ Select Best Conformer          │
        │ delta_i = |H298_cbs - hof_i|   │
        │ Pick min(delta) conformer      │
        │ H298_pm7 = selected hof        │
        └────────────────┬───────────────┘
                         │
                         ▼
        ┌────────────────────────────────┐
        │ Calculate Metrics              │
        │ abs_diff = |H298_cbs-H298_pm7| │
        │ abs_diff_% = ...               │
        └────────────────┬───────────────┘
                         │
                         ▼
        ┌────────────────────────────────┐
        │ Write CSV (49 columns)         │
        │ Include settings used          │
        │ Include all metrics            │
        └────────────────────────────────┘
```

---

## PART 9: BACKWARD COMPATIBILITY

### Old Column Names → New Names

When loading existing CSV:

```python
Column Rename Mappings:
  "most_stable_hof" → "H298_pm7"
  "delta_e_12" → "delta_1"
  "delta_e_13" → "delta_2"
  "delta_e_15" → "delta_3"

Settings Mapping (old toggle → new dropdown):
  If "crest_gfnff" = true → "crest_method" = "gfnff"
  If "crest_gfnff" = false → "crest_method" = "gfn2"
  If "crest_quick" = true → "crest_quick_mode" = "quick"
  If "crest_quick" = false → "crest_quick_mode" = "off"

Remove columns:
  "crest_timeout"  (already have assigned_crest_timeout)
  "mopac_timeout"  (already have assigned_mopac_timeout)
```

---

## PART 10: NEXT STEPS

After Phase A completion with this spec:

### Phase B: ML Model Preparation
- Train on 29k molecules with this schema
- Build baseline: Delta-learning models
- Models: XGBoost + KRR ensemble
- Target: Predict PM7 → CBS correction

### Future Enhancements
- Add more RDKit descriptors (H-bond donors, rotatable bonds, etc)
- Molecular fingerprints (Morgan, MACCS)
- Extended Hückel theory predictions
- Automated feature selection

---

## REFERENCES

**RDKit Documentation:**
- https://www.rdkit.org/
- Descriptor module: https://www.rdkit.org/docs/GettingStartedInPython.html

**CREST Documentation:**
- https://crest-lab.github.io/crest-docs/

**MOPAC Documentation:**
- https://openmopac.net/

**Delta-Learning Concept:**
- ML model trained to predict correction: ML(X) = CBS - PM7
- More stable than direct CBS prediction
- Used in materials science, computational chemistry

---

**Version:** 1.0  
**Status:** ✅ Complete  
**Last Updated:** 2026-01-20  
**For:** Claude AI Implementation



Perfeito! Continuando com os outros documentos:

***

