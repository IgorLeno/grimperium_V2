# 📋 DOCUMENTO 3: GRIMPERIUM CONTEXT FOR CLAUDE

```markdown
# GRIMPERIUM V2 - CONTEXT FOR CLAUDE (Context7)

**Prepared for:** Claude AI via Context Window  
**Date:** 2026-01-20  
**Language:** English / Portuguese  
**Goal:** Provide complete context to implement all Phase A modifications  

---

## PROJECT OVERVIEW

**Grimperium V2** is a computational chemistry pipeline that:

1. **Generates Conformers** using CREST (GFN-xTB semiempirical)
2. **Optimizes Conformers** with MOPAC/PM7 (semiempirical method)
3. **Compares to Reference** CBS-QB3 data to build training dataset
4. **Produces Dataset** for ML models (delta-learning approach)

**Goal:** Build ML models to predict PM7 → CBS corrections for millions of molecules

### Phase A: Validation (Current)
- ✅ 3 test molecules processed
- ✅ Pipeline functions correctly
- ⚠️ Schema refinement needed
- ⚠️ Settings UI improvement needed
- ⚠️ Documentation incomplete

### Phase B: Scaling (Next)
- 📊 29,000 molecules with final schema
- 🤖 ML model training (XGBoost + KRR)
- 📈 Performance benchmarking

---

## WHAT YOU'RE IMPLEMENTING

### 1. Settings System Upgrade (30 min)
**From:** Boolean toggles for method selection  
**To:** Dropdown menus with multiple options  
**Impact:** Better UX, clearer semantics

**Example:**
- Old: `gfnff: bool` (on/off only)
- New: `crest_method: str` (gfn2, gfnff, gfn2//gfnff)

### 2. CSV Schema Refinement (45 min)
**From:** ~35 columns (incomplete)  
**To:** 49 columns (complete with metrics)  
**Impact:** Full data structure for ML training

**Key additions:**
- Metrics: `abs_diff`, `abs_diff_%`
- Deltas: `delta_1`, `delta_2`, `delta_3`
- Settings tracking: `v3`, `qm`, `c_method`, `precise_scf`, etc.
- Rerun counter: `reruns` (0-3)

### 3. Data Validation (20 min)
**Current:** Some calculations incomplete/wrong  
**Fix:** Proper aggregation and delta logic  
**Impact:** Accurate training data

**Examples:**
- `mopac_time`: Aggregate across conformers (not empty)
- `delta_*`: Select best conformer by minimum delta
- `H298_pm7`: Use best conformer's HOF

### 4. Documentation (25 min)
**Create:** RDKit integration guide  
**Template:** Follow CREST/MOPAC docs pattern  
**Impact:** Knowledge capture for Phase B

### 5. Status/Rerun System (15 min)
**Verify:** Status enum exists  
**Implement:** Retry logic (max 3 attempts)  
**Impact:** Error handling for Phase B scaling

---

## HOW THE PIPELINE WORKS

```
┌─────────────────────────────────────────┐
│ Input: SMILES (from database)           │
│ Reference: H298_cbs (CBS-QB3 from DB)   │
└──────────────────┬──────────────────────┘
                   │
         ┌─────────▼─────────┐
         │ RDKit: Descriptors│
         │ - nrotbonds       │
         │ - tpsa            │
         │ - aromatic_rings  │
         └─────────┬─────────┘
                   │ (< 1ms per molecule)
         ┌─────────▼─────────────────┐
         │ CREST: Conformers         │
         │ Settings: method, qm      │
         │ Output: ~5 conformers     │
         │ Time: ~1-2 minutes        │
         └─────────┬─────────────────┘
                   │
         ┌─────────▼──────────────────────┐
         │ MOPAC: Optimize Conformers     │
         │ Settings: precise_scf          │
         │ Output: HOF for each conformer │
         │ Time: ~1-2 minutes total       │
         └─────────┬──────────────────────┘
                   │
         ┌─────────▼──────────────────────┐
         │ Select Best Conformer          │
         │ Algorithm:                     │
         │ 1. Collect HOF values          │
         │ 2. Sort by energy              │
         │ 3. Take top 3                  │
         │ 4. Calculate delta_i =         │
         │    |H298_cbs - hof_i|          │
         │ 5. Select min(delta) conformer │
         │ 6. Use that HOF as H298_pm7    │
         └─────────┬──────────────────────┘
                   │
         ┌─────────▼──────────────────────┐
         │ Calculate Metrics              │
         │ abs_diff = |H298_cbs-H298_pm7| │
         │ abs_diff_% = (abs_diff /       │
         │              |H298_cbs|) * 100 │
         └─────────┬──────────────────────┘
                   │
         ┌─────────▼──────────────────────┐
         │ Write CSV (49 columns)         │
         │ - Molecular properties         │
         │ - Results (H298_cbs, H298_pm7) │
         │ - Metrics & deltas             │
         │ - Settings used                │
         │ - Status & tracking            │
         └────────────────────────────────┘
```

---

## CURRENT TEST DATA (3 MOLECULES)

### mol_00001: C=Cc1ccc(C)o1

```
SMILES: C=Cc1ccc(C)o1
Molecular Weight: ~136 g/mol
Structure: Vinyl furan with methyl substituent

Descriptors (RDKit):
  nrotbonds: 1 (C=C double bond can't rotate)
  tpsa: ~12.03 Ų (minimal polar surface)
  aromatic_rings: 1 (furan ring)

Reference Data:
  H298_cbs (CBS-QB3): 0.15234 kcal/mol

CREST Results:
  conformers_generated: 5
  time: ~32 seconds

MOPAC Results:
  conformers_optimized: 5 successful
  HOF values: [-3.68, -3.52, -3.45, -2.90, -2.80]
  time: ~625 seconds

Best Conformer Selection:
  Top 3 HOFs: [-3.68, -3.52, -3.45]
  delta_1 = |0.15234 - (-3.68)| = 3.83234
  delta_2 = |0.15234 - (-3.52)| = 3.67234
  delta_3 = |0.15234 - (-3.45)| = 3.60234  ← MINIMUM
  
  Selected: Conformer 3
  H298_pm7: -3.45 kcal/mol

Metrics:
  abs_diff: 3.60234 kcal/mol
  abs_diff_%: 2364.32%
```

### mol_00002: Cc1ccoc1

```
SMILES: Cc1ccoc1
Molecular Weight: ~98 g/mol
Structure: Methyl furan

Descriptors:
  nrotbonds: 0 (no rotatable bonds)
  aromatic_rings: 1

Reference Data:
  H298_cbs: -12.64087 kcal/mol

CREST Results:
  conformers_generated: 1 (very rigid)
  time: ~12 seconds

MOPAC Results:
  conformers_optimized: 1
  H298_pm7: -14.66 kcal/mol

Metrics:
  abs_diff: 2.01913 kcal/mol
  abs_diff_%: 15.97%
```

### mol_00003: Cc1cccc2occc12

```
SMILES: Cc1cccc2occc12
Molecular Weight: ~132 g/mol
Structure: Fused benzofuran system

Descriptors:
  nrotbonds: 0
  aromatic_rings: 2

Reference Data:
  H298_cbs: -0.57999 kcal/mol

Results:
  H298_pm7: -0.14 kcal/mol
  abs_diff: 0.43999 kcal/mol
  abs_diff_%: 75.86%
```

---

## UNDERSTANDING DELTA SELECTION

The algorithm to find the best conformer is **critical**:

### Wrong Way (❌ Don't do this)
```
Select conformer with MAXIMUM |delta|
This selects the conformer WORST at matching CBS
```

### Right Way (✅ Do this)
```
Select conformer with MINIMUM |delta|
This selects the conformer BEST at matching CBS
```

### Example Walkthrough

```
H298_cbs = 0.15234 (reference)
HOF values (from MOPAC):
  Conformer 1: -3.68 → delta_1 = |0.15234 - (-3.68)| = 3.83234
  Conformer 2: -3.52 → delta_2 = |0.15234 - (-3.52)| = 3.67234
  Conformer 3: -3.45 → delta_3 = |0.15234 - (-3.45)| = 3.60234 ← MIN
  Conformer 4: -2.90 → delta_4 = |0.15234 - (-2.90)| = 3.05234
  Conformer 5: -2.80 → delta_5 = |0.15234 - (-2.80)| = 2.95234

Top 3 deltas: [3.83234, 3.67234, 3.60234]

Best conformer = one with minimum delta = Conformer 3

H298_pm7 = -3.45 (use Conformer 3's HOF)

Store in CSV:
  delta_1 = 3.83234
  delta_2 = 3.67234
  delta_3 = 3.60234
  conformer_selected = 3
  H298_pm7 = -3.45
```

---

## KEY DECISION POINTS

### 1. CSV Column Order Matters
User verified exact sequence → Don't rearrange  
Must be precisely: mol_id, status, smiles, ... (49 total)

### 2. Backward Compatibility is Required
Old CSV with `crest_gfnff=true` must load  
Must auto-convert to `crest_method=gfnff`

### 3. Settings Persistence is New
Every molecule row includes settings used  
Allows: "This molecule was calculated with these settings"  
Enables: Historical tracking across batches

### 4. RDKit Already Works
Just needs documentation  
Verify: Can it compute nrotbonds, tpsa, aromatic_rings?  
Answer: Yes, all are standard RDKit functions

### 5. No Major Refactoring
Just surgical fixes  
Don't improve or rewrite unrelated code  
Just fix what's broken

---

## WHAT'S ALREADY WORKING

✅ **CREST Integration**
- Generates conformers
- Settings system basic
- Execution tracking

✅ **MOPAC Integration**
- Optimizes conformers
- Extracts HOF values
- Timeout handling

✅ **RDKit Integration**
- Descriptor computation
- SMILES parsing
- 3 descriptors working

✅ **Database Integration**
- Loads H298_cbs
- Stores results
- CSV output basic

✅ **Timeout System**
- Predicts times
- Allocates resources
- Tracks execution

---

## WHAT NEEDS FIXING

❌ **Settings UI**
- Toggles don't scale
- Need dropdowns

❌ **CSV Schema**
- Missing columns (metrics, settings)
- Wrong names (most_stable_hof)
- Wrong deltas (unclear naming)
- Old columns duplicate

❌ **Data Validation**
- mopac_time empty
- Deltas unclear
- No conformer selection tracking

❌ **Documentation**
- RDKit not documented
- Help text outdated

❌ **Status System**
- Rerun logic unclear
- Need counter

---

## SUCCESS = THIS CSV

After implementation, thermo_pm7.csv will have:

```
mol_id,status,smiles,multiplicity,charge,nheavy,H298_cbs,H298_pm7,abs_diff,abs_diff_%,batch_id,timestamp,reruns,nrotbonds,tpsa,aromatic_rings,crest_status,xtb,v3,qm,nci,c_method,energy_window,rmsd_threshold,threads,crest_conformers_generated,crest_time,num_conformers_selected,mopac_status,precise_scf,scf_threshold,mopac_time,delta_1,delta_2,delta_3,conformer_selected,error_message,batch_order,batch_failure_policy,assigned_crest_timeout,assigned_mopac_timeout

mol_00001,OK,C=Cc1ccc(C)o1,1,0,8,0.15234,-3.45,3.60234,2364.32,batch_0001,20/01-01:31,0,1,12.03,1,SUCCESS,true,true,off,false,gfn2,5.0,1.5,4,5,32,1,SUCCESS,false,1e-4,625.0,3.83234,3.67234,3.60234,3,,1.0,partial_ok,20.0,20.0

mol_00002,OK,Cc1ccoc1,1,0,6,-12.64087,-14.66,2.01913,15.97,batch_0001,20/01-01:42,0,0,X,1,SUCCESS,true,true,off,false,gfn2,5.0,1.5,4,1,12,1,SUCCESS,false,1e-4,248.8,2.01913,,,1,,2.0,partial_ok,20.0,20.0

mol_00003,OK,Cc1cccc2occc12,1,0,10,-0.57999,-0.14,0.43999,75.86,batch_0001,20/01-01:46,0,0,Y,2,SUCCESS,true,true,off,false,gfn2,5.0,1.5,4,1,50,1,SUCCESS,false,1e-4,597.7,0.43999,,,1,,3.0,partial_ok,20.0,20.0
```

**49 columns total**, all populated, ready for ML training.

---

## FILE CHANGES SUMMARY

| File | Type | Impact | Time |
|------|------|--------|------|
| settings_manager.py | Modify | Settings UI | 30 min |
| [CSV writer] | Modify | Schema | 45 min |
| molecule_processor.py | Modify | Data fix | 20 min |
| config.py | Verify | Status enum | 10 min |
| RDKIT_INTEGRATION.md | Create | Documentation | 25 min |
| Help text | Update | UX | 10 min |

**Total: ~2.5 hours development + 1 hour testing**

---

## TROUBLESHOOTING CHECKLIST

### "Where's the CSV writer?"
```bash
grep -r "\.to_csv\|thermo_pm7" src/ --include="*.py"
grep -r "DataFrame" src/grimperium/crest_pm7/ --include="*.py"
```

### "What's the delta formula?"
```
delta = |H298_cbs - hof_conformer|
Select: Minimum delta
Use: That conformer's HOF as H298_pm7
```

### "Should I refactor?"
No. Surgical fixes only. Keep existing patterns.

### "Do I need tests?"
Yes. Use grimperium_final_checklist.md Testing section.

### "What about old CSVs?"
Handle in from_dict() with backward compatibility logic.

---

## NEXT STEPS AFTER PHASE A

### Phase B Starts When:
✅ 3-molecule batch validated  
✅ CSV schema finalized  
✅ Settings improved  
✅ All tests pass  

### Phase B Goals:
- Scale to 29,000 molecules
- Build ML models (XGBoost + KRR)
- Predict PM7 → CBS corrections
- Achieve <5% mean absolute error

### Phase B Inputs:
- Finalized CSV (49 columns)
- Settings persistence
- Rerun/failure handling
- Performance benchmarks

---

**Status:** ✅ Complete Context Package  
**Ready for:** Claude Implementation  
**Timeline:** ~4-6 hours total  
```
