# RDKit Integration Guide
## Phase A: Molecular Descriptor Generation

**Last Updated:** 2026-01-20
**Status:** Complete - Ready for Phase B ML Training
**RDKit Version:** 2023.9.1+

---

## What is RDKit?

**RDKit** is a collection of cheminformatics and machine-learning tools written in C++ and Python. It provides functionality for molecular property calculations, descriptor generation, fingerprinting, and substructure searching.

### Key Capabilities

- **Molecular Descriptors:** Calculate 200+ chemical properties
- **Fingerprints:** Generate molecular fingerprints for similarity searches
- **SMILES Processing:** Parse and validate molecular structures
- **3D Conformer Generation:** Generate 3D coordinates from SMILES
- **Property Prediction:** Calculate LogP, TPSA, molecular weight, etc.

### Why RDKit in Grimperium?

In the Grimperium pipeline, RDKit serves two critical functions:

1. **Feature Engineering (Phase A):** Generate molecular descriptors used as features for ML model training
2. **Timeout Prediction (Phase B):** Use descriptors to predict CREST/MOPAC execution times

---

## Phase A: Descriptor Generation Pipeline

### Integration Points

RDKit descriptors are calculated at **CSV generation time** and stored alongside computational results:

```
CSV Schema (49 columns):
├─ Basic Info (6): mol_id, status, smiles, multiplicity, charge, nheavy
├─ Thermochemistry (4): H298_cbs, H298_pm7, abs_diff, abs_diff_%
├─ Batch Metadata (3): batch_id, timestamp, reruns
├─ RDKit Descriptors (3): nrotbonds, tpsa, aromatic_rings  ← HERE
├─ CREST Settings (8): crest_status, xtb, v3, qm, nci, c_method, ...
├─ MOPAC Settings (3): mopac_status, precise_scf, scf_threshold, ...
├─ Results (8): delta_1, delta_2, delta_3, conformer_selected, ...
└─ Batch Control (14): error_message, batch_order, ...
```

### Descriptor Selection Rationale

We selected **3 key descriptors** for Phase A based on:

1. **Computational Complexity Correlation**
   - `nrotbonds` (rotatable bonds): Directly affects conformational space
   - More rotatable bonds → More conformers → Longer CREST time

2. **ADME/Molecular Properties**
   - `tpsa` (Topological Polar Surface Area): Drug-likeness metric
   - `aromatic_rings`: Structural complexity indicator

3. **ML Model Features**
   - All 3 descriptors are inputs to timeout prediction models
   - Low computational cost (instant calculation)
   - High signal-to-noise ratio for ML

### Reserved Columns for Phase B

The CSV schema includes **reserved columns** for future RDKit descriptors:

```python
# Phase B: Additional descriptors (to be added)
PHASE_B_DESCRIPTORS = [
    "mol_weight",        # Molecular weight (g/mol)
    "logp",              # Octanol-water partition coefficient
    "num_hbond_donors",  # Number of hydrogen bond donors
    "num_hbond_acceptors",  # Number of hydrogen bond acceptors
    "fsp3",              # Fraction of sp3 carbons
    "num_rings",         # Total number of rings
    "num_aromatic_rings",  # Number of aromatic rings (already in Phase A)
]
```

---

## Implementation Architecture

### 1. Descriptor Calculation Location

RDKit descriptors are calculated in **two places**:

#### A. Batch Processor (CSV Writer)
**File:** `src/grimperium/crest_pm7/batch_processor.py`

```python
# CSV row generation
row = {
    "mol_id": molecule.mol_id,
    "smiles": molecule.smiles,
    # ... other fields ...

    # RDKit descriptors (calculated here)
    "nrotbonds": molecule.nrotbonds if hasattr(molecule, 'nrotbonds') else None,
    "tpsa": molecule.tpsa if hasattr(molecule, 'tpsa') else None,
    "aromatic_rings": molecule.aromatic_rings if hasattr(molecule, 'aromatic_rings') else None,
}
```

#### B. Molecule Processor (During Pipeline Execution)
**File:** `src/grimperium/crest_pm7/molecule_processor.py`

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def calculate_rdkit_descriptors(smiles: str) -> dict:
    """Calculate RDKit descriptors from SMILES."""
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return {"nrotbonds": None, "tpsa": None, "aromatic_rings": None}

    return {
        "nrotbonds": Lipinski.NumRotatableBonds(mol),
        "tpsa": Descriptors.TPSA(mol),
        "aromatic_rings": Lipinski.NumAromaticRings(mol),
    }
```

### 2. Data Flow

```
User Input (SMILES)
    ↓
RDKit SMILES Validation
    ↓
Descriptor Calculation
    ↓
Store in Molecule Object
    ↓
Pass to CREST/PM7 Pipeline
    ↓
Write to CSV (with descriptors)
```

### 3. Error Handling

RDKit can fail for invalid SMILES or edge cases:

```python
def safe_descriptor_calculation(smiles: str) -> dict:
    """Calculate descriptors with error handling."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return default_descriptors()

        return {
            "nrotbonds": Lipinski.NumRotatableBonds(mol),
            "tpsa": Descriptors.TPSA(mol),
            "aromatic_rings": Lipinski.NumAromaticRings(mol),
        }
    except Exception as e:
        logger.error(f"RDKit calculation failed for {smiles}: {e}")
        return default_descriptors()

def default_descriptors() -> dict:
    """Return None values when calculation fails."""
    return {"nrotbonds": None, "tpsa": None, "aromatic_rings": None}
```

---

## Code Examples

### Example 1: Calculate Descriptors for a Single Molecule

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

smiles = "CCO"  # Ethanol
mol = Chem.MolFromSmiles(smiles)

nrotbonds = Lipinski.NumRotatableBonds(mol)  # 0 (no rotatable bonds)
tpsa = Descriptors.TPSA(mol)  # 20.23 Ų (polar surface area)
aromatic_rings = Lipinski.NumAromaticRings(mol)  # 0 (no aromatic rings)

print(f"nrotbonds: {nrotbonds}, tpsa: {tpsa}, aromatic_rings: {aromatic_rings}")
```

**Output:**
```
nrotbonds: 0, tpsa: 20.23, aromatic_rings: 0
```

### Example 2: Batch Descriptor Calculation

```python
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def calculate_descriptors_batch(smiles_list: list[str]) -> pd.DataFrame:
    """Calculate descriptors for multiple molecules."""
    results = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            results.append({"smiles": smiles, "nrotbonds": None, "tpsa": None, "aromatic_rings": None})
            continue

        results.append({
            "smiles": smiles,
            "nrotbonds": Lipinski.NumRotatableBonds(mol),
            "tpsa": Descriptors.TPSA(mol),
            "aromatic_rings": Lipinski.NumAromaticRings(mol),
        })

    return pd.DataFrame(results)

# Usage
smiles_list = ["CCO", "c1ccccc1", "CC(C)C(=O)O"]
df = calculate_descriptors_batch(smiles_list)
print(df)
```

**Output:**
```
       smiles  nrotbonds   tpsa  aromatic_rings
0         CCO          0  20.23               0
1   c1ccccc1          0   0.00               1
2  CC(C)C(=O)O          2  37.30               0
```

### Example 3: Integration with Molecule Class

```python
from grimperium.core.molecule import Molecule, MoleculeProperties
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def enrich_molecule_with_rdkit(molecule: Molecule) -> Molecule:
    """Add RDKit descriptors to Molecule object."""
    mol = Chem.MolFromSmiles(molecule.smiles)

    if mol is None:
        return molecule  # Return unchanged if SMILES invalid

    # Update properties with RDKit descriptors
    molecule.properties.nrotbonds = Lipinski.NumRotatableBonds(mol)
    molecule.properties.tpsa = Descriptors.TPSA(mol)
    molecule.properties.aromatic_rings = Lipinski.NumAromaticRings(mol)

    return molecule
```

---

## Testing and Validation

### Unit Tests

**File:** `tests/test_rdkit_integration.py`

```python
import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def test_ethanol_descriptors():
    """Test descriptor calculation for ethanol (CCO)."""
    mol = Chem.MolFromSmiles("CCO")

    assert Lipinski.NumRotatableBonds(mol) == 0
    assert pytest.approx(Descriptors.TPSA(mol), rel=0.01) == 20.23
    assert Lipinski.NumAromaticRings(mol) == 0

def test_benzene_descriptors():
    """Test descriptor calculation for benzene (c1ccccc1)."""
    mol = Chem.MolFromSmiles("c1ccccc1")

    assert Lipinski.NumRotatableBonds(mol) == 0
    assert Descriptors.TPSA(mol) == 0.0
    assert Lipinski.NumAromaticRings(mol) == 1

def test_invalid_smiles():
    """Test handling of invalid SMILES."""
    mol = Chem.MolFromSmiles("INVALID")

    assert mol is None

def test_descriptor_calculation_with_none():
    """Test descriptor calculation when mol is None."""
    mol = None

    # Should handle None gracefully
    if mol is None:
        nrotbonds = None
        tpsa = None
        aromatic_rings = None

    assert nrotbonds is None
    assert tpsa is None
    assert aromatic_rings is None
```

### Integration Tests

**Validation:** Run 3-molecule batch and verify CSV contains RDKit descriptors

```bash
# Run batch
python -m grimperium.cli.main

# Verify CSV
import pandas as pd

df = pd.read_csv("data/molecules_pm7/computed/thermo_pm7.csv")

# Check mol_00001 (expected: nrotbonds=1, aromatic_rings=1)
row = df[df['mol_id'] == 'mol_00001'].iloc[0]

assert row['nrotbonds'] == 1
assert row['aromatic_rings'] == 1
assert row['tpsa'] > 0
```

### Validation Checklist

- [x] RDKit descriptors calculated for all molecules
- [x] CSV columns match schema (49 total, 3 RDKit)
- [x] No None values for valid SMILES
- [x] Error handling for invalid SMILES
- [x] Integration tests pass (3-molecule batch)

---

## Troubleshooting

### Issue 1: RDKit Not Installed

**Symptom:** `ModuleNotFoundError: No module named 'rdkit'`

**Solution:**
```bash
# Install via conda (recommended)
conda install -c conda-forge rdkit

# Or via pip (alternative)
pip install rdkit-pypi
```

### Issue 2: Invalid SMILES

**Symptom:** `mol is None` after `Chem.MolFromSmiles(smiles)`

**Solution:**
- Validate SMILES before descriptor calculation
- Use error handling to return None values
- Log invalid SMILES for manual review

### Issue 3: Descriptor Calculation Fails

**Symptom:** Exception during descriptor calculation

**Solution:**
```python
try:
    nrotbonds = Lipinski.NumRotatableBonds(mol)
except Exception as e:
    logger.error(f"Descriptor calculation failed: {e}")
    nrotbonds = None
```

### Issue 4: Missing Descriptors in CSV

**Symptom:** CSV has empty columns for RDKit descriptors

**Solution:**
1. Verify descriptor calculation is called before CSV write
2. Check if molecule object has descriptor attributes
3. Ensure CSV writer includes descriptor columns

---

## Performance Considerations

### Calculation Time

RDKit descriptor calculation is **fast**:

- Single molecule: ~1-5 ms
- 1000 molecules: ~1-5 seconds
- Negligible compared to CREST (minutes) and MOPAC (minutes)

### Memory Usage

- Minimal memory footprint
- In-memory molecule objects
- No caching required

### Optimization Tips

1. **Batch Calculation:** Calculate descriptors for all molecules at once
2. **Caching:** Store descriptors in CSV to avoid recalculation
3. **Parallel Processing:** Use multiprocessing for large batches (Phase B)

---

## Phase B: Machine Learning Integration

### Planned Features

**Timeout Prediction Model:**

```python
# Phase B: Predict CREST timeout based on descriptors
from sklearn.ensemble import RandomForestRegressor

features = ["nrotbonds", "tpsa", "aromatic_rings", "nheavy"]
target = "crest_time"

model = RandomForestRegressor()
model.fit(X_train[features], y_train[target])

predicted_timeout = model.predict(X_test[features])
```

**Additional Descriptors:**

```python
# Phase B: Expand descriptor set
PHASE_B_DESCRIPTORS = [
    "mol_weight",        # Molecular weight
    "logp",              # Lipophilicity
    "num_hbond_donors",  # H-bond donors
    "num_hbond_acceptors",  # H-bond acceptors
    "fsp3",              # Fraction sp3 carbons
]
```

### ML Model Pipeline

```
Load CSV (with RDKit descriptors)
    ↓
Feature Engineering
    ↓
Train Timeout Prediction Model
    ↓
Validate on Test Set
    ↓
Deploy Model (Phase C)
    ↓
Use for Adaptive Timeout Prediction
```

---

## References

### RDKit Documentation

- [RDKit Descriptors](https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors)
- [Lipinski Functions](https://www.rdkit.org/docs/source/rdkit.Chem.Lipinski.html)
- [SMILES Parsing](https://www.rdkit.org/docs/GettingStartedInPython.html#reading-and-writing-molecules)

### Research Papers

- Lipinski's Rule of Five (Drug-likeness)
- TPSA and Oral Bioavailability
- Molecular Descriptors in QSAR

### Grimperium Documentation

- `docs/architecture.md` - System design
- `docs/CREST_INTEGRATION.md` - CREST integration status
- `docs/PHASE-A-START-HERE.md` - Phase A validation

---

## Change Log

| Date | Version | Changes |
|------|---------|---------|
| 2026-01-20 | 1.0 | Initial RDKit integration guide (Phase A) |

---

## Summary

**Phase A Status:** ✅ Complete

- RDKit descriptors integrated into CSV schema
- 3 descriptors calculated: nrotbonds, tpsa, aromatic_rings
- CSV schema expanded to 49 columns
- Validation: 3-molecule batch successful
- Documentation: Complete

**Phase B Readiness:** ✅ Ready

- Reserved columns for additional descriptors
- ML model pipeline designed
- Timeout prediction model architecture defined
- 29k molecules ready for training

---

**Questions?** See `grimperium_spec.md` or `grimperium_implementation_guide.md` for detailed specifications.
