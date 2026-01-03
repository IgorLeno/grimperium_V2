# Feature Engineering Guide

This guide explains the molecular feature engineering approach used in Grimperium.

## Overview

Grimperium uses a **hybrid feature** approach combining three types:

```
Features = [Tabular] + [Morgan Fingerprints] + [RDKit Descriptors]
         = [3]       + [256]                + [10+]
         ≈ 270 dimensions total
```

## Feature Types

### 1. Tabular Features (3 dimensions)

Basic molecular properties from the dataset:

| Feature | Description | Range |
|---------|-------------|-------|
| `nheavy` | Number of heavy atoms (non-H) | 1-50+ |
| `charge` | Total molecular charge | -2, -1, 0, +1, +2 |
| `multiplicity` | Spin multiplicity | 1, 2, 3, ... |

**Why include?**
- Direct physical meaning
- Strong correlation with thermodynamic properties
- nheavy correlates with computational cost

### 2. Morgan Fingerprints (256 bits)

Circular fingerprints capturing molecular substructures.

**Algorithm:**
1. Assign initial identifiers to atoms
2. Iteratively update based on neighbors
3. Hash to fixed-size bit vector

**Parameters:**
- `n_bits = 256` (default, can use 512 or 1024)
- `radius = 2` (captures up to 4-bond environment)

**Example:**
```python
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles("CCO")  # Ethanol
fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=256)
```

**Why Morgan FP?**
- Captures substructure patterns
- Invariant to atom numbering
- Good for similarity searches
- Works well with tree-based models

### 3. RDKit Descriptors (10+ dimensions)

Physicochemical properties computed from SMILES:

| Descriptor | Description | Typical Range |
|------------|-------------|---------------|
| `MolWt` | Molecular weight | 16-500 g/mol |
| `TPSA` | Topological polar surface area | 0-200 Å² |
| `MolLogP` | Wildman-Crippen LogP | -5 to +10 |
| `NumRotatableBonds` | Rotatable bonds | 0-20 |
| `NumHDonors` | H-bond donors | 0-10 |
| `NumHAcceptors` | H-bond acceptors | 0-15 |
| `NumHeteroatoms` | Non-C/H atoms | 0-20 |
| `NumAromaticRings` | Aromatic rings | 0-5 |
| `FractionCSP3` | sp³ carbon fraction | 0-1 |
| `HeavyAtomMolWt` | MW without H | 12-400 |

**Example:**
```python
from rdkit import Chem
from rdkit.Chem import Descriptors

mol = Chem.MolFromSmiles("CCO")
mw = Descriptors.MolWt(mol)        # 46.07
tpsa = Descriptors.TPSA(mol)       # 20.23
logp = Descriptors.MolLogP(mol)    # -0.0014
```

## Feature Engineering Pipeline

```python
from grimperium.utils import FeatureEngineer

# Initialize
fe = FeatureEngineer(
    morgan_bits=256,
    morgan_radius=2,
    rdkit_descriptors=["MolWt", "TPSA", "MolLogP"],
    tabular_features=["nheavy", "charge", "multiplicity"],
)

# Fit and transform
X = fe.fit_transform(smiles_list, dataframe)

# Get feature names
feature_names = fe.get_feature_names()
```

## Feature Scaling

### Morgan Fingerprints
- Binary (0/1), no scaling needed
- Works well with tree-based models
- May benefit from scaling for KRR

### RDKit Descriptors
- Continuous, different scales
- StandardScaler recommended for KRR
- Optional for XGBoost

```python
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
```

## Feature Selection

### Correlation Analysis

Remove highly correlated features:

```python
import pandas as pd

corr_matrix = pd.DataFrame(X).corr()
highly_corr = (corr_matrix.abs() > 0.95).sum() > 1
```

### Feature Importance

From XGBoost:

```python
importance = model.xgb.feature_importances()
top_features = np.argsort(importance)[-20:]
```

### Recursive Feature Elimination

```python
from sklearn.feature_selection import RFE

rfe = RFE(estimator, n_features_to_select=50)
rfe.fit(X, y)
X_selected = rfe.transform(X)
```

## Alternative Features (Not Used)

### Coulomb Matrix
- 3D geometry required
- O(n²) size
- Good for quantum properties

### SOAP Descriptors
- Full 3D structure
- Very high dimensional
- Best for energies

### Graph Neural Networks
- Learn features automatically
- Requires more data
- State-of-the-art for large datasets

## Why This Feature Set?

### 1. No 3D Geometry Required

- SMILES-only features
- Faster computation
- No conformer generation needed

### 2. Balanced Information

- Structure (Morgan FP)
- Properties (RDKit)
- Size/charge (Tabular)

### 3. Computational Efficiency

- ~1000 molecules/second
- Single-threaded
- No external dependencies beyond RDKit

### 4. Interpretability

- Known descriptors
- Can analyze feature importance
- Chemically meaningful

## Performance Comparison

| Feature Set | RMSE (kcal/mol) | Training Time |
|-------------|-----------------|---------------|
| Tabular only | 2.5-3.0 | Fast |
| Morgan only | 1.5-2.0 | Fast |
| RDKit only | 1.8-2.2 | Fast |
| **Hybrid (all)** | **1.0-1.5** | Fast |
| Coulomb Matrix | 0.8-1.2 | Slow |
| SOAP | 0.5-1.0 | Very Slow |

## Best Practices

### 1. Data Preprocessing

```python
# Validate SMILES first
from grimperium.utils import validate_smiles

valid_smiles = [s for s in smiles_list if validate_smiles(s)]
```

### 2. Handle Missing Values

```python
# RDKit may fail on some molecules
X = fe.transform(smiles_list)
mask = ~np.isnan(X).any(axis=1)
X_clean = X[mask]
```

### 3. Feature Normalization

```python
# For KRR, scale features
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)
```

### 4. Dimensionality

- Start with default (256 bits Morgan)
- Increase if underfitting (512, 1024)
- Reduce if overfitting or slow

## Code Reference

### Compute Morgan Fingerprints

```python
from grimperium.utils import compute_morgan_fingerprints

fps = compute_morgan_fingerprints(
    smiles_list,
    n_bits=256,
    radius=2,
)
# Shape: (n_molecules, 256)
```

### Compute RDKit Descriptors

```python
from grimperium.utils import compute_rdkit_descriptors

descs = compute_rdkit_descriptors(
    smiles_list,
    descriptors=["MolWt", "TPSA", "MolLogP"],
)
# Shape: (n_molecules, 3)
```

### Full Pipeline

```python
from grimperium.utils import FeatureEngineer

fe = FeatureEngineer()
X = fe.fit_transform(smiles_list, df)

print(f"Feature shape: {X.shape}")
print(f"Feature names: {fe.get_feature_names()[:10]}...")
```

## Troubleshooting

### RDKit Fails on Some SMILES

```python
# Use error handling
for smiles in smiles_list:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: {smiles}")
    except Exception as e:
        print(f"Error: {e}")
```

### NaN Values in Features

```python
# Check for NaN
nan_mask = np.isnan(X).any(axis=1)
print(f"Molecules with NaN: {nan_mask.sum()}")

# Option 1: Remove
X_clean = X[~nan_mask]

# Option 2: Impute
from sklearn.impute import SimpleImputer
imputer = SimpleImputer(strategy='mean')
X_imputed = imputer.fit_transform(X)
```

### Memory Issues with Large Datasets

```python
# Use batch processing
batch_size = 10000
features = []
for i in range(0, len(smiles_list), batch_size):
    batch = smiles_list[i:i+batch_size]
    features.append(fe.transform(batch))
X = np.vstack(features)
```
