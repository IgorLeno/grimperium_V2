# Grimperium Architecture

**Version:** 2.1
**Last Updated:** 2026-01-18
**Status:** Production (BATCH 12 - CLI fixes in progress)

This document describes the architecture of Grimperium, a delta-learning framework for molecular thermodynamic property prediction.

## System Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                         GrimperiumAPI                                │
│                    (High-level orchestration)                        │
└────────────────────────────┬────────────────────────────────────────┘
                             │
         ┌───────────────────┼───────────────────┐
         ▼                   ▼                   ▼
┌─────────────────┐ ┌─────────────────┐ ┌─────────────────┐
│      Data       │ │     Models      │ │      Core       │
│  (loader.py)    │ │  (ensemble.py)  │ │ (delta_learn)   │
│  (fusion.py)    │ │  (kernel_ridge) │ │  (metrics.py)   │
│  (semiemp.py)   │ │  (xgboost)      │ │                 │
└─────────────────┘ └─────────────────┘ └─────────────────┘
         │                   │                   │
         └───────────────────┼───────────────────┘
                             ▼
                    ┌─────────────────┐
                    │     Utils       │
                    │  (features.py)  │
                    │  (logging.py)   │
                    │  (validation)   │
                    └─────────────────┘
```

## Datasets

Grimperium uses two primary datasets for delta-learning:

### Primary: thermo_cbs_chon.csv
- **Size:** 29,568 molecules (CHON only: C, H, O, N)
- **Origin:** CBS-QB3 level thermodynamic properties
- **Removed:** Halogenes, sulfur, rare heteroatoms (52.8k → 30k → 29.6k)
- **Purpose:** High-accuracy reference for delta-learning

**Columns:**
| Column | Type | Description | Range |
|--------|------|-------------|-------|
| mol_id | string | Unique molecule identifier | - |
| smiles | string | SMILES molecular structure | - |
| charge | int | Total molecular charge | -1, 0, +1 |
| multiplicity | int | Spin multiplicity | 1-3 |
| nheavy | int | Number of heavy atoms | 1-50 |
| **H298_cbs** | float | CBS-QB3 enthalpy at 298K (kcal/mol) | -325,407 to 164,949 |
| H298_b3lyp | float | B3LYP enthalpy at 298K (kcal/mol) | - |

### Secondary: thermo_pm7.csv
- **Size:** PM7 semiempirical optimization results
- **Origin:** CREST conformer search + MOPAC PM7 optimization
- **Purpose:** Semiempirical baseline for delta-learning validation

**Columns:**
| Column | Type | Description |
|--------|------|-------------|
| mol_id | string | Unique molecule identifier |
| smiles | string | SMILES molecular structure |
| **H298_pm7** | float | PM7 enthalpy at 298K (kcal/mol) |
| conformer_details | object | Conformer generation metadata |
| quality_grade | string | Quality assessment (A, B, C) |

### Dataset Statistics

#### thermo_cbs_chon.csv (Primary)
| Metric | Value |
|--------|-------|
| Total molecules | 29,568 (CHON only) |
| Composition | C, H, O, N (no halogenes/sulfur) |
| H298_cbs mean | -320.37 kcal/mol |
| H298_cbs std | 7,230.27 kcal/mol |
| Molecular size | 1-50 heavy atoms (mean: 9.73) |
| Quality | No missing values, no duplicates |
| Creation | 52.8k → 30k → 29.6k (CHON filter) |

#### thermo_pm7.csv (Secondary)
| Metric | Value |
|--------|-------|
| Source | PM7 semiempirical optimization |
| Properties | H298_pm7, conformer data, quality grades |
| Purpose | Baseline for delta-learning validation |

📖 **Full dataset guide:** [docs/DATASETS.md](DATASETS.md)

## Module Responsibilities

### 1. Data Module (`src/grimperium/data/`)

Handles all data loading, preprocessing, and fusion operations.

#### ChemperiumLoader (`loader.py`)

Loads and validates the Grimperium thermochemistry datasets.

**Datasets Available:**

### Primary: thermo_cbs_chon.csv
- **Molecules:** 29,568 (CHON only)
- **Origin:** CBS-QB3 level thermochemistry
- **Removed:** Halogenes, sulfur, rare heteroatoms (52.8k → 30k → 29.6k)
- **Columns:** mol_id, smiles, charge, multiplicity, nheavy, H298_cbs, H298_b3lyp
- **Purpose:** High-accuracy reference for delta-learning

### Secondary: thermo_pm7.csv
- **Origin:** PM7 semiempirical optimization
- **Columns:** mol_id, smiles, H298_pm7, conformer_details, quality_grade
- **Purpose:** Semiempirical baseline predictions

**Key methods:**
- `load_thermo_cbs_chon(max_nheavy)` - Load CHON dataset (primary)
- `load_thermo_pm7(max_nheavy)` - Load PM7 results (secondary)
- `split()` - 2-way train/test split
- `train_val_test_split()` - 3-way train/val/test split
- `get_features()` / `get_targets()` - Extract X, y for ML

**Example (Primary Dataset):**
```python
from grimperium.data import ChemperiumLoader

loader = ChemperiumLoader()

# Load CHON dataset (primary)
df_chon = loader.load_thermo_cbs_chon(max_nheavy=20)

# Load PM7 results (secondary)
df_pm7 = loader.load_thermo_pm7(max_nheavy=20)

# Train/test split
train, test = loader.split(df_chon, test_size=0.2)
```

**Migration Note:**
- ❌ `load_thermo_cbs_clean()` → ✅ `load_thermo_cbs_chon()`
- ❌ `load_thermo_batch_final()` → ✅ `load_thermo_pm7()`

📖 **Full dataset guide:** [docs/DATASETS.md](DATASETS.md)

#### DataFusion (`fusion.py`)

Combines data sources and creates task-specific views for ML training.

**Key methods:**
- `merge()` - Combine Chemperium + PM7 data
- `compute_deltas()` - Calculate CBS - PM7 corrections
- `select_task_view(task)` - Get X, y for specific task:
  - `"enthalpy"` → targets H298_cbs
  - `"entropy"` → targets S298
  - `"heat_capacity"` → targets cp_1...cp_45 (multioutput)
- `analyze_deltas()` - Compute delta statistics

```python
fusion = DataFusion()
merged = fusion.merge(chemperium_df, pm7_df, on="smiles")
merged = fusion.compute_deltas(merged)
stats = fusion.analyze_deltas()  # {"mean": ..., "std": ..., ...}

# Task-specific views
X, y = fusion.select_task_view(df, task="enthalpy")
X_cp, Y_cp = fusion.select_task_view(df, task="heat_capacity")  # multioutput
```

#### SemiempiricalHandler (`semiempirical.py`)

- Interface for PM7 calculations via MOPAC
- Supports CREST conformational search
- Batch processing with caching
- **Status:** Stub (will be implemented in Batch 6)

```python
handler = SemiempiricalHandler(method="PM7")
h298_pm7 = handler.calculate_single("CCO")
```

### 2. Models Module (`src/grimperium/models/`)

Implements ML models following scikit-learn conventions.

#### BaseModel (`base.py`)

Abstract base class defining the interface:

```python
class BaseModel(ABC):
    @abstractmethod
    def fit(self, X, y) -> "BaseModel": ...

    @abstractmethod
    def predict(self, X) -> np.ndarray: ...

    def score(self, X, y) -> float: ...
```

#### KernelRidgeRegressor (`kernel_ridge.py`)

- Wraps scikit-learn's KernelRidge
- Default RBF kernel (good for smooth corrections)
- Hyperparameter tuning support

#### XGBoostRegressor (`xgboost_model.py`)

- Wraps XGBoost regressor
- Early stopping support
- Feature importance analysis

#### DeltaLearningEnsemble (`delta_ensemble.py`)

- Combines KRR + XGBoost
- Weighted averaging (default 50/50)
- Uncertainty estimation via model disagreement

```python
ensemble = DeltaLearningEnsemble(weights=(0.6, 0.4))
ensemble.fit(X_train, y_train)
predictions = ensemble.predict(X_test)
```

### 3. Core Module (`src/grimperium/core/`)

Core algorithms and orchestration.

#### DeltaLearner (`delta_learning.py`)

Main orchestrator for the delta-learning workflow:

1. Load and fuse data
2. Compute features
3. Train ensemble model
4. Evaluate performance

```python
learner = DeltaLearner()
learner.load_data(chemperium_path, pm7_path)
learner.compute_features()
learner.train()
metrics = learner.evaluate()
```

#### Metrics (`metrics.py`)

Evaluation metrics:

- RMSE (primary metric)
- MAE
- R²
- MAPE
- Method comparison utilities

### 4. Utils Module (`src/grimperium/utils/`)

Supporting utilities.

#### FeatureEngineer (`feature_engineering.py`)

Computes molecular features:

```
Features = [Tabular] + [Morgan FP] + [RDKit]
         = [3]       + [256]       + [10+]
         ≈ 270 dimensions
```

- **Tabular**: nheavy, charge, multiplicity
- **Morgan FP**: Circular fingerprints (ECFP-like)
- **RDKit**: MolWt, TPSA, LogP, NumRotatableBonds, etc.

#### Logging (`logging.py`)

- Configurable logging levels
- File and console handlers
- Progress tracking

#### Validation (`validation.py`)

- SMILES validation via RDKit
- DataFrame schema validation
- Feature array validation

## Data Flow

```
                    ┌──────────────────────────────┐
                    │      Data Layer              │
                    │  PRIMARY:   thermo_cbs_chon  │
                    │  (29,568 molecules)          │
                    │  SECONDARY: thermo_pm7       │
                    └──────────────┬───────────────┘
                           │
                    ┌──────▼───────┐
                    │   Loader     │
                    │ (validation) │
                    └──────┬───────┘
                           │
        ┌──────────────────┼──────────────────┐
        │                  │                  │
        ▼                  ▼                  ▼
┌───────────────┐  ┌───────────────┐  ┌───────────────┐
│   H298_CBS    │  │    SMILES     │  │   H298_PM7    │
│  (reference)  │  │  (features)   │  │  (semiemp)    │
└───────┬───────┘  └───────┬───────┘  └───────┬───────┘
        │                  │                  │
        │          ┌───────▼───────┐          │
        │          │    Feature    │          │
        │          │  Engineering  │          │
        │          │ Morgan+RDKit  │          │
        │          └───────┬───────┘          │
        │                  │                  │
        │          ┌───────▼───────┐          │
        │          │   Features    │          │
        │          │  (n, ~270)    │          │
        │          └───────┬───────┘          │
        │                  │                  │
┌───────▼──────────────────┼──────────────────▼───────┐
│                    DataFusion                        │
│            delta = H298_CBS - H298_PM7               │
└──────────────────────────┬──────────────────────────┘
                           │
                    ┌──────▼───────┐
                    │   Training   │
                    │ KRR + XGBoost│
                    └──────┬───────┘
                           │
                    ┌──────▼───────┐
                    │   Ensemble   │
                    │  Prediction  │
                    └──────┬───────┘
                           │
                    ┌──────▼───────┐
                    │  Evaluation  │
                    │ RMSE, MAE, R²│
                    └──────────────┘
```

## Design Decisions

### 1. Why PM7?

PM7 (Stewart, 2013) was chosen over alternatives because:

- **Best for enthalpy**: Specifically parametrized for thermochemistry
- **Good geometries**: Excellent equilibrium structures
- **Radical support**: Good spin multiplicity handling
- **Wide coverage**: C, H, O, N, S, halogens supported

### 2. Why Hybrid Features?

Combining feature types provides:

- **Tabular**: Direct physical meaning (size, charge)
- **Morgan FP**: Substructure patterns
- **RDKit**: Physicochemical properties

This captures both structural and physical information.

### 3. Why KRR + XGBoost Ensemble?

- **KRR**: Smooth, continuous corrections via RBF kernel
- **XGBoost**: Complex patterns and interactions
- **Ensemble**: Reduces variance, improves robustness

### 4. Why Delta-Learning?

Training on deltas instead of absolute values:

- Semiempirical carries ~80% of information
- ML only needs to learn the ~20% correction
- Faster convergence, better generalization
- Physically meaningful predictions

### 5. Why CHON-Only Dataset?

**Physical Reasoning:**
- Homogeneous electronic physics (no exotic valence states)
- No relativistic effects (heavy atoms)
- No hypervalence extremes
- Well-described by GFN-xTB and common DFTs

**Statistical Reasoning:**
- Rare heteroatoms act as structural outliers
- Few examples prevent generalization
- Remove without information loss

**Methodological Reasoning (Delta-Learning):**
- Delta = CBS - B3LYP correction
- Homogeneous error patterns are learnable
- Mixing regimes violates delta-learning assumptions
- Better transferability to new molecules

**Impact:**
- Original: 52,837 molecules (all heteroatoms)
- Filtered: 30,026 molecules (no halogenes/sulfur)
- Final: 29,568 molecules (no rare heteroatoms: B, P, As, Ge)
- Result: More consistent delta predictions

## Extension Points

### Adding New Models

Implement `BaseModel` interface:

```python
from grimperium.models.base import BaseModel

class MyModel(BaseModel):
    def fit(self, X, y, sample_weight=None):
        # Implementation
        self.is_fitted = True
        return self

    def predict(self, X):
        self._check_is_fitted()
        # Implementation
        return predictions
```

### Adding New Features

Extend `FeatureEngineer`:

```python
fe = FeatureEngineer(
    morgan_bits=512,  # More bits
    rdkit_descriptors=["MolWt", "TPSA", "LogP", "custom"],
)
```

### Adding New Semiempirical Methods

Extend `SemiempiricalHandler.SUPPORTED_METHODS`:

```python
handler = SemiempiricalHandler(method="PM6-D3H+")
```

## Performance Considerations

### Memory

- Features: ~270 × n_samples (float64)
- Morgan FP: Can use sparse matrices for large datasets
- RDKit: Computed on-the-fly, not stored

### Speed

- Feature computation: ~1000 molecules/second
- Training: ~10-60 seconds for 50k samples
- Prediction: ~0.1ms per molecule

### Scalability

- Supports batch processing
- Parallel feature computation (n_jobs)
- Caching for PM7 calculations
