# Grimperium Architecture

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

## Module Responsibilities

### 1. Data Module (`src/grimperium/data/`)

Handles all data loading, preprocessing, and fusion operations.

#### ChemperiumLoader (`loader.py`)

- Loads Chemperium dataset (CSV/Parquet)
- Validates required columns
- Provides train/test splitting
- Extracts features and targets

```python
loader = ChemperiumLoader()
df = loader.load("chemperium.csv")
train, test = loader.split(df, test_size=0.2)
```

#### DataFusion (`fusion.py`)

- Merges CBS and PM7 datasets
- Computes delta values
- Validates merge integrity

```python
fusion = DataFusion()
merged = fusion.merge(chemperium_df, pm7_df, on="smiles")
merged = fusion.compute_deltas(merged)
```

#### SemiempiricalHandler (`semiempirical.py`)

- Interface for PM7 calculations via MOPAC
- Supports CREST conformational search
- Batch processing with caching

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
                    ┌──────────────┐
                    │  Chemperium  │
                    │   Dataset    │
                    │  (52k mols)  │
                    └──────┬───────┘
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
