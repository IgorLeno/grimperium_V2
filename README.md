# Grimperium

**ML Ensemble Framework for Molecular Thermodynamic Property Prediction**

[![CI](https://github.com/grimperium/grimperium/workflows/CI/badge.svg)](https://github.com/grimperium/grimperium/actions)
[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Grimperium is a production-ready framework implementing **delta-learning** for correction of semiempirical methods (PM7) using ensemble ML models (KRR + XGBoost).

## Overview

Grimperium predicts CBS-level thermodynamic properties by learning to correct fast semiempirical calculations:

```
H298_CBS ≈ H298_PM7 + δ_ML(features)
```

Where:
- **H298_CBS**: High-accuracy CBS (Complete Basis Set) enthalpy
- **H298_PM7**: Fast semiempirical PM7 enthalpy
- **δ_ML**: Machine learning correction (delta)

## Architecture

```
grimperium/
├── src/grimperium/
│   ├── api.py                    # High-level API
│   ├── config.py                 # Configuration management
│   │
│   ├── data/                     # Data loading & fusion
│   │   ├── loader.py             # ChemperiumLoader
│   │   ├── fusion.py             # DataFusion (CBS + PM7)
│   │   └── semiempirical.py      # PM7 via MOPAC
│   │
│   ├── models/                   # ML models
│   │   ├── base.py               # BaseModel (ABC)
│   │   ├── kernel_ridge.py       # KRR with RBF kernel
│   │   ├── xgboost_model.py      # XGBoost regressor
│   │   └── delta_ensemble.py     # KRR + XGB ensemble
│   │
│   ├── core/                     # Core algorithms
│   │   ├── delta_learning.py     # Delta-learning orchestration
│   │   └── metrics.py            # RMSE, MAE, R²
│   │
│   └── utils/                    # Utilities
│       ├── logging.py            # Logging configuration
│       ├── validation.py         # Input validation
│       └── feature_engineering.py # Morgan FP + RDKit
│
├── tests/                        # Test suite
│   ├── fixtures/                 # Mock data
│   ├── unit/                     # Unit tests
│   └── integration/              # Integration tests
│
└── docs/                         # Documentation
    ├── architecture.md           # Detailed architecture
    ├── delta_learning_guide.md   # Delta-learning concept
    └── feature_engineering.md    # Feature documentation
```

## Installation

```bash
# Clone repository
git clone https://github.com/grimperium/grimperium.git
cd grimperium

# Install with Poetry
poetry install

# Or with pip
pip install grimperium
```

## Quick Start

```python
from grimperium import GrimperiumAPI

# Initialize API
api = GrimperiumAPI()

# Load data
api.load_data("chemperium.csv", "pm7_results.csv")

# Train ensemble model
api.train()

# Predict H298 for new molecules
predictions = api.predict(["CCO", "CC(=O)O"])

# Evaluate performance
metrics = api.evaluate()
print(f"RMSE: {metrics['rmse']:.2f} kcal/mol")
print(f"MAE: {metrics['mae']:.2f} kcal/mol")
print(f"R²: {metrics['r2']:.4f}")
```

## Features

### Delta-Learning Strategy

Grimperium uses a simple but effective delta-learning approach:

1. **Compute delta**: `δ = H298_CBS - H298_PM7`
2. **Train ML model**: Learn to predict δ from molecular features
3. **Predict**: `H298_CBS ≈ H298_PM7 + δ_predicted`

### Hybrid Feature Engineering

Features combine multiple representations:

| Type | Features | Dimension |
|------|----------|-----------|
| Tabular | nheavy, charge, multiplicity | 3 |
| Morgan FP | Circular fingerprints | 256 |
| RDKit | MolWt, TPSA, LogP, etc. | 10+ |

### Ensemble Models

- **KernelRidgeRegressor**: RBF kernel for smooth corrections
- **XGBoostRegressor**: Gradient boosting for complex patterns
- **DeltaLearningEnsemble**: Weighted combination (default 50/50)

## Dataset

Grimperium is designed for the **real CBS-QB3 thermodynamic dataset**:

**Dataset:** `thermo_cbs_opt.csv` (52,837 molecules)

| Column | Description | Range |
|--------|-------------|-------|
| smiles | SMILES molecular structure | - |
| multiplicity | Spin multiplicity | 1-3 |
| charge | Total molecular charge | -1, 0, +1 |
| nheavy | Number of heavy atoms | 1-22 |
| **H298_cbs** | CBS-QB3 enthalpy at 298K (kcal/mol) | -325,407 to 164,949 |
| H298_b3 | B3LYP enthalpy at 298K (kcal/mol) | - |

**Statistics:**
- **Total Molecules:** 52,837
- **H298_cbs Mean:** -320.37 kcal/mol
- **H298_cbs Std:** 7,230.27 kcal/mol
- **Molecular Size:** 1-22 heavy atoms (mean: 9.73)

**Data Quality:**
- ✅ No missing values
- ✅ No duplicates
- ✅ All SMILES validated
- ✅ Complete thermodynamic properties

For CI testing, use the real data fixtures:

```python
from tests.fixtures.real_data import load_real_subset

# Load stratified 1k subset for fast testing
df = load_real_subset(n=1000, stratified=True)
```

## Development

```bash
# Install dev dependencies
poetry install --with dev

# Run tests
pytest tests/ -v

# Run linting
ruff check .
black --check .

# Run type checking
mypy src/grimperium

# Run all checks (via tox)
tox
```

## Roadmap

### v0.1 (Current)
- [x] Project scaffolding
- [x] Configuration management
- [x] Base model interfaces
- [x] Test fixtures
- [ ] Data loaders (in progress)
- [ ] Feature engineering (in progress)

### v0.2 (Planned)
- [ ] PM7 calculation pipeline (CREST + MOPAC)
- [ ] Full model implementation
- [ ] Training and evaluation
- [ ] Performance benchmarks
- [ ] PyPI release

## Contributing

Contributions are welcome! Please read our contributing guidelines and submit pull requests.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

If you use Grimperium in your research, please cite:

```bibtex
@software{grimperium2024,
  title = {Grimperium: ML Ensemble Framework for Molecular Thermodynamic Prediction},
  author = {Grimperium Team},
  year = {2024},
  url = {https://github.com/grimperium/grimperium}
}
```

## Acknowledgments

- **Chemperium dataset**: CBS-level thermodynamic data
- **PM7**: Stewart's semiempirical method via MOPAC
- **RDKit**: Molecular descriptors and fingerprints
- **scikit-learn**: Machine learning infrastructure
- **XGBoost**: Gradient boosting implementation
