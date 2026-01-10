
# Grimperium - Project Overview

## Purpose
Grimperium is an **ML Ensemble Framework** for molecular thermodynamic property prediction using **Delta-Learning**. It combines semi-empirical quantum chemistry methods (PM7) with machine learning to predict CBS (Complete Basis Set) quality results efficiently.

**Core Hypothesis**: Learning the delta correction (y_cbs - y_pm7) is easier than learning y_cbs directly.

## Tech Stack
- **Language**: Python 3.9-3.12
- **Build System**: Poetry
- **Core Dependencies**:
  - NumPy, Pandas (data handling)
  - scikit-learn (ML base)
  - XGBoost (gradient boosting models)
  - RDKit (optional, molecular chemistry)

## Architecture

```text
grimperium/
├── src/grimperium/config.py    # Global configuration (ConfigManager)
src/grimperium/
├── core/              # Delta Learning algorithm
│   ├── delta_learning.py   # DeltaLearner class
│   └── metrics.py          # Evaluation metrics (RMSE, MAE, R², MAPE)
├── data/              # Data loading & fusion
│   ├── loader.py           # DataLoader class
│   ├── semiempirical.py    # PM6/PM7 data interface
│   └── fusion.py           # Data fusion system
├── models/            # ML models (Phase B - Available but not yet integrated)
│   ├── base_model.py       # BaseModel abstract class
│   ├── kernel_ridge.py     # KernelRidgeModel (KRR)
│   ├── xgboost_model.py    # XGBoostModel
│   └── delta_ensemble.py   # DeltaEnsemble system
├── crest_pm7/         # CREST + PM7 pipeline (Phase A - Production Ready)
│   ├── pipeline.py                # Main orchestration
│   ├── conformer_generator.py     # Conformer generation
│   ├── conformer_selector.py      # Conformer selection logic
│   ├── mopac_optimizer.py         # MOPAC optimization wrapper
│   ├── threshold_monitor.py       # Threshold monitoring system
│   ├── timeout_predictor.py       # Timeout prediction logic
│   ├── energy_extractor.py        # Energy extraction utilities
│   ├── molecule_processor.py      # Molecule processing pipeline
│   ├── result_evaluator.py        # Result evaluation logic
│   ├── config.py                  # CREST PM7 configuration
│   ├── logging_utils.py           # Logging utilities
│   └── validation.py              # Validation logic
├── utils/             # Utilities
│   ├── validation.py
│   ├── logging.py
│   └── feature_engineering.py
└── api.py             # CLI entrypoint
```

## Entry Points
- **CLI**: `grimperium` command (via poetry scripts)
- **API**: `grimperium.api:main`
- **Library**: Import modules directly

## Status
- Version: 0.2.0
- Status: Production-Ready
- Test Coverage: ~95%
