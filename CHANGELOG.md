# Changelog

All notable changes to Grimperium will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **CI/CD: Error Summary Report** (2026-01-07)
  - Fixed CI Error Summary generating contradictory status reports (PASSED ✅ but job result: failure ❌)
  - Fixed report running tools again instead of capturing original logs from failed steps
  - Implemented log capture system: each job now saves outputs to artifacts
  - Created Python script (`.github/scripts/generate_ci_report.py`) to parse captured logs
  - Report now shows exact errors from original execution with file:line context
  - Added comprehensive summary table showing status of all CI components
  - Added artifact upload for full error report (retention: 30 days)
  - Report displays in GitHub Step Summary UI with collapsible error details

- **CI/CD: Lint, Type Checks, and Module Imports** (2026-01-07)
  - Fixed lint error: Removed unused `h298_pm7` variable in `tests/integration/test_pipeline.py:270`
  - Fixed 16 mypy type errors across 5 files:
    - `src/grimperium/core/metrics.py`: Added type annotations to array conversions (9 errors)
    - `src/grimperium/core/delta_learning.py`: Added explicit `np.ndarray` type hints (1 error)
    - `src/grimperium/utils/logging.py`: Added type annotations to `__exit__` method (1 error)
    - `src/grimperium/models/delta_ensemble.py`: Added type annotations to ensemble predictions (1 error)
    - `src/grimperium/models/kernel_ridge.py`: Added `Optional[KernelRidge]` type and assertions (4 errors)
    - `src/grimperium/models/xgboost_model.py`: Fixed signature compatibility with BaseModel (6 errors)
  - Verified module structure: `grimperium.models` package imports correctly

### Breaking Changes
- Dropped support for Python 3.9 (EOL October 2025)
  - Minimum required version is now Python 3.10
  - Recommended versions: 3.10 LTS, 3.11, 3.12

### Changed
- **CI/CD: Python Version Optimization** (2026-01-07)
  - Reduced CI test matrix from 4 to 3 Python versions
  - Retained Python 3.10, 3.11, 3.12 (LTS and active versions)
  - Estimated CI time reduction: ~25%

### Added
- **BATCH 3: Hypothesis Validation Test Suite** (2026-01-07)
  - `tests/experiments/conftest.py`: New fixtures with filtered and extreme data
    - `real_data_1k_filtered()`: Realistic distribution [-1000, +1000] kcal/mol (99.1% of data)
    - `real_data_1k_extreme()`: Pathological distribution for stress testing (includes outliers)
    - `synthetic_data_1k()`: Fast synthetic data for CI/fallback tests
  - `tests/experiments/test_validate_hypothesis.py`: Main hypothesis validation tests
    - `test_decision_gate_delta_vs_direct()`: Primary test with filtered data
    - `test_synthetic_fallback()`: Synthetic data fallback test
  - `tests/experiments/test_stress_distribution_shift.py`: Robustness stress tests
    - `test_stress_distribution_shift_extreme()`: Tests with severe distribution shift
    - `test_distribution_shift_detection()`: Validates distribution shift detection
  - Comprehensive documentation of methodological decisions in test docstrings

### Changed
- **BATCH 3: Fixture Methodology Correction** (2026-01-07)
  - Replaced unfiltered data approach with filtered realistic distribution
  - Split validation testing (OPTION B) from stress testing (OPTION A)
  - Updated mock PM7 generator with configurable magnitude bias
  - Improved fixture logging and statistics reporting

### Fixed
- **BATCH 3: Distribution Shift Artifact** (2026-01-07)
  - Fixed misleading RMSE=1008 caused by severe distribution shift in unfiltered data
  - Corrected hypothesis validation to use realistic data regime (std~70, not 7230)
  - Resolved train/test distribution mismatch (6.1 vs 615 kcal/mol mean difference)
  - Fixed Direct model comparison to be fair (61.11 RMSE vs 1008.88 artifact)

### Deprecated
- **BATCH 3** (2026-01-07)
  - `tests/fixtures/conftest.py::real_data_1k()`: Now deprecated with warning
    - Reason: Uses unfiltered data causing distribution shift artifacts
    - Replacement: Use `real_data_1k_filtered` or `real_data_1k_extreme` from `tests/experiments/conftest.py`
    - Warning: Added DeprecationWarning to guide users to new fixtures

### Test Results - BATCH 3
- **Hypothesis Validation (Realistic Regime)**
  - Filter: [-1000, +1000] kcal/mol (removes 0.9% outliers)
  - Distribution shift: 6.1 kcal/mol (minimal)
  - RMSE Delta: 9.31 kcal/mol ✓
  - RMSE Direct: 61.11 kcal/mol ✓
  - Improvement: 6.6x (84.8%)
  - R² Delta: 0.9768
  - **DECISION GATE: PASS** ✅

- **Stress Test (Extreme Regime)**
  - Unfiltered data (outliers up to -325407 kcal/mol)
  - Distribution shift: 615.3 kcal/mol (severe)
  - RMSE Delta: 13.83 kcal/mol (robust)
  - RMSE Direct: 1008.88 kcal/mol (catastrophic failure expected)
  - Robustness ratio: 73x
  - **STRESS TEST: PASS** ✅

### Planned
- PM7 calculation pipeline (CREST + MOPAC integration)
- Full model implementation (fit/predict)
- Performance benchmarks on Chemperium dataset
- PyPI release

## [0.2.0] - 2024-12-29

### Added
- Complete project scaffolding
- Module structure with docstrings and type hints
- Configuration management via dataclasses
- Base model interface following scikit-learn conventions
- Test fixtures with mock Chemperium data
- Comprehensive documentation (architecture, guides)
- CI/CD pipeline with GitHub Actions
- Multi-Python testing with tox (3.9-3.12)
- Pre-commit hooks (ruff, black, mypy)

### Project Structure
```
grimperium/
├── src/grimperium/
│   ├── api.py              # High-level API (stub)
│   ├── config.py           # Configuration dataclasses
│   ├── data/
│   │   ├── loader.py       # ChemperiumLoader (stub)
│   │   ├── fusion.py       # DataFusion (stub)
│   │   └── semiempirical.py # PM7 handler (stub)
│   ├── models/
│   │   ├── base.py         # BaseModel ABC
│   │   ├── kernel_ridge.py # KRR wrapper (stub)
│   │   ├── xgboost_model.py # XGB wrapper (stub)
│   │   └── delta_ensemble.py # Ensemble (stub)
│   ├── core/
│   │   ├── delta_learning.py # DeltaLearner (stub)
│   │   └── metrics.py      # Evaluation metrics (stub)
│   └── utils/
│       ├── logging.py      # Logging utilities (stub)
│       ├── validation.py   # Input validation (stub)
│       └── feature_engineering.py # Features (stub)
├── tests/
│   ├── conftest.py         # Shared fixtures
│   ├── fixtures/
│   │   └── mock_data.py    # Mock data generators
│   ├── unit/               # Unit tests
│   └── integration/        # Integration tests
└── docs/
    ├── architecture.md     # System architecture
    ├── delta_learning_guide.md # Delta-learning concept
    └── feature_engineering.md # Feature documentation
```

### Dependencies
- numpy ^1.24.0
- pandas ^2.0.0
- scikit-learn ^1.3.0
- xgboost ^2.0.0
- rdkit ^2023.9.0

### Development Dependencies
- pytest ^7.4.0
- pytest-cov ^4.1.0
- mypy ^1.5.0
- ruff ^0.1.0
- black ^23.9.0
- pre-commit ^3.4.0
- tox ^4.11.0

## [0.1.0] - 2024-12-XX (Planned)

### Added
- Initial concept and design
- Decision documentation
- Dataset analysis (Chemperium)
- Delta-learning strategy definition
- PM7 method selection

---

## Version History Summary

| Version | Date | Description |
|---------|------|-------------|
| 0.2.0 | 2024-12-29 | Project scaffolding |
| 0.1.0 | TBD | Initial concept |
| 0.3.0 | TBD | Data loaders implementation |
| 0.4.0 | TBD | Model implementation |
| 0.5.0 | TBD | PM7 pipeline |
| 1.0.0 | TBD | Production release |
