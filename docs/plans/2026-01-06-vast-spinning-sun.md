# Critical Analysis: PM7 + Delta-Learning Pipeline Plan

## Executive Summary

The original plan has **significant structural issues** that will lead to wasted effort. The core problem: **you're building UI and infrastructure before validating the core hypothesis**.

---

## Critical Issues Identified

### Issue 1: Inverted Dependency Chain (CRITICAL)

**Problem**: Batch 3 (Models ML) is scheduled before Batch 4 (H298_PM7 Calculation), but:
- Delta-learning formula: `H298_pred = H298_PM7 + delta_ML`
- Models need `delta = H298_CBS - H298_PM7` to train
- Without H298_PM7, models can only train on mocks

**Impact**: You'll build models that work on fake data, then potentially need to redesign when real PM7 data has different error distributions.

**Solution**: Use mock PM7 with **realistic error patterns** OR acquire a small subset of real PM7 data first.

---

### Issue 2: Over-Engineering Before Validation (HIGH)

**Problem**: The plan creates:
- `calculations/` module (5 files)
- `cli/` module (4 files)
- `dashboard/` module (4+ files)
- `scripts/` directory
- Model Registry system

...before proving that delta-learning actually improves predictions.

**Impact**: 3-4 days building infrastructure that may need redesign.

**Solution**: **Prove the hypothesis first** with minimal code:
1. Generate mock PM7 (realistic noise: ~5-15 kcal/mol from CBS)
2. Train KRR + XGBoost on delta
3. Compare RMSE: direct ML vs delta-learning
4. If delta-learning wins, then build infrastructure

---

### Issue 3: CREST/MOPAC Bottleneck Ignored (HIGH)

**Problem**: 52,837 molecules × 57s = **~35 days on single core**.
Even with 16 cores: **~3 days** just for conformer search.

The plan acknowledges this but doesn't address it strategically.

**Solution**:
1. Start with subset (1k-5k molecules) for development
2. Use pre-computed PM7 if available (check literature/databases)
3. Consider GFN2-xTB as faster alternative to full CREST+MOPAC
4. Run GridUNESP job in parallel with development work

---

### Issue 4: Interface Before Function (MEDIUM)

**Problem**: Batch 6 (CLI/Dashboard) is planned before core ML works.

**Impact**: You'll build progress bars for processes that don't exist yet.

**Solution**: CLI/Dashboard should be the **last** batch, after:
1. ML models work
2. Metrics are validated
3. End-to-end pipeline tested

---

### Issue 5: Missing Baseline Comparison (MEDIUM)

**Problem**: Plan doesn't establish clear baseline to beat.

**Solution**: Before any ML implementation:
1. Calculate: What's the RMSE of H298_PM7 vs H298_CBS? (This is your baseline)
2. Target: Delta-ML must beat this baseline significantly (>20% improvement)
3. Document this in the plan

---

## Recommended Batch Reordering

### Phase 0: Validate Hypothesis (NEW - 2h)
```
Goal: Prove delta-learning works before building infrastructure

1. Mock PM7 Generator:
   - H298_PM7_mock = H298_CBS + N(0, sigma) where sigma ∈ [5, 15] kcal/mol
   - Add systematic bias for heavy atoms (PM7 struggles with conjugated systems)

2. Quick Prototype:
   - Load real CBS data (52k molecules)
   - Generate mock PM7
   - Train simple KRR on delta
   - Compare RMSE: KRR(direct) vs PM7 + KRR(delta)

3. Decision Gate:
   - If delta-learning improves RMSE by <10%: reconsider approach
   - If improvement >20%: proceed with full implementation

Files: tests/experiments/validate_delta_learning.py (disposable)
```

### Phase 1: Core ML (Current Batch 3, but enhanced)
```
Files to implement:
- src/grimperium/models/kernel_ridge.py     # Real implementation
- src/grimperium/models/xgboost_model.py    # Real implementation
- src/grimperium/models/delta_ensemble.py   # Real implementation
- src/grimperium/core/metrics.py            # RMSE, MAE, R²
- src/grimperium/core/delta_learning.py     # DeltaLearner orchestration

Test with: Mock PM7 data (realistic noise patterns)
```

### Phase 2: Feature Engineering
```
Files to implement:
- src/grimperium/utils/feature_engineering.py  # Morgan FP, RDKit descriptors

Note: Currently scaffolded. Implement actual feature extraction.
Dependencies: rdkit
```

### Phase 3: API Layer + Model Registry
```
Files to implement:
- src/grimperium/api.py                    # Complete implementation
- src/grimperium/models/model_registry.py  # NEW - save/load/version

Note: Merge Model Registry into models/ rather than creating separate module
```

### Phase 4: PM7 Calculation Wrappers
```
Files to implement:
- src/grimperium/data/semiempirical.py     # Complete the existing stub
- src/grimperium/calculations/            # NEW if needed

Strategy:
1. Start with RDKit → XYZ (already have rdkit dependency)
2. MOPAC wrapper (simplest: just PM7, no CREST conformer search initially)
3. Add CREST later only if conformer search proves necessary

Simplification: Skip CREST initially. Many molecules have only 1-2 conformers.
```

### Phase 5: CLI Interface (Only after core works)
```
Recommendation: Rich CLI only (no Streamlit initially)

Files:
- src/grimperium/cli/__init__.py
- src/grimperium/cli/main.py

Keep minimal:
- grimperium train --data file.csv
- grimperium predict --model model.joblib --smiles "CCO"
- grimperium evaluate --model model.joblib --data test.csv
```

### Phase 6: Production Calculations
```
Only after everything else works:
- GridUNESP scripts
- Full 52k molecule processing
- Dashboard (if still needed)
```

---

## Simplified Architecture

Instead of 4 new directories, consolidate:

```
src/grimperium/
├── api.py                    # Complete orchestration
├── config.py                 # Already done
├── data/
│   ├── loader.py            # Done
│   ├── fusion.py            # Done
│   └── semiempirical.py     # Complete PM7 handling here
├── models/
│   ├── base.py              # Complete ABC
│   ├── kernel_ridge.py      # Implement
│   ├── xgboost_model.py     # Implement
│   ├── delta_ensemble.py    # Implement
│   └── registry.py          # NEW - save/load/list
├── core/
│   ├── delta_learning.py    # Implement
│   └── metrics.py           # Implement
├── utils/
│   └── feature_engineering.py  # Implement
└── cli/                     # NEW (minimal)
    └── main.py
```

**What we're NOT creating:**
- ❌ `calculations/` directory (use `data/semiempirical.py`)
- ❌ `dashboard/` directory (premature)
- ❌ `scripts/` directory (put in project root or later)
- ❌ Multiple CLI submodules (one file is enough initially)

---

## Questions for Clarification

Before finalizing the plan, I need to understand:

1. **Do you have pre-computed PM7 data?**
   - If yes: Skip Phase 4 (PM7 wrappers), use existing data
   - If no: Phase 4 becomes critical path

2. **What's your timeline constraint?**
   - 1 week: Focus on mock validation + core ML only
   - 2-3 weeks: Can include PM7 wrappers + CLI
   - 1 month+: Full plan including real calculations

3. **Is CREST conformer search strictly necessary?**
   - For initial validation: Probably not
   - For publication-quality: Yes
   - Recommendation: Skip CREST initially, add later

4. **Interface priority?**
   - A: Python API only (fastest, for experienced users)
   - B: Minimal CLI (grimperium train/predict/evaluate)
   - C: Rich CLI with progress bars
   - D: Web dashboard (Streamlit)
   - Recommendation: Start with A+B, add C/D later

---

## Recommended Immediate Next Steps

1. **Phase 0 validation** (2h): Prove delta-learning works with mock PM7
2. **Implement core ML models** (4-6h): KRR + XGBoost + Ensemble
3. **Implement metrics** (1h): RMSE, MAE, R²
4. **Integration test** (2h): End-to-end with real CBS data + mock PM7
5. **Decision point**: If it works, continue; if not, investigate why

---

## Files to Modify (Immediate)

| File | Action | Priority |
|------|--------|----------|
| `src/grimperium/models/kernel_ridge.py` | Implement fit/predict | P0 |
| `src/grimperium/models/xgboost_model.py` | Implement fit/predict | P0 |
| `src/grimperium/models/delta_ensemble.py` | Implement orchestration | P0 |
| `src/grimperium/core/metrics.py` | Implement RMSE/MAE/R² | P0 |
| `src/grimperium/core/delta_learning.py` | Implement DeltaLearner | P1 |
| `tests/fixtures/conftest.py` | Add realistic mock PM7 generator | P1 |
| `tests/integration/test_delta_learning.py` | End-to-end test | P1 |

---

## Final Plan (Based on User Answers)

**Constraints:**
- No PM7 data available (must calculate)
- Timeline: 2-3 weeks
- Interface: Rich CLI

---

## WEEK 1: Core ML + Validation (Days 1-5)

### Day 1-2: Mock PM7 + Validate Hypothesis
```
Goal: Prove delta-learning works before investing in infrastructure

Create: tests/experiments/validate_delta_learning.py
- Mock PM7 generator with realistic error patterns
- Quick prototype: KRR on direct vs delta
- Document baseline RMSE to beat

Modify: tests/fixtures/conftest.py
- Add create_mock_pm7() with realistic noise (5-15 kcal/mol std)
- Add systematic bias for heavy atoms
```

### Day 2-3: Implement Core Models
```
Modify: src/grimperium/models/kernel_ridge.py
- Implement fit() using sklearn.kernel_ridge.KernelRidge
- Implement predict()
- Add hyperparameter tuning (GridSearchCV)

Modify: src/grimperium/models/xgboost_model.py
- Implement fit() using xgboost.XGBRegressor
- Implement predict()
- Add feature_importances()

Modify: src/grimperium/models/delta_ensemble.py
- Implement ensemble logic (weighted average)
- tune_weights() via optimization
```

### Day 3-4: Metrics + Delta Learning
```
Modify: src/grimperium/core/metrics.py
- RMSE, MAE, R², MAPE, max_error (all 5 functions)
- compute_all_metrics() returning dict
- compare_methods() for baseline comparison

Modify: src/grimperium/core/delta_learning.py
- DeltaLearner.fit() orchestrating ensemble
- DeltaLearner.predict()
- DeltaLearner.evaluate() with metrics
```

### Day 4-5: Integration Tests
```
Create: tests/integration/test_delta_learning.py
- End-to-end test with real CBS + mock PM7
- Compare: direct ML vs delta-learning
- K-fold cross-validation test

Validate: All 95+ tests pass, >79% coverage
```

---

## WEEK 2: PM7 Wrappers + CLI (Days 6-10)

### Day 6-7: PM7 Calculation (Local)
```
Modify: src/grimperium/data/semiempirical.py
- RDKit SMILES → XYZ conversion
- MOPAC input file generation (.mop)
- MOPAC output parsing (extract H298)
- calculate_single() real implementation
- calculate_batch() with ProcessPoolExecutor

Dependencies: Add ase or directly call MOPAC subprocess

Test with: 100-molecule subset locally
```

### Day 8: Rich CLI
```
Create: src/grimperium/cli/__init__.py
Create: src/grimperium/cli/main.py

Commands:
- grimperium calculate-h298 --input smiles.csv --output h298.csv --parallel 4
- grimperium train --data merged.csv --model-type ensemble --output models/
- grimperium predict --model models/ensemble.joblib --smiles "CCO"
- grimperium evaluate --model models/ensemble.joblib --data test.csv

Rich features:
- Progress bars for batch operations
- Colored output (success/warning/error)
- Tables for results display

Update: pyproject.toml [tool.poetry.scripts]
Add dependency: rich
```

### Day 9-10: Model Registry + API
```
Create: src/grimperium/models/registry.py
- save_model() with metadata (RMSE, R², timestamp)
- load_model() by name or latest
- list_models() returning formatted table

Modify: src/grimperium/api.py
- Complete all NotImplementedError stubs
- load_data() → ChemperiumLoader integration
- train() → DeltaLearner orchestration
- predict() → Model loading + prediction
- save()/load() → Registry integration
```

---

## WEEK 3: Real PM7 + Polish (Days 11-15)

### Day 11-12: GridUNESP Submission
```
Create: scripts/gridunesp/
- submit_pm7.sh (SLURM array job)
- monitor.py (check job status)
- collect_results.py (merge output files)

Strategy:
- Split 52k into 500 batches of ~100 molecules
- Array job: #SBATCH --array=1-500
- Each job: ~100 molecules × 57s = 1.5h per job
- Total: ~3-4 days with queue wait
```

### Day 13-14: Integration + Polish
```
Tasks:
- Merge PM7 results with CBS data
- Train final models on real delta
- Tune hyperparameters (Optuna optional)
- Document results (RMSE before/after)

Update CLI:
- Add progress tracking for long operations
- Add --resume flag for interrupted jobs
- Improve error messages
```

### Day 15: Final Validation
```
Tests:
- Full integration test with real data
- Performance benchmarks
- Edge case handling

Documentation:
- Update README with usage examples
- Add CLI --help documentation
```

---

## File Modification Summary

### Week 1 (Core ML)
| File | Action |
|------|--------|
| `src/grimperium/models/kernel_ridge.py` | Implement fit/predict |
| `src/grimperium/models/xgboost_model.py` | Implement fit/predict |
| `src/grimperium/models/delta_ensemble.py` | Implement ensemble |
| `src/grimperium/core/metrics.py` | Implement all metrics |
| `src/grimperium/core/delta_learning.py` | Implement DeltaLearner |
| `tests/fixtures/conftest.py` | Add mock PM7 generator |
| `tests/integration/test_delta_learning.py` | Create E2E tests |

### Week 2 (PM7 + CLI)
| File | Action |
|------|--------|
| `src/grimperium/data/semiempirical.py` | Implement PM7 calc |
| `src/grimperium/cli/__init__.py` | Create |
| `src/grimperium/cli/main.py` | Create Rich CLI |
| `src/grimperium/models/registry.py` | Create |
| `src/grimperium/api.py` | Complete implementation |
| `pyproject.toml` | Add rich, update scripts |

### Week 3 (Production)
| File | Action |
|------|--------|
| `scripts/gridunesp/submit_pm7.sh` | Create SLURM script |
| `scripts/gridunesp/collect_results.py` | Create |

---

## Success Criteria

1. **Week 1 Gate**: Delta-learning RMSE < Direct ML RMSE (with mock PM7)
2. **Week 2 Gate**: CLI functional, PM7 calculation works for 100 molecules
3. **Week 3 Gate**: Full pipeline with real PM7 data, documented results

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| PM7 calc too slow locally | Use GridUNESP from Week 2 |
| MOPAC not installed | Provide Docker container |
| Delta-learning doesn't improve | Investigate error patterns, adjust features |
| GridUNESP queue delays | Start subset locally, scale on cluster |

---

## NOT Doing (Deferred to v0.3)

- Streamlit dashboard
- CREST conformer search (start with single conformer)
- Advanced hyperparameter tuning (Optuna)
- Docker deployment
- Multi-property prediction (entropy, heat capacity)

---

# REFINAMENTOS CRÍTICOS (5 Bloqueadores Resolvidos)

## Bloqueador 1: DeltaLearner.fit() EXPLÍCITO

```python
# src/grimperium/core/delta_learning.py

class DeltaLearner:
    """
    Delta-Learning Model: aprende APENAS a correção (delta) entre CBS e PM7.

    Arquitetura:
    - Input: X (features), y_cbs (H298 CBS real), y_pm7 (H298 PM7 approx)
    - Internamente: y_delta = y_cbs - y_pm7  ← SÓ ISSO que treina
    - Treina: KRR(X, y_delta) + XGB(X, y_delta)  ← AMBOS em delta
    - Ensemble: weighted_avg = 0.5*KRR + 0.5*XGB
    - Predição: y_pred = y_pm7 + delta_pred  ← Composição
    """

    def __init__(self, w_krr: float = 0.5, w_xgb: float = 0.5):
        self.w_krr = w_krr
        self.w_xgb = w_xgb
        self.krr = KernelRidgeRegressor(alpha=1.0, gamma=0.1)
        self.xgb = XGBoostRegressor(max_depth=5, learning_rate=0.1)
        self.scaler = StandardScaler()
        self.is_fitted = False

    def fit(self, X: np.ndarray, y_cbs: np.ndarray, y_pm7: np.ndarray) -> "DeltaLearner":
        """
        Fit DeltaLearner.

        Process (EXPLÍCITO para evitar bugs):
        1. Calcula delta: y_delta = y_cbs - y_pm7
        2. Normaliza X: StandardScaler
        3. Treina KRR: krr.fit(X_scaled, y_delta)
        4. Treina XGB: xgb.fit(X_scaled, y_delta)
        5. Ensemble pronto para weighted avg
        """
        # STEP 1: DELTA (EXPLÍCITO!)
        y_delta = y_cbs - y_pm7

        # STEP 2: Normalize features
        X_scaled = self.scaler.fit_transform(X)

        # STEP 3-4: Train both models on DELTA (not y_cbs!)
        self.krr.fit(X_scaled, y_delta)
        self.xgb.fit(X_scaled, y_delta)

        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray, y_pm7: np.ndarray) -> np.ndarray:
        """
        Prediz H298_CBS via delta-learning.

        Process (EXPLÍCITO):
        1. Normaliza X com scaler fitted
        2. Prediz delta: delta_pred = w_krr*KRR + w_xgb*XGB
        3. Composição: y_pred = y_pm7 + delta_pred
        """
        X_scaled = self.scaler.transform(X)

        krr_pred = self.krr.predict(X_scaled)
        xgb_pred = self.xgb.predict(X_scaled)

        # Weighted average ensemble
        delta_pred = self.w_krr * krr_pred + self.w_xgb * xgb_pred

        # COMPOSIÇÃO (EXPLÍCITA!)
        y_pred = y_pm7 + delta_pred
        return y_pred

    def evaluate(self, X: np.ndarray, y_cbs: np.ndarray, y_pm7: np.ndarray) -> dict:
        """Evaluate with all metrics."""
        y_pred = self.predict(X, y_pm7)
        return compute_all_metrics(y_cbs, y_pred)
```

---

## Bloqueador 2: DeltaEnsemble Week 1 = Weighted Avg Simples

```python
# src/grimperium/models/delta_ensemble.py

class DeltaLearningEnsemble:
    """
    Simple Ensemble para Delta-Learning.

    Week 1: Weighted average (w_krr=0.5, w_xgb=0.5)
    - Porquê simples? Week 1 é VALIDAÇÃO de hipótese, não otimização.
    - Se gate passar, mantém simples.

    Week 2+ (SE NECESSÁRIO):
    - Stacking com meta-learner
    - Optuna para tunar w_krr/w_xgb
    - Mas SÓ SE weighted avg não funcionar bem.
    """

    def __init__(
        self,
        w_krr: float = 0.5,
        w_xgb: float = 0.5,
        krr_params: Optional[dict] = None,
        xgb_params: Optional[dict] = None,
    ):
        self.w_krr = w_krr
        self.w_xgb = w_xgb

        # Week 1 defaults
        krr_defaults = {"alpha": 1.0, "gamma": 0.1}
        xgb_defaults = {"max_depth": 5, "learning_rate": 0.1, "n_estimators": 100}

        self.krr = KernelRidgeRegressor(**(krr_params or krr_defaults))
        self.xgb = XGBoostRegressor(**(xgb_params or xgb_defaults))
        self.is_fitted = False

    def fit(self, X: np.ndarray, y: np.ndarray) -> "DeltaLearningEnsemble":
        """Fit both models on same target (delta)."""
        self.krr.fit(X, y)
        self.xgb.fit(X, y)
        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Simple weighted average: w_krr*KRR + w_xgb*XGB"""
        krr_pred = self.krr.predict(X)
        xgb_pred = self.xgb.predict(X)
        return self.w_krr * krr_pred + self.w_xgb * xgb_pred

    def predict_individual(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Return individual predictions for analysis."""
        return self.krr.predict(X), self.xgb.predict(X)
```

---

## Bloqueador 3: Fixtures Separadas (synthetic vs real)

```python
# tests/fixtures/conftest.py

import pytest
import numpy as np
from typing import Tuple

# ============================================================
# MOCK PM7 FUNCTIONS
# ============================================================

def create_realistic_mock_pm7(
    y_cbs: np.ndarray,
    X_basic: np.ndarray,
    seed: int = 42
) -> np.ndarray:
    """
    Gerar H298_PM7 mock com padrões realistas.

    Components:
    1. Base bias: PM7 overestimates ~5 kcal/mol
    2. Size-dependent: error scales with sqrt(nheavy)
    3. Magnitude bias: 2% of |y_cbs|
    4. Gaussian noise: ~7 kcal/mol std
    """
    np.random.seed(seed)
    nheavy = X_basic[:, 0]

    base_bias = -5.0
    size_error = (1 + np.sqrt(np.maximum(nheavy, 1))) * np.random.normal(0, 1.5, len(y_cbs))
    magnitude_bias = (np.abs(y_cbs) / 100) * np.random.normal(0, 3, len(y_cbs))
    gaussian_noise = np.random.normal(0, 7, len(y_cbs))

    total_error = base_bias + size_error + magnitude_bias + gaussian_noise
    return y_cbs - total_error


def create_enriched_features(X_basic: np.ndarray) -> np.ndarray:
    """
    Expand [nheavy, charge, mult] to 10 dimensions.

    Output: [nheavy, charge, mult, nheavy², sqrt(nheavy),
             nheavy×charge, nheavy×mult, charge², bias]
    """
    nheavy = X_basic[:, 0:1]
    charge = X_basic[:, 1:2]
    mult = X_basic[:, 2:3]

    return np.hstack([
        nheavy, charge, mult,                    # Original (3)
        nheavy ** 2,                             # Polynomial (1)
        np.sqrt(np.abs(nheavy) + 1),             # Polynomial (1)
        nheavy * charge,                         # Interaction (1)
        nheavy * mult,                           # Interaction (1)
        charge ** 2,                             # Polynomial (1)
        mult ** 2,                               # Polynomial (1)
        np.ones_like(nheavy)                     # Bias term (1)
    ])  # Total: 10 dimensions


# ============================================================
# FIXTURES
# ============================================================

@pytest.fixture
def synthetic_data_1k() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Dados SINTÉTICOS com features enriquecidas.

    Uso: Testes rápidos (<1s), validação controlada
    Returns: (X_enriched, y_cbs, y_pm7_mock, y_delta)
    """
    np.random.seed(42)
    n = 1000

    nheavy = np.random.randint(1, 25, size=n)
    charge = np.random.choice([-2, -1, 0, 1, 2], size=n)
    mult = np.random.choice([1, 2, 3], size=n)

    X_basic = np.column_stack([nheavy, charge, mult]).astype(float)
    y_cbs = np.random.normal(loc=-40, scale=50, size=n)

    X = create_enriched_features(X_basic)
    y_pm7 = create_realistic_mock_pm7(y_cbs, X_basic, seed=42)
    y_delta = y_cbs - y_pm7

    return X, y_cbs, y_pm7, y_delta


@pytest.fixture
def real_data_1k() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Dados REAIS CBS com features enriquecidas + mock PM7.

    Uso: Testes de integração, validação com dados verdadeiros
    Timeline: ~2-3s para load + process
    Returns: (X_enriched, y_cbs, y_pm7_mock, y_delta)
    """
    from grimperium.data.loader import ChemperiumLoader

    loader = ChemperiumLoader()
    df = loader.load_thermo_cbs_opt(max_rows=1000)

    X_basic = df[['nheavy', 'charge', 'multiplicity']].values.astype(float)
    y_cbs = df['H298_cbs'].values.astype(float)

    X = create_enriched_features(X_basic)
    y_pm7 = create_realistic_mock_pm7(y_cbs, X_basic, seed=42)
    y_delta = y_cbs - y_pm7

    return X, y_cbs, y_pm7, y_delta
```

---

## Bloqueador 4: Normalização Encapsulada

```python
# src/grimperium/models/kernel_ridge.py

from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import StandardScaler

class KernelRidgeRegressor(BaseModel):
    """KRR with encapsulated scaling."""

    def __init__(self, alpha: float = 1.0, gamma: float = 0.1):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.scaler = StandardScaler()  # ← Encapsulado
        self._model = None

    def fit(self, X: np.ndarray, y: np.ndarray, **kwargs) -> "KernelRidgeRegressor":
        X_scaled = self.scaler.fit_transform(X)  # ← Automático
        self._model = KernelRidge(kernel='rbf', alpha=self.alpha, gamma=self.gamma)
        self._model.fit(X_scaled, y)
        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        if not self.is_fitted:
            raise ValueError("Model not fitted. Call fit() first.")
        X_scaled = self.scaler.transform(X)  # ← Usa scaler fitted
        return self._model.predict(X_scaled)


# src/grimperium/models/xgboost_model.py

import xgboost as xgb
from sklearn.preprocessing import StandardScaler

class XGBoostRegressor(BaseModel):
    """XGBoost with encapsulated scaling + early stopping."""

    def __init__(
        self,
        n_estimators: int = 100,
        max_depth: int = 5,
        learning_rate: float = 0.1,
        early_stopping_rounds: int = 10,
    ):
        super().__init__()
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.early_stopping_rounds = early_stopping_rounds
        self.scaler = StandardScaler()  # ← Encapsulado
        self._model = None

    def fit(
        self,
        X: np.ndarray,
        y: np.ndarray,
        eval_set: Optional[Tuple] = None,
        **kwargs
    ) -> "XGBoostRegressor":
        X_scaled = self.scaler.fit_transform(X)

        self._model = xgb.XGBRegressor(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            learning_rate=self.learning_rate,
            random_state=42,
        )

        fit_params = {}
        if eval_set is not None:
            X_eval_scaled = self.scaler.transform(eval_set[0])
            fit_params['eval_set'] = [(X_eval_scaled, eval_set[1])]
            fit_params['verbose'] = False

        self._model.fit(X_scaled, y, **fit_params)
        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        if not self.is_fitted:
            raise ValueError("Model not fitted. Call fit() first.")
        X_scaled = self.scaler.transform(X)
        return self._model.predict(X_scaled)
```

---

## Bloqueador 5: Hyperparameters Documentados

```python
# tests/experiments/test_validate_hypothesis.py

"""
HYPERPARAMETER STRATEGY - Week 1 vs Week 2

══════════════════════════════════════════════════════════════
WEEK 1 (Fast Validation - Defaults)
══════════════════════════════════════════════════════════════

KernelRidgeRegressor:
  alpha = 1.0       # L2 regularization (sklearn default)
  gamma = 0.1       # RBF kernel scale (moderate)

  Rationale: Standard defaults that work for most problems.
  If fails: Try alpha=[0.01, 0.1, 10], gamma=[0.01, 1.0]

XGBoostRegressor:
  n_estimators = 100      # Enough trees for learning
  max_depth = 5           # Shallow (prevent overfitting)
  learning_rate = 0.1     # Standard rate
  early_stopping = 10     # Stop if no improvement

  Rationale: Conservative settings for small dataset (1k).
  If fails: Increase n_estimators=200, max_depth=7

DeltaLearner:
  w_krr = 0.5, w_xgb = 0.5  # Equal weighting

  Rationale: Week 1 is validation, not optimization.
  If fails: Try w_krr=0.3, w_xgb=0.7 (XGB usually better)

══════════════════════════════════════════════════════════════
GATE PASS/FAIL CRITERIA
══════════════════════════════════════════════════════════════

✅ PASS IF ALL:
  - RMSE_delta < RMSE_direct (core hypothesis)
  - RMSE_delta < 20 kcal/mol (reasonable accuracy)
  - R² > 0.6 (explains >60% variance)
  - Improvement > 10% (meaningful gain)

❌ FAIL IF ANY:
  - RMSE_delta >= RMSE_direct
  - R² < 0.5
  - Improvement < 5%

══════════════════════════════════════════════════════════════
IF GATE FAILS - Debug Checklist
══════════════════════════════════════════════════════════════

1. Check mock PM7 distribution:
   - Is std(y_delta) much smaller than std(y_cbs)?
   - Should be ~10-15 vs ~50 kcal/mol

2. Check features:
   - Add more polynomial terms
   - Try log(nheavy) if nheavy distribution is skewed

3. Tune hyperparameters:
   - KRR: GridSearchCV over alpha/gamma
   - XGB: Increase n_estimators, tune learning_rate

4. Check for data leakage:
   - y_pm7_mock should NOT correlate perfectly with y_cbs
   - Verify: corr(y_pm7, y_cbs) should be ~0.85-0.95, not 1.0
"""
```

---

## Arquitetura Refinada - Resumo

```
tests/fixtures/conftest.py
├── create_realistic_mock_pm7()     # Mock com 4 components de erro
├── create_enriched_features()      # 3 dims → 10 dims
├── synthetic_data_1k              # Fixture sintética
└── real_data_1k                   # Fixture com dados CBS reais

src/grimperium/models/
├── kernel_ridge.py
│   └── KernelRidgeRegressor       # fit/predict com StandardScaler encapsulado
├── xgboost_model.py
│   └── XGBoostRegressor           # fit/predict com early_stopping
└── delta_ensemble.py
    └── DeltaLearningEnsemble      # Week 1 = weighted avg simples

src/grimperium/core/
├── metrics.py
│   ├── rmse(), mae(), r2_score(), mape(), max_error()
│   └── compute_all_metrics()
└── delta_learning.py
    └── DeltaLearner               # fit(X, y_cbs, y_pm7) → treina em y_delta

tests/experiments/
└── test_validate_hypothesis.py    # Phase 0 validation (pytest + standalone)
```
