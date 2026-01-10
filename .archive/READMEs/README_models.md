# Grimperium Models Module

O módulo `models` contém implementações de modelos de Machine Learning otimizados para química computacional e Delta Learning.

## Componentes

### Base Model

Classe base abstrata para todos os modelos.

```python
from grimperium.models.base import BaseModel

# Todos os modelos herdam de BaseModel
# Interface comum para fit(), predict(), score()
```

**Interface padronizada:**
- `fit(X, y)` - Treinar modelo
- `predict(X)` - Fazer predições
- `score(X, y)` - Avaliar desempenho
- `save(path)` / `load(path)` - Persistência

### Kernel Ridge Regression

Modelo de regressão com kernel para propriedades moleculares.

```python
from grimperium.models.kernel_ridge import KernelRidgeModel

# Criar e treinar modelo
model = KernelRidgeModel(
    kernel='rbf',
    alpha=0.1,
    gamma='scale'
)
model.fit(X_train, y_train)

# Fazer predições
predictions = model.predict(X_test)
```

**Kernels disponíveis:**
- RBF (Radial Basis Function)
- Linear
- Polynomial
- Custom kernels

**Aplicações:**
- Predição de energias
- Propriedades eletrônicas
- Pequenos datasets (~1k amostras)

### XGBoost Model

Modelo de boosting otimizado para datasets grandes.

```python
from grimperium.models.xgboost_model import XGBoostModel

# Criar e treinar modelo
model = XGBoostModel(
    n_estimators=1000,
    max_depth=6,
    learning_rate=0.1
)
model.fit(X_train, y_train)

# Predições com importância de features
predictions = model.predict(X_test)
feature_importance = model.get_feature_importance()
```

**Características:**
- Alta performance em datasets grandes
- Feature importance nativa
- Otimização GPU (opcional)
- Early stopping

**Aplicações:**
- Datasets grandes (>10k amostras)
- Feature selection
- Propriedades complexas

### Delta Ensemble

Ensemble de modelos para Delta Learning.

```python
from grimperium.models.delta_ensemble import DeltaEnsemble

# Criar ensemble
ensemble = DeltaEnsemble(
    models=[
        KernelRidgeModel(),
        XGBoostModel(),
    ],
    weights=[0.5, 0.5]
)

# Treinar ensemble
ensemble.fit(X_train, y_train)

# Predições combinadas
predictions = ensemble.predict(X_test)
```

**Estratégias de ensemble:**
- Weighted averaging
- Stacking
- Voting
- Custom aggregation

**Vantagens:**
- Reduz overfitting
- Melhora generalização
- Combina pontos fortes de diferentes modelos

## Arquivos

- `base.py` - Classe base para todos os modelos
- `kernel_ridge.py` - Implementação Kernel Ridge Regression
- `xgboost_model.py` - Implementação XGBoost
- `delta_ensemble.py` - Sistema de ensemble

## Comparação de Modelos

| Modelo | Dataset Size | Speed | Accuracy | Feature Importance |
|--------|-------------|-------|----------|-------------------|
| Kernel Ridge | Small-Medium | Medium | High | ❌ |
| XGBoost | Large | Fast | High | ✅ |
| Delta Ensemble | Any | Slow | Very High | ✅ |

## Workflow de Seleção de Modelo


## Uso com Delta Learning

```python
from grimperium.core import DeltaLearner
from grimperium.models import KernelRidgeModel, XGBoostModel

# Criar DeltaLearner com modelos específicos
learner = DeltaLearner(
    base_model=KernelRidgeModel(),  # Modelo base (rápido)
    delta_model=XGBoostModel()      # Modelo delta (preciso)
)

# Treinar
learner.fit(X_train, y_train)

# Predições hierárquicas
predictions = learner.predict(X_test)
```

## Status

✅ **Production Ready**
- Todos os modelos testados
- Hyperparameter tuning disponível
- GPU support (XGBoost)

## Dependências

- scikit-learn
- XGBoost
- NumPy
- SciPy

## Ver também

- [Core Module](README_core.md) - Delta Learning algorithm
- [Data Module](README_data.md) - Data loading
- [Utils Module](README_utils.md) - Feature engineering
- [Documentation](../docs/build/html/grimperium.models.html) - API completa
