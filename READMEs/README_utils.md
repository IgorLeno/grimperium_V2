# Grimperium Utils Module

O módulo `utils` fornece utilitários e ferramentas auxiliares para validação, logging e feature engineering.

## Componentes

### Validation

Sistema de validação de dados e modelos.

```python
from grimperium.utils.validation import validate_data, validate_model

# Validar dados de entrada
is_valid, errors = validate_data(
    X_train,
    expected_features=['smiles', 'descriptors'],
    check_nan=True,
    check_duplicates=True
)

# Validar modelo treinado
is_valid, report = validate_model(
    model,
    X_test,
    y_test,
    min_r2=0.7
)
```

**Validações disponíveis:**
- Formato de dados
- Valores ausentes (NaN)
- Duplicatas
- Ranges de features
- Qualidade do modelo
- Cross-validation

### Logging

Sistema de logging estruturado para experimentos.

```python
from grimperium.utils.logging import setup_logger, log_experiment

# Configurar logger
logger = setup_logger(
    name='grimperium',
    level='INFO',
    log_file='experiments.log'
)

# Registrar experimento
log_experiment(
    experiment_name='delta_learning_v1',
    params={'alpha': 0.1, 'kernel': 'rbf'},
    metrics={'rmse': 0.05, 'r2': 0.95}
)
```

**Características:**
- Logging estruturado (JSON)
- Rastreamento de experimentos
- Métricas temporais
- Integração com MLflow (opcional)

### Feature Engineering

Ferramentas para criação e transformação de features moleculares.

```python
from grimperium.utils.feature_engineering import (
    compute_descriptors,
    create_fingerprints,
    normalize_features,
    select_features
)

# Computar descritores moleculares
descriptors = compute_descriptors(
    molecules,
    descriptor_types=['molecular_weight', 'logP', 'TPSA']
)

# Criar fingerprints
fingerprints = create_fingerprints(
    molecules,
    fp_type='morgan',
    radius=2,
    n_bits=2048
)

# Normalizar features
X_normalized = normalize_features(
    X_train,
    method='standard',
    clip_outliers=True
)

# Seleção de features
X_selected, selected_indices = select_features(
    X_train,
    y_train,
    method='variance_threshold',
    n_features=100
)
```

**Tipos de descritores:**
- Molecular weight
- LogP (lipophilicity)
- TPSA (Topological Polar Surface Area)
- Rotatable bonds
- Aromatic rings
- Custom descriptors

**Tipos de fingerprints:**
- Morgan (circular)
- MACCS keys
- Topological
- ECFP
- Custom fingerprints

**Métodos de normalização:**
- Standard scaling (z-score)
- Min-max scaling
- Robust scaling
- Quantile transformation

**Métodos de seleção:**
- Variance threshold
- Mutual information
- Correlation-based
- L1 regularization
- Recursive feature elimination

## Arquivos

- `validation.py` - Sistema de validação
- `logging.py` - Sistema de logging
- `feature_engineering.py` - Feature engineering tools

## Workflow Completo

```python
from grimperium.utils import (
    validate_data,
    setup_logger,
    compute_descriptors,
    normalize_features,
    select_features
)

# 1. Setup logging
logger = setup_logger('experiment_1')

# 2. Validar dados
is_valid, errors = validate_data(raw_data)
if not is_valid:
    logger.error(f"Data validation failed: {errors}")
    from grimperium.exceptions import DataValidationError
    raise DataValidationError(f"Data validation failed: {errors}")

# 3. Feature engineering
descriptors = compute_descriptors(molecules)
X_normalized = normalize_features(descriptors)
X_selected, indices = select_features(X_normalized, y_train)

# 4. Treinar modelo
from grimperium.core import DeltaLearner
learner = DeltaLearner()
learner.fit(X_selected, y_train)

# 5. Validar modelo
from grimperium.utils.validation import validate_model
is_valid, report = validate_model(learner, X_test, y_test)
logger.info(f"Model validation: {report}")
```

## Integração com Pipeline

```python
from sklearn.pipeline import Pipeline
from grimperium.utils.feature_engineering import (
    DescriptorComputer,
    FeatureNormalizer,
    FeatureSelector
)

# Criar pipeline completo
pipeline = Pipeline([
    ('descriptors', DescriptorComputer()),
    ('normalizer', FeatureNormalizer(method='standard')),
    ('selector', FeatureSelector(n_features=100)),
    ('model', DeltaLearner())
])

# Treinar pipeline
pipeline.fit(X_train, y_train)

# Fazer predições
predictions = pipeline.predict(X_test)
```

## Logging de Experimentos

```python
from grimperium.utils.logging import ExperimentLogger

# Criar logger de experimento
exp_logger = ExperimentLogger('experiments/')

# Registrar experimento
with exp_logger.experiment('delta_v1') as exp:
    # Registrar parâmetros
    exp.log_params({
        'base_model': 'kernel_ridge',
        'delta_model': 'xgboost',
        'alpha': 0.1
    })

    # Treinar modelo
    model.fit(X_train, y_train)

    # Registrar métricas
    exp.log_metrics({
        'train_rmse': 0.03,
        'test_rmse': 0.05,
        'r2': 0.95
    })

    # Salvar artefatos
    exp.log_artifact(model, 'model.pkl')
```

## Status

✅ **Production Ready**
- Validação robusta
- Logging estruturado
- Feature engineering completo

## Dependências

- NumPy
- scikit-learn
- RDKit (para descritores moleculares)
- Pandas

## Ver também

- [Core Module](README_core.md) - Delta Learning
- [Data Module](README_data.md) - Data loading
- [Models Module](README_models.md) - ML models
- [Documentation](https://grimperium.readthedocs.io/en/latest/grimperium.utils.html) - API completa (para gerar localmente: `cd docs && make html`, depois abra `docs/build/html/index.html`)
