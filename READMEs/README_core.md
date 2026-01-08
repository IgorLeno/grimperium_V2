# Grimperium Core Module

O módulo `core` contém a lógica fundamental do framework Grimperium, incluindo o algoritmo de Delta Learning e métricas de avaliação.

## Componentes

### DeltaLearner

Implementação do algoritmo Delta Learning para aprendizado hierárquico.

```python
from grimperium.core.delta_learning import DeltaLearner

# Criar instância do DeltaLearner
learner = DeltaLearner(
    base_model=my_base_model,
    delta_model=my_delta_model
)

# Treinar modelo
learner.fit(X_train, y_train)

# Fazer predições
predictions = learner.predict(X_test)
```

**Características principais:**
- Aprendizado hierárquico em dois níveis (base + delta)
- Suporte para múltiplos tipos de modelos base e delta
- Validação cruzada integrada
- Métricas de desempenho automáticas

### Metrics

Funções para avaliação de desempenho dos modelos.

```python
from grimperium.core.metrics import compute_rmse, compute_mae, compute_r2

# Calcular métricas
rmse = compute_rmse(y_true, y_pred)
mae = compute_mae(y_true, y_pred)
r2 = compute_r2(y_true, y_pred)
```

**Métricas disponíveis:**
- RMSE (Root Mean Squared Error)
- MAE (Mean Absolute Error)
- R² (Coefficient of Determination)
- Custom metrics support

## Arquivos

- `delta_learning.py` - Implementação do DeltaLearner
- `metrics.py` - Funções de avaliação e métricas

## Status

✅ **Production Ready**
- Testes completos
- Documentação atualizada
- API estável

## Dependências

- NumPy
- scikit-learn
- XGBoost (opcional, para modelos delta)

## Ver também

- [Data Module](README_data.md) - Carregamento e processamento de dados
- [Models Module](README_models.md) - Modelos de ML implementados
- [Documentation](../docs/build/html/grimperium.core.html) - API completa
