# BATCH 3: Validação da Hipótese Delta-Learning

## Objetivo

Criar teste experimental automatizado que valida a hipótese central do Grimperium:

> **Delta-learning** (treinar em `y_delta = y_cbs - y_pm7`) é estatisticamente melhor que **direct learning** (treinar direto em `y_cbs`).

## Contexto

### Componentes Disponíveis

- **Fixtures** (`tests/fixtures/conftest.py`):
  - `real_data_1k()` → retorna `(X, y_cbs, y_pm7, y_delta)` - 1000 amostras, 10D features
  - `synthetic_data_1k()` → mesma estrutura (fallback se real falhar)

- **DeltaLearner** (`src/grimperium/core/delta_learning.py`):
  - `fit(X, y_cbs, y_pm7)` → treina ensemble em delta internamente
  - `evaluate(X, y_cbs, y_pm7)` → retorna dict com todas as métricas

- **DeltaLearningEnsemble** (`src/grimperium/models/delta_ensemble.py`):
  - `fit(X, y)` → treina em qualquer target
  - `predict(X)` → retorna predições
  - Defaults: `w_krr=0.5, w_xgb=0.5`

- **Metrics** (`src/grimperium/core/metrics.py`):
  - `compute_all_metrics(y_true, y_pred)` → retorna dict com 'rmse', 'mae', 'r2', etc.

## Arquivos a Criar

1. **`tests/experiments/`** (diretório - não existe)
2. **`tests/experiments/__init__.py`** (arquivo vazio)
3. **`tests/experiments/test_validate_hypothesis.py`** (teste principal)

## Estrutura do Teste

### Imports
```python
import numpy as np
import pytest
from sklearn.model_selection import train_test_split

from grimperium.core.delta_learning import DeltaLearner
from grimperium.models.delta_ensemble import DeltaLearningEnsemble
from grimperium.core.metrics import compute_all_metrics
```

### Função de Teste

**Nome**: `test_decision_gate_delta_vs_direct(real_data_1k, synthetic_data_1k)`

## Passos de Implementação

### 1. Data Loading com Fallback
```python
try:
    X, y_cbs, y_pm7, y_delta = real_data_1k
    data_source = "real_data_1k"
except Exception as e:
    X, y_cbs, y_pm7, y_delta = synthetic_data_1k
    data_source = "synthetic_data_1k"
```
- Preferir dados reais para validação realista
- Fallback garante que teste nunca falha por falta de dados
- Logar fonte de dados e shapes básicos

### 2. Train/Test Split Consistente
```python
X_train, X_test, y_cbs_train, y_cbs_test, y_pm7_train, y_pm7_test = train_test_split(
    X, y_cbs, y_pm7,
    test_size=0.2,
    random_state=42
)
```
- **CRÍTICO**: Mesmo split para ambos os modelos
- 80/20 split
- `random_state=42` para reprodutibilidade

### 3. Treinar Modelo Delta (Referência)
```python
model_delta = DeltaLearner()  # Defaults: w_krr=0.5, w_xgb=0.5
model_delta.fit(X_train, y_cbs_train, y_pm7_train)
metrics_delta = model_delta.evaluate(X_test, y_cbs_test, y_pm7_test)

rmse_delta = metrics_delta['rmse']
mae_delta = metrics_delta['mae']
r2_delta = metrics_delta['r2']
```

### 4. Treinar Modelo Direto (Baseline)
```python
model_direct = DeltaLearningEnsemble(w_krr=0.5, w_xgb=0.5)
model_direct.fit(X_train, y_cbs_train)  # Treina DIRETO em y_cbs

y_pred_direct = model_direct.predict(X_test)
metrics_direct = compute_all_metrics(y_cbs_test, y_pred_direct)

rmse_direct = metrics_direct['rmse']
mae_direct = metrics_direct['mae']
r2_direct = metrics_direct['r2']
```
- **CRÍTICO**: Mesma arquitetura (w_krr=0.5, w_xgb=0.5)
- **DIFERENÇA**: Target (y_cbs vs y_delta)

### 5. Comparação e Decision Gate
```python
# Calcular melhorias
rmse_improvement = rmse_direct - rmse_delta
mae_improvement = mae_direct - mae_delta
r2_improvement = r2_delta - r2_direct

# Decision Gate
criterion_1 = rmse_delta < rmse_direct  # Hipótese: delta < direct
criterion_2 = rmse_delta < 20.0          # Qualidade: RMSE < 20 kcal/mol

gate_pass = criterion_1 and criterion_2
```

### 6. Output e Asserção
```python
# Tabela de comparação
print(f"{'Metric':<15} {'Delta':<15} {'Direct':<15} {'Improvement':<15}")
print(f"{'RMSE':<15} {rmse_delta:<15.4f} {rmse_direct:<15.4f} {rmse_improvement:<15.4f}")
print(f"{'MAE':<15} {mae_delta:<15.4f} {mae_direct:<15.4f} {mae_improvement:<15.4f}")
print(f"{'R²':<15} {r2_delta:<15.4f} {r2_direct:<15.4f} {r2_improvement:<15.4f}")

# Decision Gate
if gate_pass:
    print("✅ DECISION GATE: PASS – delta-learning beats direct learning")
else:
    print("❌ DECISION GATE: FAIL")
    if not criterion_1:
        print(f"   Delta ({rmse_delta:.4f}) >= Direct ({rmse_direct:.4f})")
    if not criterion_2:
        print(f"   RMSE ({rmse_delta:.4f}) >= 20.0 kcal/mol")

# Assert (teste falha se gate não passar)
assert gate_pass, f"Decision gate failed: c1={criterion_1}, c2={criterion_2}"
```

## Prints Informativos

### Durante Execução
- 🔬 Header do experimento
- 📊 Data source e shapes
- 🔷 Status Model Delta (training + métricas)
- 🔶 Status Model Direct (training + métricas)
- ⚖️ Tabela de comparação
- 🚦 Decision Gate (critérios + resultado)
- ✅/❌ Resultado final

### Formato Visual
- Usar box-drawing chars (═, ─) para separação
- Width consistente (60 chars)
- Emojis para visual scanning
- Colunas alinhadas na tabela

## Critérios de Validação

### Comportamento Esperado (Normal)
- `rmse_delta`: **8-12 kcal/mol** (com features enriquecidos)
- `rmse_direct`: **15-30 kcal/mol** (maior target = maior erro)
- `criterion_1`: **True** (delta < direct)
- `criterion_2`: **True** (delta < 20.0)
- **Resultado**: Teste **PASSA**

### Comportamento Anormal (Bug)
- Se `rmse_delta >= rmse_direct`: **criterion_1 = False**
- Se `rmse_delta >= 20.0`: **criterion_2 = False**
- **Resultado**: Teste **FALHA** com mensagem clara

## Comando de Execução

```bash
pytest tests/experiments/test_validate_hypothesis.py -v -s
```

Flags:
- `-v`: verbose (mostra nome do teste)
- `-s`: capture=no (mostra prints)

## Qualidade do Código

### Variáveis Chave
- `model_delta`: Instância DeltaLearner
- `model_direct`: Instância DeltaLearningEnsemble (para direct learning)
- `metrics_delta`, `metrics_direct`: Dicts de métricas
- `criterion_1`, `criterion_2`: Booleanos de gate
- `gate_pass`: Decisão final

### Sem Over-Engineering
- ❌ Não fazer loops de tuning
- ❌ Não fazer múltiplas repetições
- ❌ Não introduzir novas dependências
- ✅ Um teste simples, robusto e claro

### Error Handling
- **Fixture loading**: Try/except com fallback
- **Model fitting**: Deixar exceptions propagarem (fail loudly)
- **Predictions**: Confiamos que metrics.py lida com NaN/Inf

## Success Criteria

1. ✅ Teste executa sem import/runtime errors
2. ✅ Dados carregam (real ou synthetic)
3. ✅ Ambos modelos treinam completamente
4. ✅ Métricas computam corretamente
5. ✅ Decision gate avalia e imprime resultado
6. ✅ Teste passa se hipótese validada, falha caso contrário

## Arquivos Críticos (Referência)

- `tests/fixtures/conftest.py` - Fixtures de dados
- `src/grimperium/core/delta_learning.py` - DeltaLearner
- `src/grimperium/models/delta_ensemble.py` - DeltaLearningEnsemble
- `src/grimperium/core/metrics.py` - compute_all_metrics

## Timeline Estimado

- Criar diretório + `__init__.py`: **1 min**
- Criar `test_validate_hypothesis.py`: **5 min**
- Executar teste inicial: **2-3 min** (load + training)
- **Total**: ~10 minutos
