# Grimperium Documentation Index

Documentação modular do framework Grimperium para Delta Learning em química computacional.

## 📚 Módulos

### [Core Module](README_core.md)
Algoritmo Delta Learning e métricas de avaliação.

**Componentes principais:**
- `DeltaLearner` - Implementação do algoritmo Delta Learning
- `Metrics` - Funções de avaliação (RMSE, MAE, R²)

**Status:** ✅ Production Ready

---

### [Data Module](README_data.md)
Carregamento, processamento e fusão de dados.

**Componentes principais:**
- `DataLoader` - Carregamento de múltiplos formatos
- `SemiempiricalLoader` - Interface para dados PM6/PM7
- `DataFusion` - Sistema de fusão de dados

**Status:** ✅ Production Ready

---

### [Models Module](README_models.md)
Modelos de Machine Learning para química computacional.

**Componentes principais:**
- `BaseModel` - Classe base abstrata
- `KernelRidgeModel` - Regressão com kernel
- `XGBoostModel` - Gradient boosting
- `DeltaEnsemble` - Sistema de ensemble

**Status:** ✅ Production Ready

---

### [Utils Module](README_utils.md)
Utilitários e ferramentas auxiliares.

**Componentes principais:**
- `Validation` - Validação de dados e modelos
- `Logging` - Sistema de logging estruturado
- `Feature Engineering` - Descritores e fingerprints moleculares

**Status:** ✅ Production Ready

---

## 🚀 Quick Start

```python
# 1. Carregar dados
from grimperium.data import DataLoader
loader = DataLoader()
data = loader.load_from_file("dataset.csv")

# 2. Feature engineering
from grimperium.utils.feature_engineering import compute_descriptors
descriptors = compute_descriptors(data['molecules'])

# 3. Treinar modelo Delta Learning
from grimperium.core import DeltaLearner
from grimperium.models import KernelRidgeModel, XGBoostModel

learner = DeltaLearner(
    base_model=KernelRidgeModel(),
    delta_model=XGBoostModel()
)
learner.fit(descriptors, data['targets'])

# 4. Fazer predições
predictions = learner.predict(new_molecules)
```

## 📊 Arquitetura

```
grimperium/
├── core/          # Delta Learning algorithm
├── data/          # Data loading & fusion
├── models/        # ML models
└── utils/         # Utilities
```

## 🔗 Links Úteis

- [Sphinx Documentation](../docs/build/html/index.html) - API completa
- [CHANGELOG](../CHANGELOG.md) - Histórico de mudanças
- [Technical Report](../TECHNICAL_REPORT.md) - Relatório técnico

## 📈 Status do Projeto

| Métrica | Valor |
|---------|-------|
| Módulos | 4 principais |
| Cobertura de Testes | ~95% |
| Python Versions | 3.9-3.12 |
| Status | Production Ready |

## 💡 Próximos Passos

- Explorar [Feature Engineering](README_utils.md#feature-engineering)
- Configurar [Logging](README_utils.md#logging)
- Experimentar [Ensemble Models](README_models.md#delta-ensemble)
- Ler [Delta Learning Guide](../docs/delta_learning_guide.md)

---

**Última atualização:** 2026-01-07
