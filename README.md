# Grimperium - Delta Learning Framework

**Status:** Production Ready (Phase A - Active)
**Python:** 3.10+
**License:** MIT
**Last Updated:** 2026-01-18

---

## Status do Projeto

| Fase | Status | Descrição |
|------|--------|-----------|
| **Phase A** | ✅ Completo | CREST PM7 baseline & validation |
| **Phase B** | ⏳ Pronto | Delta-Learning training (after BATCH 12) |
| **Phase C** | 🔧 BATCH 12 | CLI Critical Fixes (11 bugs) |

### Recent Updates (2026-01-18)
- ✅ Dataset refactoring: CHON (29,568) + PM7 naming clarity
- ✅ Tests passing: 242/262 (20 skipped)
- ✅ Quality gates: mypy, ruff, black ✅
- ✅ Pre-commit hooks: Active
- 🔄 BATCH 12: Ready to start (11 bugs identified)

---

## Objetivo

Framework de Delta Learning para predição de propriedades moleculares usando:
- **Delta Learning** - Correção de métodos semiempíricos (PM7)
- **Ensemble Models** - XGBoost, Kernel Ridge, e combinações
- **Feature Engineering** - Descritores moleculares otimizados
- **Threshold Monitoring** - Detecção de anomalias e falhas

---

## Estrutura do Projeto

```text
grimperium/
├── src/grimperium/          <- Código principal
│   ├── core/                <- Delta Learning, Metrics
│   ├── data/                <- Data loading, Fusion
│   ├── models/              <- ML Models (XGBoost, KRR, Ensemble)
│   ├── utils/               <- Helpers (logging, validation)
│   ├── config.py            <- Configurações globais
│   └── api.py               <- Interface pública
├── tests/                   <- Testes unitários + integração
├── data/                    <- Datasets e baselines
│   ├── thermo_cbs_chon.csv  <- 29,568 molecules (CHON only) - PRIMARY
│   └── thermo_pm7.csv       <- PM7 optimization results - SECONDARY
├── docs/                    <- Documentação
│   ├── CLAUDE.md            <- Guia para Claude Code
│   ├── PHASE-A-START-HERE.md <- Comece aqui
│   └── architecture.md      <- Visão da arquitetura
├── scripts/                 <- Scripts utilitários
└── .archive/                <- Arquivos deprecated
```

---

## Quick Start
### 1. Setup (5 min)

```bash
git clone https://github.com/IgorLeno/grimperium.git
cd grimperium
python -m venv venv && source venv/bin/activate
pip install -e .
```

### 2. Teste Rápido

```bash
pytest tests/ -v
```

### 3. Uso Básico

```python
from grimperium import DeltaLearning
from grimperium.models import XGBoostDelta

# Carregar dados
model = XGBoostDelta()
model.train(X_train, y_train)
predictions = model.predict(X_test)
```

---

## Datasets

Grimperium uses two primary datasets:

### Primary: thermo_cbs_chon.csv
- **29,568 molecules** (CHON only: C, H, O, N)
- **High-accuracy CBS enthalpies** + B3LYP baseline
- **Optimized for delta-learning:** Homogeneous chemistry, learnable corrections
- **Use case:** Model training, validation, benchmarking

### Secondary: thermo_pm7.csv
- **PM7-optimized results** from CREST conformer pipeline
- **Semiempirical baseline** for experimental validation
- **Use case:** Alternative baseline, production predictions

For detailed information, see [`docs/DATASETS.md`](docs/DATASETS.md).

**Quick load:**
```python
from grimperium.data.loader import ChemperiumLoader

loader = ChemperiumLoader()
df = loader.load_thermo_cbs_chon(max_nheavy=50)
print(f"Loaded {len(df)} CHON molecules")
```

---

## Documentação

| Documento | Propósito |
|-----------|-----------|
| `docs/CLAUDE.md` | Guia comportamento Claude Code |
| `docs/architecture.md` | Visão da arquitetura |
| `docs/DATASETS.md` | **NEW** - Dataset reference guide |
| `docs/PHASE-A-START-HERE.md` | Comece aqui para Phase A |
| `docs/PHASE-A-RESULTS.md` | modelo para resultados |
| `docs/feature_engineering.md` | Guia de features |
| `docs/delta_learning_guide.md` | Guia de Delta Learning |

---

## Tecnologias

- **Python:** 3.10+
- **ML:** scikit-learn, XGBoost
- **Data:** pandas, numpy
- **Testing:** pytest, pytest-cov
- **Code Quality:** ruff, black, mypy
- **CI/CD:** GitHub Actions

---

## Instalação Detalhada

```bash
# Clone
git clone https://github.com/IgorLeno/grimperium.git
cd grimperium

# Ambiente virtual
python -m venv venv
source venv/bin/activate  # ou: venv\Scripts\activate (Windows)

# Instalar modo desenvolvimento
pip install -e .

# Instalar dev dependencies
pip install pytest pytest-cov ruff black mypy

# Validar instalação
python -c "from grimperium import *; print('OK')"
```

---

## Testes

```bash
# Todos os testes
pytest tests/ -v

# Com coverage
pytest tests/ --cov=src/grimperium --cov-report=html

# Abrir relatório
open htmlcov/index.html
```

---

## Code Quality

```bash
# Linting
ruff check src/

# Formatting
black src/ tests/

# Type checking
mypy src/ --strict

# Pre-commit hooks
pre-commit run --all-files
```

---

## Para Claude Code

Leia: `docs/CLAUDE.md`

**Resumo:**
- Type hints 100%
- Docstrings completas
- Testes (85%+ coverage)
- Code quality (ruff, black, mypy)
- Pre-commit hooks

---

## Roadmap

- **Phase A** (NOW): Validação com moléculas reais
- **Phase B** (Next): Expansão para 50+ moléculas
- **Phase C** (Later): Otimização e produção

---

## License

MIT License - Veja LICENSE para detalhes

---

## Autor

Igor Leno - São Paulo, BR
