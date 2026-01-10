# PHASE A - COMECE AQUI

## Bem-vindo ao Grimperium Phase A

**Objetivo:** Validar sistema com moléculas REAIS (10-20 moléculas)  
**Tempo:** ~60 minutos  
**Risco:** Baixo (testes isolados)

---

## TL;DR (Quick Start - 5 min)

From the repository root, run:

```bash
# 1. Setup (Clone if needed, then enter the repo)
# git clone <repo-url> && cd grimperium
cd <repo-root>

# 2. Virtual Environment
python -m venv venv
source venv/bin/activate
pip install -e .

# 3. Teste Rápido (5 min)
python scripts/phase_a_quick_test.py

# 4. Resultado esperado
# Tests Pass: 100%
# Success Rate: >=95%
# Alertas Normais: 0 anormais
```

---

## FULL GUIDE (60 min, 6 passos)

---

## PASSO 1: Setup Environment (5 min)

From the repository root, run:

```bash
cd <repo-root>
python --version  # Expected: 3.10+

# Criar virtual env
python -m venv venv
source venv/bin/activate

# Instalar dependências
pip install -e .
pip install pytest pytest-cov

# Verificar
python -c "from grimperium import *; print('OK')"
```

---

## PASSO 2: Carregar Moléculas Teste (5 min)

```bash
# Arquivo: data/molecules_pm7/testing/baselines/phase_a_baseline.json
# Contém: 10-20 moléculas reais para teste
# Formato: {"smiles": [...], "props": {...}}

python -c "
import json
with open('data/molecules_pm7/testing/baselines/phase_a_baseline.json') as f:
    data = json.load(f)
    print(f'Loaded {len(data[\"smiles\"])} molecules')
"
```

---

## PASSO 3: Rodar Testes do Sistema (15 min)

```bash
pytest tests/ -v --cov=src/grimperium --cov-report=term-missing

# Expected output:
# Test PASSED
# Coverage: >85%
```

---

## PASSO 4: Final Validation (Phase A completion)

**Success Criteria:**
- 100% testes passam
- Success rate >=95%
- Zero alertas anormais
- Métricas baseline estabelecidas

**Note:** Após completar PASSO 1, 2 e 3, você terá validado os componentes core do Phase A. Os passos de monitoramento avançado (alert logs, metrics collection, automated validation) serão implementados durante a execução completa do Phase A.

---

## Phase A Advanced (After Core Steps Complete)

Os seguintes recursos serão implementados durante a execução completa do Phase A:

### Alert Monitoring

- Arquivo: `data/molecules_pm7/logs/phase_a_real_test/alerts.log`
- Monitoramento de threshold violations e timeout predictions
- Estrutura de logs em JSON

### Metrics Collection

- Módulo: `grimperium.core.metrics.compute_all_metrics()`
- Cálculo automático de success rate e outras métricas
- Geração de relatórios

### Automated Validation

- Script: `python scripts/phase_a_quick_test.py`
- Validação automatizada dos success criteria
- Checklist completo de Phase A

Para detalhes do checklist completo de Phase A, veja `docs/PHASE-A-RESULTS.md`.

**Full Phase A execution guide:** see `docs/PHASE-A-RESULTS.md` for success criteria and recording results.

---

## Próximos Passos

Após Phase A:
- Phase B: Expandir para 50+ moléculas, calibrar thresholds
- Phase C: Otimização e produção

---

## Troubleshooting

### Import Error

```bash
pip install -e .
```


### Testes Falhando

```bash
pytest tests/ -v --tb=long
```

### Coverage Baixo

```bash
pytest tests/ --cov=src/grimperium --cov-report=html
# Abrir htmlcov/index.html
```

---

## Recursos

- [CLAUDE.md](./CLAUDE.md) - Guia de comportamento
- [PHASE-A-RESULTS.md](./PHASE-A-RESULTS.md) - modelo para resultados
- [architecture.md](./architecture.md) - Visão da arquitetura
