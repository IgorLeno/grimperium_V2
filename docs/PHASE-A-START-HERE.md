# PHASE A - COMECE AQUI

## Bem-vindo ao Grimperium Phase A

**Objetivo:** Validar sistema com moléculas REAIS (10-20 moléculas)  
**Tempo:** ~60 minutos  
**Risco:** Baixo (testes isolados)

---

## TL;DR (Quick Start - 5 min)

```bash
# 1. Setup
cd ~/Projetos/grimperium
python -m venv venv
source venv/bin/activate
pip install -e .

# 2. Teste Rápido (5 min)
python scripts/phase_a_quick_test.py

# 3. Resultado esperado
# Tests Pass: 100%
# Success Rate: >=95%
# Alertas Normais: 0 anormais
```

---

## FULL GUIDE (60 min, 6 passos)

---

## PASSO 1: Setup Environment (5 min)

```bash
cd ~/Projetos/grimperium
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

## PASSO 4: Monitorar Alertas (10 min)

```bash
# Ver logs de alertas
tail -f data/molecules_pm7/logs/phase_a_real_test/alerts.log

# Checks:
# Threshold violations: [0 abnormal]
# Timeout predictions: [within bounds]
# Log structure: [valid JSON]
```

---

## PASSO 5: Coletar Métricas (10 min)

```bash
python -c "
from grimperium.core.metrics import compute_metrics
metrics = compute_metrics('data/molecules_pm7/logs/phase_a_real_test/')
print(f'Success rate: {metrics[\"success_rate\"]:.1%}')
"
```

---

## PASSO 6: Validar Sucesso (5 min)

**Success Criteria:**
- 100% testes passam
- Success rate >=95%
- Zero alertas anormais
- Métricas baseline estabelecidas

```bash
# Checklist final
./scripts/validate_phase_a.sh

# Expected:
# PHASE A VALIDATION: PASSED
# Ready for Phase B
```

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
- [PHASE-A-RESULTS.md](./PHASE-A-RESULTS.md) - Template para resultados
- [architecture.md](./architecture.md) - Visão da arquitetura
