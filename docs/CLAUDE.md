# CLAUDE.md - Guia de Comportamento

## Seu Papel no Projeto Grimperium

Você é um **implementador de código production-ready** com responsabilidades claras.

### Responsabilidades

#### 1. CÓDIGO
- Implementar features com qualidade production
- Python 3.10+ com type hints em 100% funções
- Docstrings Google-style em classes/funções
- Sem comentários óbvios, apenas lógica complexa

#### 2. TESTES
- Unit tests para toda lógica
- Integration tests para pipelines
- Coverage mínimo: 85%
- Fixtures reutilizáveis em conftest.py
- Pytest markers: @pytest.mark.unit, @pytest.mark.integration

#### 3. PADRÕES
- Black formatting (line-width: 88)
- Ruff linting (sem warnings críticos)
- Mypy strict type checking
- Pre-commit hooks validation

#### 4. DOCUMENTAÇÃO
- Docstrings (Google style)
- Type hints completas
- README em cada package principal
- Inline comments apenas para lógica não-óbvia

---

## CHECKLIST PRÉ-IMPLEMENTAÇÃO

Antes de começar ANY feature:

- [ ] Li o requisito completo
- [ ] Identifiquei arquivos afetados
- [ ] Rodei testes existentes localmente
- [ ] Criei branch `feature/ISSUE-XX` ou `fix/ISSUE-XX`
- [ ] Identificado escopo exato (o que NÃO fazer)

---

## CHECKLIST PÓS-IMPLEMENTAÇÃO

Antes de submeter:

- [ ] Código escrito com type hints 100%
- [ ] Docstrings em classes/funções principais
- [ ] Testes unitários criados (85%+ coverage)
- [ ] Testes integration criados (se aplicável)
- [ ] Rodou: `pytest tests/`
- [ ] Rodou: `ruff check src/`
- [ ] Rodou: `black src/ tests/`
- [ ] Rodou: `mypy src/`
- [ ] Sem erros críticos (warnings ok)
- [ ] Commit message: "feat/fix: ISSUE-XX - descrição concisa"
- [ ] Pronto para produção (nada de TODOs)

---

## NUNCA FAÇA ISTO

#### Testes Incompletos

- Não deixe testes com `@pytest.mark.skip`
- Testes devem passar 100% ou não existir

#### Lógica Não Documentada

- Sem docstrings em funções públicas
- Sem type hints
- Sem exemplos de uso em docstrings

#### Ignorar Linting

- Não commite com erros ruff/black/mypy
- Pre-commit hooks devem passar

#### Mudar Estrutura

- Não crie diretórios fora do plano
- Não renomeie pastas sem aprovação
- Não delete código sem justificar

#### Código Incompleto

- Sem TODOs ou XXX comments
- Sem placeholders
- Tudo deve funcionar 100%

---

## SEMPRE FAÇA ISTO

#### Type Hints Completas

```python
def process_molecules(
    data: pd.DataFrame,
    threshold: float = 0.95
) -> tuple[np.ndarray, dict[str, float]]:
    """Process molecules with quality threshold."""
    pass
```

#### Docstrings Google Style

```python
def threshold_alert(
    value: float,
    threshold: float
) -> bool:
    """Check if value exceeds threshold and raise alert.
    
    Args:
        value: Current measurement value
        threshold: Alert threshold value
        
    Returns:
        True if alert triggered, False otherwise
        
    Raises:
        ValueError: If threshold is negative
        
    Example:
        >>> threshold_alert(0.95, 0.90)
        True
    """
    pass
```

#### Testes Estruturados

```python
# tests/unit/test_threshold_monitor.py
import pytest
from grimperium.core import ThresholdMonitor

@pytest.mark.unit
class TestThresholdMonitor:
    @pytest.fixture
    def monitor(self):
        return ThresholdMonitor(threshold=0.95)
    
    def test_alert_triggered(self, monitor):
        assert monitor.check_threshold(0.98) is True
```

**Pre-commit Validation**
```bash
# Antes de commitar
pre-commit run --all-files
# Se falhar, fix automático:
black . && ruff check --fix .
```

---

## ESTRUTURA ATUAL

```bash
src/grimperium/
├── core/                   ← Delta Learning, Metrics
│   ├── delta_learning.py
│   └── metrics.py
├── data/                   ← Data loading, Fusion
│   ├── loader.py
│   ├── semiempirical.py
│   └── fusion.py
├── crest_pm7/              ← Phase A Pipeline (Ativo)
│   ├── pipeline.py
│   ├── conformer_generator.py
│   ├── mopac_optimizer.py
│   └── threshold_monitor.py
├── models/                 ← ML Models (OBSOLETO/Phase B)
│   └── (obsoleto: base.py, xgboost_model.py, etc.)
├── utils/                  ← Helpers
│   └── feature_engineering.py (logging/validation movidos para crest_pm7)
├── config.py               ← Configurações globais
└── api.py                  ← Interface pública

tests/
├── unit/                   ← Testes isolados
├── integration/            ← Testes pipeline
├── fixtures/               ← Mock & Real data (mock_data.py, real_data.py)
├── experiments/            ← Testes de hipóteses
└── conftest.py             ← Fixtures compartilhadas
```

---

## PHASE A: O que FAZER e O que NÃO FAZER

### FAZER em Phase A
- Validar sistema com moléculas reais (10-20)
- Implementar threshold monitoring
- Testes unitários completos
- Integração com logging
- Documentação clara

### NÃO FAZER em Phase A
- Otimização prematura
- Features extras não solicitadas
- Refatorações fora do escopo

---

## VALIDAÇÕES FINAIS

### Run Antes de Submeter

```bash
# 1. Testes
pytest tests/ -v --cov=src/grimperium --cov-report=term-missing

# 2. Linting
ruff check src/ tests/
black --check src/ tests/
mypy src/ --strict

# 3. Pre-commit
pre-commit run --all-files

# 4. Imports
python -c "from grimperium import *; print('OK')"
```

---

## PERGUNTAS FREQUENTES

**P: Posso usar biblioteca X?**  
R: Só se estiver em `pyproject.toml`. Não adicione deps sem pedir.

**P: Como estruturo testes?**  
R: Veja `tests/conftest.py` para fixtures. Unit tests isolados, integration tests com dados reais.

**P: Qual é a convention de nomes?**  
R: snake_case para funções/variáveis, PascalCase para classes, UPPER_CASE para constantes.

**P: Posso refatorar código antigo?**  
R: Só se for relevante ao escopo atual. Deprecated fica em `.archive/`.

---

## Resumo

**Você é production-ready quando:**
- Código escrito (type hints 100%)
- Testes verdes (85%+ coverage)
- Linting passa (ruff, black, mypy)
- Documentado (docstrings, README)
- Pronto para produção (sem TODOs)
