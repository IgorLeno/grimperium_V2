# 🧬 CLAUDE.md - GRIMPERIUM PROJECT (MERGED)

**Claude Code Agent | Grimperium Delta Learning Framework**  
**Project Focus:** Phase A (Threshold Monitoring + Timeout Prediction)  
**Last Updated:** Jan 10, 2026

***

## 🎯 Seu Papel no Projeto Grimperium

Você é um **implementador de código production-ready** com responsabilidades claras e escopo bem definido.

### Responsabilidades Principais

#### 1. 💻 CÓDIGO
- Implementar features com qualidade production
- Python 3.10+ com **type hints em 100%** das funções
- Docstrings **Google-style** em classes/funções
- Apenas comentários para lógica complexa (sem óbvios)
- Sem TODOs, XXX, FIXME ou placeholders

#### 2. ✅ TESTES
- Unit tests para toda lógica nova
- Integration tests para pipelines Phase A
- Coverage **mínimo 85%** (objetivo: 95%+)
- Fixtures reutilizáveis em `conftest.py`
- Pytest markers: `@pytest.mark.unit`, `@pytest.mark.integration`
- Testes devem passar **100%** ou não existir

#### 3. 🔍 PADRÕES
- **Black** formatting (line-width: 88)
- **Ruff** linting (**ZERO WARNINGS** obrigatório)
- **Mypy** strict type checking (sem errors)
- **Pre-commit** hooks validation (antes de cada commit)

#### 4. 📚 DOCUMENTAÇÃO
- Docstrings Google-style em functions/classes públicas
- Type hints completas (100% obrigatório)
- README em cada package principal
- Exemplos de uso em docstrings (`Example:` section)
- Comentários apenas para lógica não-óbvia

***

## 📋 CHECKLIST PRÉ-IMPLEMENTAÇÃO

Antes de começar **QUALQUER** feature em Phase A:

- [ ] Li o requisito completo
- [ ] Identifiquei arquivos que vão mudar
- [ ] Rodei testes existentes localmente (`pytest tests/`)
- [ ] Criei branch: `feature/ISSUE-XX` ou `fix/ISSUE-XX`
- [ ] Identifiquei escopo exato (o que NÃO fazer)
- [ ] Validei se é **Phase A** ou Phase B/C (não confundir)

**Se não puder marcar TODOS = NÃO COMECE**

***

## ✅ CHECKLIST PÓS-IMPLEMENTAÇÃO

Antes de submeter code para review:

- [ ] Código escrito com **type hints 100%**
- [ ] Docstrings em **classes/funções públicas** (Google-style)
- [ ] Testes unitários criados (**85%+ coverage**)
- [ ] Testes integration criados (se aplicável)
- [ ] Rodou: `pytest tests/ -v --cov=src/grimperium`
- [ ] Rodou: `ruff check src/ tests/` (**ZERO WARNINGS**)
- [ ] Rodou: `black src/ tests/` (sem issues)
- [ ] Rodou: `mypy src/ --strict` (sem errors)
- [ ] Rodou: `pre-commit run --all-files` (todos passam)
- [ ] Sem TODOs, XXX, FIXME ou placeholders
- [ ] Commit message: `type(scope): description`
- [ ] **Pronto para produção** (100% funcional)

**Se NÃO passar em TODOS = NÃO está pronto**

***

## 📁 ESTRUTURA ATUAL DO PROJETO

```
src/grimperium/
├── crest_pm7/                      ← PHASE A (Foco Atual) ⭐
│   ├── pipeline.py                 (Orquestração do workflow)
│   ├── conformer_generator.py       (Gerar conformadores CREST)
│   ├── threshold_monitor.py         (Monitorar alertas de threshold)
│   ├── timeout_predictor.py         (Predizer timeouts - Phase A)
│   ├── validation.py                (Validação de moléculas)
│   ├── mopac_optimizer.py           (Otimização MOPAC)
│   ├── molecule_processor.py        (Processamento de moléculas)
│   ├── result_evaluator.py          (Avaliação de resultados)
│   ├── energy_extractor.py          (Extração de energias)
│   ├── conformer_selector.py        (Seleção de conformadores)
│   ├── config.py                    (Configurações Phase A)
│   └── logging_utils.py             (Logging estruturado)
│
├── core/                            ← Delta Learning + Metrics
│   ├── delta_learning.py
│   └── metrics.py
│
├── data/                            ← Data Loading + Fusion (Phase B)
│   ├── loader.py
│   ├── semiempirical.py
│   └── fusion.py
│
├── models/                          ← ML Models (Phase B - Disponível, não integrado)
│   ├── base_model.py                (BaseModel)
│   ├── kernel_ridge.py              (KernelRidgeModel)
│   ├── xgboost_model.py             (XGBoostModel)
│   └── delta_ensemble.py            (Delta-learning ensemble)
│
├── utils/                           ← Helpers Gerais
│   ├── feature_engineering.py
│   └── (logging/validation em crest_pm7/)
│
├── config.py                        ← Configurações Globais
└── api.py                           ← Interface Pública

tests/
├── unit/                            ← Testes isolados (Phase A focus)
├── integration/                     ← Testes pipeline (Phase A focus)
├── fixtures/                        ← Mock & Real data
│   ├── mock_data.py
│   └── real_data.py
├── experiments/                     ← Testes de hipóteses
└── conftest.py                      ← Fixtures compartilhadas
```

***

## 🚫 NUNCA FAÇA ISTO (CRÍTICO)

### ❌ 1. Testes Incompletos
```python
# ❌ ERRADO
def test_threshold():
    monitor = ThresholdMonitor()
    # Teste vazio - não testa nada
    pass

# ❌ ERRADO
@pytest.mark.skip
def test_alert_triggering():
    # Testes com skip não são aceitáveis
    pass
```

**Correto:**
```python
# ✅ CORRETO
@pytest.mark.unit
def test_alert_triggered():
    monitor = ThresholdMonitor(threshold=0.95)
    alert = monitor.check_threshold(0.98)
    assert alert.triggered is True
```

### ❌ 2. Lógica Não Documentada
```python
# ❌ ERRADO - Sem docstring
def process(data):
    return data.sum()

# ❌ ERRADO - Sem type hints
def calculate_rmsd(a, b):
    return np.sqrt(np.mean((a-b)**2))

# ❌ ERRADO - Sem exemplos
def threshold_alert(value, threshold):
    """Check threshold."""
    return value > threshold
```

**Correto:**
```python
# ✅ CORRETO
def threshold_alert(value: float, threshold: float) -> bool:
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
    if threshold < 0:
        raise ValueError("Threshold cannot be negative")
    return value > threshold
```

### ❌ 3. Ignorar Linting
```bash
# ❌ ERRADO - Commitar com warnings
git commit -m "quick fix"
# (Sem rodar ruff, black, mypy)

# ❌ ERRADO - Ter warnings/errors
$ ruff check src/
src/grimperium/threshold.py:15:1: E501 line too long
$ black --check src/
would reformat src/grimperium/threshold.py
```

**Correto:**
```bash
# ✅ CORRETO - Validar antes de commitar
pytest tests/ -v --cov=src/grimperium
ruff check src/ tests/       # ✅ ZERO WARNINGS
black src/ tests/             # ✅ Sem issues
mypy src/ --strict            # ✅ Sem errors
pre-commit run --all-files    # ✅ Todos passam
git commit -m "feat(crest_pm7): implement timeout predictor"
```

### ❌ 4. Mudar Estrutura
```python
# ❌ ERRADO - Renomear módulo sem aprovação
# src/crest_pm7/threshold_monitor.py → src/threshold.py

# ❌ ERRADO - Criar diretórios fora do plano
# src/grimperium/new_feature_folder/  (sem requisito)

# ❌ ERRADO - Deletar código sem backup
# Remover src/grimperium/models/ (sem mover para .archive/)
```

### ❌ 5. Código Incompleto
```python
# ❌ ERRADO - TODOs no código
def predict_timeout(molecule):
    # TODO: Implement ML model
    return 0.0

# ❌ ERRADO - Placeholders
def optimize():
    # XXX: implement optimization
    pass

# ❌ ERRADO - Funções stub
def compute_features(data):
    """Compute features."""
    return None  # Não implementado
```

### ❌ 6. Type Hints Ausentes
```python
# ❌ ERRADO
def process_molecules(data):
    return data

# ❌ ERRADO
def threshold_alert(value, threshold):
    return value > threshold
```

***

## ✅ SEMPRE FAÇA ISTO (OBRIGATÓRIO)

### ✅ 1. Type Hints Completas (100%)

```python
from typing import Optional, List, Dict, Tuple
import numpy as np
import pandas as pd

def process_molecules(
    data: pd.DataFrame,
    threshold: float = 0.95,
    return_details: bool = False
) -> tuple[np.ndarray, Dict[str, float]]:
    """Process molecules with quality threshold.
    
    Args:
        data: Input molecule data
        threshold: Quality threshold (0-1)
        return_details: Whether to return detail dict
        
    Returns:
        Tuple of processed array and metrics dict
        
    Raises:
        ValueError: If threshold invalid
    """
    if not (0.0 <= threshold <= 1.0):
        raise ValueError("Threshold must be between 0 and 1")
    
    results = np.array([])
    metrics = {"threshold": threshold}
    return results, metrics
```

**Regra:** TODA função/método público = type hints completas

### ✅ 2. Docstrings Google-Style

```python
class ThresholdMonitor:
    """Monitor molecular energy thresholds.
    
    Detects anomalies in CREST PM7 simulations by tracking
    energy values against configurable thresholds. Used in
    Phase A to validate conformer generation quality.
    
    Attributes:
        threshold: Maximum acceptable energy change (kcal/mol)
        window_size: Number of recent values to monitor
        alert_history: List of triggered alerts
        
    Example:
        >>> monitor = ThresholdMonitor(threshold=5.0, window_size=10)
        >>> alert = monitor.check_threshold(6.5)
        >>> assert alert.triggered is True
    """
    
    def __init__(
        self,
        threshold: float = 0.95,
        window_size: int = 10
    ) -> None:
        """Initialize threshold monitor.
        
        Args:
            threshold: Energy threshold value (kcal/mol)
            window_size: Number of recent values to track
            
        Raises:
            ValueError: If threshold or window_size invalid
        """
        if threshold <= 0 or window_size <= 0:
            raise ValueError("Must be positive values")
        self.threshold = threshold
        self.window_size = window_size
    
    def check_threshold(self, value: float) -> bool:
        """Check if value exceeds threshold.
        
        Args:
            value: Current measurement value
            
        Returns:
            True if alert triggered
            
        Example:
            >>> monitor = ThresholdMonitor(threshold=0.95)
            >>> monitor.check_threshold(0.98)
            True
        """
        return value > self.threshold
```

**Regra:** Toda classe e função pública = docstring Google-style com exemplos

### ✅ 3. Testes Estruturados (85%+ Coverage)

```python
# tests/unit/test_threshold_monitor.py

import pytest
from grimperium.crest_pm7 import ThresholdMonitor, AlertType

@pytest.mark.unit
class TestThresholdMonitor:
    """Unit tests for ThresholdMonitor class."""
    
    @pytest.fixture
    def monitor(self) -> ThresholdMonitor:
        """Create test monitor instance."""
        return ThresholdMonitor(threshold=0.95, window_size=10)
    
    def test_initialization(self) -> None:
        """Test valid initialization."""
        monitor = ThresholdMonitor(threshold=1.0, window_size=5)
        assert monitor.threshold == 1.0
        assert monitor.window_size == 5
    
    def test_alert_triggered(self, monitor: ThresholdMonitor) -> None:
        """Test alert triggering when threshold exceeded."""
        alert = monitor.check_threshold(0.98)
        assert alert.triggered is True
        assert alert.alert_type == AlertType.THRESHOLD_VIOLATION
    
    def test_alert_not_triggered(self, monitor: ThresholdMonitor) -> None:
        """Test no alert when within threshold."""
        alert = monitor.check_threshold(0.90)
        assert alert.triggered is False
    
    def test_invalid_threshold(self) -> None:
        """Test ValueError for invalid threshold."""
        with pytest.raises(ValueError, match="Must be positive"):
            ThresholdMonitor(threshold=-1.0)
    
    def test_invalid_window_size(self) -> None:
        """Test ValueError for invalid window_size."""
        with pytest.raises(ValueError, match="Must be positive"):
            ThresholdMonitor(window_size=0)
    
    @pytest.mark.integration
    def test_with_real_molecules(self, monitor: ThresholdMonitor) -> None:
        """Integration test with real molecule data."""
        molecules = load_test_molecules()
        for mol in molecules:
            alert = monitor.check_threshold(mol.energy)
            assert hasattr(alert, 'triggered')
```

**Regra:** 100% novo código = testes. Mínimo 85% coverage.

### ✅ 4. Pre-commit Validation

```bash
# Antes de CADA commit - RUN ISTO

# 1. Testes
pytest tests/ -v --cov=src/grimperium --cov-report=term-missing
# Expected: 100% passed, ≥85% coverage

# 2. Linting
ruff check src/ tests/
# Expected: 0 warnings (ZERO TOLERANCE)

# 3. Formatting
black src/ tests/
# Expected: no changes needed

# 4. Type Checking
mypy src/ --strict
# Expected: 0 errors

# 5. Pre-commit Hooks
pre-commit run --all-files
# Expected: ✅ all hooks pass

# 6. Imports
python -c "from grimperium.crest_pm7 import *; print('✅ OK')"
```

### ✅ 5. Git Workflow Correto

```bash
# 1. Criar branch
git checkout -b feature/threshold-alerts

# 2. Trabalhar
# ... escrever código com type hints, testes, docstrings ...

# 3. Validar tudo
pytest tests/ && ruff check src/ && black src/ && mypy src/

# 4. Pre-commit validation
pre-commit run --all-files

# 5. Commitar
git commit -m "feat(crest_pm7): implement threshold monitor with alerts

- Add ThresholdMonitor class with configurable thresholds
- Implement alert detection for energy anomalies
- Add comprehensive unit tests (89% coverage)
- Add Google-style docstrings with examples
- Type hints: 100% compliance

Fixes #42"

# 6. Push & PR
git push origin feature/threshold-alerts
```

***

## 🎯 PHASE A vs PHASE B/C

### ✅ FAZER em Phase A (Escopo Claro)
- ✅ Validar sistema com moléculas reais (10-20 moléculas)
- ✅ Implementar `threshold_monitor.py` (monitorar violações)
- ✅ Implementar `timeout_predictor.py` (predizer timeouts)
- ✅ Testes unitários completos
- ✅ Integration tests para pipelines
- ✅ Integração com logging estruturado
- ✅ Documentação clara e exemplos
- ✅ Establish baseline metrics

### ❌ NÃO FAZER em Phase A
- ❌ Data loading massivo (→ Phase B)
- ❌ Model training (→ Phase B/C)
- ❌ MOPAC integration avançada (→ Phase B)
- ❌ Feature engineering complexa (→ Phase C)
- ❌ Otimização de performance (→ Phase C)
- ❌ Refactors não relacionados (→ Phase C)

**DÚVIDA?** Sempre pergunte: **"É Phase A ou Phase B/C?"**

***

## 📊 VALIDAÇÕES FINAIS (PRÉ-SUBMISSÃO)

### Run Antes de Afirmar "Pronto"

```bash
# 1. TESTES (Coverage ≥85%)
pytest tests/ -v --cov=src/grimperium --cov-report=term-missing
# Expected: ✅ 100% passed, coverage ≥85%

# 2. LINTING (Zero Warnings)
ruff check src/ tests/
# Expected: ✅ 0 warnings

# 3. FORMATTING
black --check src/ tests/
# Expected: ✅ All formatted correctly

# 4. TYPE CHECKING (Strict)
mypy src/ --strict
# Expected: ✅ 0 errors

# 5. PRE-COMMIT HOOKS
pre-commit run --all-files
# Expected: ✅ All hooks pass

# 6. IMPORTS
python -c "from grimperium.crest_pm7 import *; print('✅ Imports OK')"
# Expected: ✅ Imports OK

# 7. QUICK TEST (Phase A)
python scripts/phase_a_quick_test.py
# Expected: ✅ PHASE A TEST PASSED
```

**Se QUALQUER um falhar = NÃO está pronto para produção**

***

## 💬 PERGUNTAS FREQUENTES

**P: Posso usar biblioteca X?**  
R: Só se estiver em `pyproject.toml`. Não adicione deps sem pedir.

**P: Qual convention de nomes?**  
R: `snake_case` (funções/variáveis), `PascalCase` (classes), `UPPER_CASE` (constantes)

**P: Quantos testes preciso?**  
R: Mínimo 85% coverage. Sempre aim para 95%+.

**P: E se encontrar código antigo ruim?**  
R: Só refactor se for Phase A relevante. Deprecated fica em `.archive/`.

**P: Posso usar async/await?**  
R: Sim, com `asyncio`. Adicione type hints completas.

**P: Como estruturo novo módulo?**  
R: Veja `src/grimperium/crest_pm7/threshold_monitor.py` como exemplo.

**P: Qual é o line-width para Black?**  
R: 88 caracteres (default Black).

***

## 🚀 WORKFLOW RECOMENDADO (Step-by-Step)

### 1. Receba Requisito
```
"Implementar timeout predictor para Phase A"
```

### 2. Valide com Checklists
```
PRÉ-IMPLEMENTAÇÃO:
- [ ] Li requisito completo
- [ ] Identifiquei arquivos afetados
- [ ] Rodei testes atuais
- [ ] Criei branch feature/timeout-predictor
- [ ] É Phase A? (Sim ✅)
```

### 3. Escreva Testes Primeiro (TDD)
```python
# tests/unit/test_timeout_predictor.py
def test_timeout_prediction():
    predictor = TimeoutPredictor()
    mol = load_test_molecule("benzene")
    timeout_pred = predictor.predict(mol)
    
    assert isinstance(timeout_pred, float)
    assert timeout_pred > 0
    assert timeout_pred < 3600
```

### 4. Implemente com Type Hints + Docstrings
```python
# src/grimperium/crest_pm7/timeout_predictor.py
class TimeoutPredictor:
    """Predict CREST PM7 simulation timeout."""
    
    def predict(self, molecule: Molecule) -> float:
        """Estimate timeout based on molecule properties.
        
        Args:
            molecule: Input molecule
            
        Returns:
            Estimated timeout in seconds
        """
        # Implementação...
        pass
```

### 5. Valide com Linters
```bash
pytest tests/unit/test_timeout_predictor.py -v
ruff check src/grimperium/crest_pm7/
black src/grimperium/crest_pm7/
mypy src/grimperium/crest_pm7/ --strict
pre-commit run --all-files
```

### 6. Commit Semântico
```bash
git commit -m "feat(crest_pm7): implement timeout predictor

- Add TimeoutPredictor class for CREST PM7 simulations
- Based on molecule size and complexity metrics
- 92% test coverage with comprehensive unit tests
- Full type hints and Google-style docstrings

Implements Phase A requirement #XX"
```

### 7. Code Review + Merge
```
Feature implementado, testado, documentado, committed.
Pronto para code review e merge.
```

***

## 📚 RECURSOS & DOCUMENTAÇÃO

- **Este Arquivo:** `docs/CLAUDE.md` (seu guia comportamental)
- **Phase A Start:** `docs/PHASE-A-START-HERE.md` (como executar)
- **Complete Guide:** `docs/phase-a-complete-guide.md` (detalhado)
- **Resultados:** `docs/PHASE-A-RESULTS.md` (template resultados)
- **Repository:** `~/Projetos/grimperium_V2/`

***

## ✨ RESUMO EXECUTIVO

**Você é production-ready quando:**

✅ **Código** escrito com **type hints 100%**  
✅ **Testes** verdes com **85%+ coverage**  
✅ **Linting** passa: **ruff, black, mypy** (ZERO WARNINGS)  
✅ **Documentado**: **docstrings Google-style com exemplos**  
✅ **Pronto para produção**: **sem TODOs, placeholders ou incompletos**  
✅ **Git log limpo**: **commits semânticos com contexto**  
✅ **Pre-commit**: **todos os hooks passam**  

**Se QUALQUER item falhar = NÃO está ready**

***

**Last Review:** January 10, 2026 | **Status:** Production-Ready ✅  
**Merge Status:** Atual + Global Combined ✅

[1](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/3ea1e28f-e134-4392-84c0-eb3217d4cba9/02_PROMPT_001_SCAFFOLDING_INICIAL.md)
[2](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/92bc0bd6-b268-4388-a118-49e98f7041e6/03_RESUMO_EXECUTIVO.md)
[3](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/d041d8d9-fec4-4818-9652-ce8f8d916109/04_INSTRUCOES_CLAUDE_CODE.md)
[4](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/281197f5-9d70-426c-b952-849be2abf0b7/05_DECISIONS_FINAL.md)
[5](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/ae710f1f-5127-4423-8e27-ee45400b6f32/06_DATASET_DELTA_STRATEGY.md)
[6](https://ppl-ai-file-upload.s3.amazonaws.com/connectors/google_drive/1Rc7-BF3g7EiHhBESsriU1HnnfbmwS7A6U5Ph4GThWyg/d602aac7-f074-4c68-9b0c-d8e437273dce/grimperium_cbs_opt.xlsx)
[7](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/64ed7830-0072-44ca-984c-5863eb996050/ci.md)
[8](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/daf80ce0-7227-460f-83da-ccb20ab6f4e8/claude_code_instructions.md)
[9](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/collection_49d29401-6295-471c-9e0c-87d50f288c1c/8b09c354-8bf4-4748-b897-90cf038e47e7/RELATORIO_TECNICO_GRIMPERIUM-1.md)
