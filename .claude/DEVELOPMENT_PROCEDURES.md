# Grimperium Development Procedures

Este documento estabelece os procedimentos padrão de desenvolvimento para o projeto Grimperium quando trabalhando com Claude Code.

## Índice

1. [Fluxo de Trabalho por Batch](#fluxo-de-trabalho-por-batch)
2. [Atualização do CHANGELOG](#atualização-do-changelog)
3. [Testes e Validação](#testes-e-validação)
4. [Commits e Versionamento](#commits-e-versionamento)
5. [Documentação](#documentação)

---

## Fluxo de Trabalho por Batch

Cada "batch" de trabalho deve seguir este fluxo:

### 1. Planejamento (Plan Mode)
- [ ] Ler especificações/requisitos
- [ ] Explorar codebase relevante
- [ ] Escrever plano de implementação em arquivo dedicado
- [ ] Obter aprovação do plano antes de implementar

### 2. Implementação
- [ ] Seguir plano aprovado
- [ ] Usar test-driven development quando apropriado
- [ ] Documentar decisões críticas nos docstrings
- [ ] Usar TodoWrite para tracking de progresso

### 3. Testes
- [ ] Escrever testes antes ou durante implementação
- [ ] Executar todos os testes relacionados
- [ ] Verificar coverage (se aplicável)
- [ ] Executar testes de integração

### 4. Documentação
- [ ] Atualizar docstrings
- [ ] Adicionar comentários onde lógica não é óbvia
- [ ] Atualizar documentação de arquitetura (se mudou)
- [ ] Atualizar guides/tutoriais (se necessário)

### 5. **CHANGELOG Update (OBRIGATÓRIO)** 📝
**⚠️ NOVO HÁBITO - SEMPRE EXECUTAR AO FINAL DE CADA BATCH**

- [ ] Abrir `CHANGELOG.md`
- [ ] Adicionar entradas na seção `[Unreleased]`
- [ ] Categorizar mudanças:
  - **Added**: Novas funcionalidades, arquivos, testes
  - **Changed**: Mudanças em funcionalidades existentes
  - **Fixed**: Correções de bugs, problemas
  - **Deprecated**: Funcionalidades marcadas para remoção
  - **Removed**: Funcionalidades removidas
  - **Security**: Correções de segurança
- [ ] Incluir data no formato `(YYYY-MM-DD)`
- [ ] Referenciar arquivos/funções modificados
- [ ] Incluir resultados de testes (se relevante)

**Exemplo de entrada:**
```markdown
### Added
- **BATCH X: [Título]** (2026-01-07)
  - `path/to/file.py`: Descrição da funcionalidade
  - `tests/test_feature.py`: Testes cobrindo X, Y, Z
  - Documentação detalhada em docstrings

### Fixed
- **BATCH X** (2026-01-07)
  - Corrigido problema X causando Y
  - Resolvido bug Z em `module.py:123`
```

### 6. Review & Commit
- [ ] Revisar código implementado
- [ ] Executar linters (ruff, mypy)
- [ ] Criar commit com mensagem descritiva
- [ ] Referenciar BATCH no commit message

---

## Atualização do CHANGELOG

### Formato Padrão

Seguimos [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

### Seções Principais

#### Added (Adições)
Novas funcionalidades, arquivos, módulos, testes.

```markdown
### Added
- **BATCH 3: Hypothesis Validation Test Suite** (2026-01-07)
  - `tests/experiments/conftest.py`: Fixtures com dados filtrados
    - `real_data_1k_filtered()`: Distribuição realista
    - `real_data_1k_extreme()`: Dados extremos para stress test
```

#### Changed (Mudanças)
Modificações em funcionalidades existentes, refatorações.

```markdown
### Changed
- **BATCH 3: Fixture Methodology** (2026-01-07)
  - Substituído approach não-filtrado por distribuição realista
  - Split entre testes de validação e stress tests
```

#### Fixed (Correções)
Bugs corrigidos, problemas resolvidos.

```markdown
### Fixed
- **BATCH 3: Distribution Shift Artifact** (2026-01-07)
  - Corrigido RMSE=1008 causado por distribution shift
  - Resolvido mismatch entre treino/teste
```

#### Deprecated (Depreciado)
Funcionalidades marcadas para remoção futura.

```markdown
### Deprecated
- **BATCH 3** (2026-01-07)
  - `tests/fixtures/conftest.py::real_data_1k()`: Deprecado
    - Reason: Distribution shift artifacts
    - Replacement: `real_data_1k_filtered`
```

### Quando Atualizar

**SEMPRE ao final de cada batch**, antes do commit final.

### Checklist de CHANGELOG

- [ ] Seção `[Unreleased]` atualizada
- [ ] Todas as mudanças categorizadas corretamente
- [ ] Data incluída no formato `(YYYY-MM-DD)`
- [ ] Arquivos/funções referenciados com path completo
- [ ] Resultados de testes incluídos (se aplicável)
- [ ] Decisões metodológicas documentadas (se relevante)
- [ ] Deprecations têm motivo e replacement documentados

---

## Testes e Validação

### Tipos de Testes

1. **Unit Tests** (`tests/unit/`)
   - Testes isolados de funções/métodos
   - Mocks para dependências externas
   - Rápidos (<1s cada)

2. **Integration Tests** (`tests/integration/`)
   - Testes de módulos integrados
   - Dados reais ou realistas
   - Moderados (1-5s cada)

3. **Experiment Tests** (`tests/experiments/`)
   - Validação de hipóteses científicas
   - Comparações de modelos
   - Podem ser lentos (30s-2min)

### Executando Testes

```bash
# Todos os testes
pytest

# Batch específico
pytest tests/experiments/test_validate_hypothesis.py -v -s

# Com coverage
pytest --cov=grimperium --cov-report=html

# Apenas testes rápidos
pytest -m "not slow"
```

### Critérios de Aceitação

- [ ] Todos os testes passam
- [ ] Coverage não diminuiu (meta: >80%)
- [ ] Nenhum warning crítico
- [ ] Performance não regrediu significativamente

---

## Commits e Versionamento

### Formato de Commit Message

```
type(scope): short description

Longer description if needed.

- Bullet point 1
- Bullet point 2

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
```

**Types:**
- `feat`: Nova funcionalidade
- `fix`: Correção de bug
- `refactor`: Refatoração de código
- `test`: Adição/modificação de testes
- `docs`: Documentação
- `chore`: Manutenção, config, etc.

**Scopes:**
- `data`: Data loading/processing
- `models`: Model implementations
- `core`: Core algorithms
- `tests`: Test suite
- `experiments`: Hypothesis validation

### Exemplo

```
feat(tests): add BATCH 3 hypothesis validation suite

Implemented comprehensive test suite for delta-learning validation:
- Realistic data regime (filtered [-1000, +1000] kcal/mol)
- Extreme data stress test (unfiltered, severe distribution shift)
- Proper separation of validation vs robustness testing

Results:
- RMSE Delta: 9.31 (realistic) vs 13.83 (extreme)
- RMSE Direct: 61.11 (realistic) vs 1008.88 (extreme)
- Decision Gate: PASS ✅

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
```

---

## Documentação

### Docstrings

Seguir formato Google/NumPy:

```python
def function(arg1: type, arg2: type) -> return_type:
    """
    Short description.

    Longer description if needed, explaining:
    - What the function does
    - Why certain decisions were made
    - Important caveats or edge cases

    Args:
        arg1: Description of arg1
        arg2: Description of arg2

    Returns:
        Description of return value

    Raises:
        ExceptionType: When this happens

    Example:
        >>> result = function(value1, value2)
        >>> assert result == expected
    """
```

### Comentários

- Prefira código auto-explicativo
- Use comentários para **WHY**, não **WHAT**
- Marque TODOs com `# TODO(batch-X): description`
- Marque decisões críticas com `# DESIGN DECISION:`

### Arquitetura

Atualizar `docs/architecture.md` quando:
- Adicionar novos módulos
- Mudar fluxo de dados
- Modificar interfaces públicas

---

## Quick Reference Card

### Final de Batch - Checklist Completo

1. ✅ Todos os testes passam
2. ✅ Linters executados (ruff, mypy)
3. ✅ **CHANGELOG.md atualizado** 📝
4. ✅ Commit criado com mensagem descritiva
5. ✅ Push para repositório (se aplicável)

### Comandos Rápidos

```bash
# Executar testes
pytest tests/experiments/ -v

# Linter
ruff check .
mypy src/

# Atualizar CHANGELOG (manual)
vim CHANGELOG.md  # Adicionar entradas na seção [Unreleased]

# Commit
git add .
git commit -m "feat(scope): description"
```

---

## Notas Importantes

### ⚠️ Sempre Lembrar

1. **CHANGELOG é obrigatório** ao final de cada batch
2. **Testes antes de commit** - nunca commitar código quebrado
3. **Documentação junto com código** - não deixar para depois
4. **Deprecations precisam de migration path** - documentar replacement

### 🎯 Meta de Qualidade

- Coverage: >80%
- Linter warnings: 0
- Type errors: 0
- Documentação: 100% das APIs públicas
- CHANGELOG: Sempre atualizado

---

**Última atualização:** 2026-01-07
**Versão:** 1.0
**Estabelecido em:** BATCH 3
