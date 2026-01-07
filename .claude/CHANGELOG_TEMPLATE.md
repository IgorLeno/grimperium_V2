# CHANGELOG Entry Template

Use este template ao final de cada batch para atualizar o CHANGELOG.md.

## Template de Entrada

```markdown
### Added
- **BATCH X: [Título Descritivo]** (YYYY-MM-DD)
  - `path/to/new_file.py`: Descrição da funcionalidade
    - `function_name()`: Descrição específica
    - `class_name`: Descrição específica
  - `tests/test_new_feature.py`: Testes cobrindo X, Y, Z
  - Documentação detalhada em docstrings

### Changed
- **BATCH X: [Título da Mudança]** (YYYY-MM-DD)
  - Substituído X por Y para melhorar Z
  - Refatorado `module.py` para usar pattern ABC
  - Atualizado interface de `function()` para aceitar novo parâmetro

### Fixed
- **BATCH X: [Problema Corrigido]** (YYYY-MM-DD)
  - Corrigido bug X que causava Y em `file.py:123`
  - Resolvido problema de Z quando condição W
  - Fixed distribution shift artifact (RMSE=1008 → 61)

### Deprecated
- **BATCH X** (YYYY-MM-DD)
  - `old_function()`: Deprecado em favor de `new_function()`
    - Reason: [Explicação clara do motivo]
    - Replacement: Use `new_function()` from `module.py`
    - Migration: [Passos para migração]
    - Warning: DeprecationWarning adicionado

### Removed
- **BATCH X** (YYYY-MM-DD)
  - Removido `obsolete_module.py` (depreciado em BATCH Y)
  - Deletado testes obsoletos de `old_feature`

### Test Results - BATCH X
- **[Nome do Teste]**
  - Métrica 1: valor ✓
  - Métrica 2: valor ✓
  - Resultado: PASS/FAIL ✅/❌
```

## Exemplo Real (BATCH 3)

```markdown
### Added
- **BATCH 3: Hypothesis Validation Test Suite** (2026-01-07)
  - `tests/experiments/conftest.py`: New fixtures with filtered and extreme data
    - `real_data_1k_filtered()`: Realistic distribution [-1000, +1000] kcal/mol (99.1% of data)
    - `real_data_1k_extreme()`: Pathological distribution for stress testing (includes outliers)
    - `synthetic_data_1k()`: Fast synthetic data for CI/fallback tests
  - `tests/experiments/test_validate_hypothesis.py`: Main hypothesis validation tests
    - `test_decision_gate_delta_vs_direct()`: Primary test with filtered data
    - `test_synthetic_fallback()`: Synthetic data fallback test
  - `tests/experiments/test_stress_distribution_shift.py`: Robustness stress tests
    - `test_stress_distribution_shift_extreme()`: Tests with severe distribution shift
    - `test_distribution_shift_detection()`: Validates distribution shift detection
  - Comprehensive documentation of methodological decisions in test docstrings

### Changed
- **BATCH 3: Fixture Methodology Correction** (2026-01-07)
  - Replaced unfiltered data approach with filtered realistic distribution
  - Split validation testing (OPTION B) from stress testing (OPTION A)
  - Updated mock PM7 generator with configurable magnitude bias
  - Improved fixture logging and statistics reporting

### Fixed
- **BATCH 3: Distribution Shift Artifact** (2026-01-07)
  - Fixed misleading RMSE=1008 caused by severe distribution shift in unfiltered data
  - Corrected hypothesis validation to use realistic data regime (std~70, not 7230)
  - Resolved train/test distribution mismatch (6.1 vs 615 kcal/mol mean difference)
  - Fixed Direct model comparison to be fair (61.11 RMSE vs 1008.88 artifact)

### Deprecated
- **BATCH 3** (2026-01-07)
  - `tests/fixtures/conftest.py::real_data_1k()`: Now deprecated with warning
    - Reason: Uses unfiltered data causing distribution shift artifacts
    - Replacement: Use `real_data_1k_filtered` or `real_data_1k_extreme` from `tests/experiments/conftest.py`
    - Warning: Added DeprecationWarning to guide users to new fixtures

### Test Results - BATCH 3
- **Hypothesis Validation (Realistic Regime)**
  - Filter: [-1000, +1000] kcal/mol (removes 0.9% outliers)
  - Distribution shift: 6.1 kcal/mol (minimal)
  - RMSE Delta: 9.31 kcal/mol ✓
  - RMSE Direct: 61.11 kcal/mol ✓
  - Improvement: 6.6x (84.8%)
  - R² Delta: 0.9768
  - **DECISION GATE: PASS** ✅

- **Stress Test (Extreme Regime)**
  - Unfiltered data (outliers up to -325407 kcal/mol)
  - Distribution shift: 615.3 kcal/mol (severe)
  - RMSE Delta: 13.83 kcal/mol (robust)
  - RMSE Direct: 1008.88 kcal/mol (catastrophic failure expected)
  - Robustness ratio: 73x
  - **STRESS TEST: PASS** ✅
```

## Checklist de Validação

Antes de finalizar a entrada do CHANGELOG:

- [ ] Todas as mudanças estão categorizadas?
- [ ] Data está no formato correto `(YYYY-MM-DD)`?
- [ ] Paths de arquivo estão corretos e completos?
- [ ] Funções/classes referenciadas têm backticks?
- [ ] Resultados de testes incluídos (se aplicável)?
- [ ] Deprecations têm replacement documentado?
- [ ] Fixes têm contexto suficiente (file:line, condição)?
- [ ] Seção está na ordem correta (Added → Changed → Fixed → Deprecated → Removed)?

## Dicas

1. **Seja específico**: "Fixed bug" → "Fixed RMSE=1008 caused by distribution shift"
2. **Referencie código**: Sempre use `backticks` para código
3. **Inclua métricas**: Se há testes, mostre resultados
4. **Documente decisões**: Especialmente em "Changed" e "Fixed"
5. **Migration path**: Sempre para Deprecated/Removed

---

**Última atualização:** 2026-01-07 (BATCH 3)
