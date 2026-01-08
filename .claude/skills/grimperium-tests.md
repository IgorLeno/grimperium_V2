---
name: grimperium-tests
description: Roda teste completa em background sem bloquear seu desenvolvimento
tools: [bash]
context: fork
user-invocable: true
allowed-tools:
  - Bash(poetry *)
  - Bash(pytest *)
---

# 🧪 Skill: Run Full Test Suite (Background)

**Propósito:** Executar testes completos em background enquanto você desenvolve

**Quando usar:** Quando você quer validar seu código mas não quer esperar

## O que esta skill faz

1. ✅ Instala dependências (se necessário)
2. ✅ Roda pytest com coverage
3. ✅ Gera relatório HTML de coverage
4. ✅ Exibe summary de resultados
5. ✅ Roda em background (context: fork)

## Como Usar

### Básico: Rodar Testes Simples

```
@claude /grimperium-tests
```

Resultado:
- Terminal output com pass/fail
- Summary: "88 passed, 11 xfailed, 20 skipped"
- Você continua desenvolvendo sem esperar

### Avançado: Com Coverage Detalhado

```
@claude /grimperium-tests --with-html

# Gera:
# - Terminal report
# - htmlcov/index.html com coverage por arquivo
# - Relatório de linhas não cobertas
```

## Workflow Recomendado

```
┌─────────────────────────────────────────────┐
│ Você está desenvolvendo (editing code)      │
├─────────────────────────────────────────────┤
│ Você: @claude /grimperium-tests             │
│ ↓                                           │
│ Claude: Inicia testes em background (fork)  │
│ ↓                                           │
│ Você: Continua codificando (não bloqueado)  │
│ ↓                                           │
│ [5 minutos depois]                          │
│ ↓                                           │
│ Notificação: "Tests completed! 88 passed"   │
│ ↓                                           │
│ Você: Clica em htmlcov/index.html           │
│ ↓                                           │
│ Vê coverage por arquivo                     │
└─────────────────────────────────────────────┘
```

## Commands Disponíveis

| Command | O Que Faz |
|---------|-----------|
| `/grimperium-tests` | Teste rápido com summary |
| `/grimperium-tests --with-html` | Teste + coverage HTML |
| `/grimperium-tests --verbose` | Output detalhado |
| `/grimperium-tests src/` | Só testa módulo específico |

## Output

```
============================= test session starts ==============================
platform linux -- Python 3.11.x, pytest-7.4.4
rootdir: /home/user/grimperium, configfile: pyproject.toml
collected 119 items

tests/unit/test_loader.py ..................                          [ 15%]
tests/unit/test_fusion.py ....................                        [ 30%]
tests/integration/test_pipeline.py ................                   [ 45%]
tests/experiments/test_validate_hypothesis.py ....................... [ 100%]

---------- coverage: platform linux, python 3.11.x -----------
Name                          Stmts   Miss  Cover
────────────────────────────────────────────────
src/grimperium/__init__.py       12      2    83%
src/grimperium/core/metrics.py   145      8    94%
src/grimperium/data/loader.py    203     15    93%
...
────────────────────────────────────────────────
TOTAL                          1842     87    95%

Coverage HTML written to dir `htmlcov`

=============================== 88 passed in 5.23s ==============================
```

## Notas Importantes

- Roda em **background** (context: fork) — você não espera
- Se houver falhas, você recebe notificação
- Coverage report em HTML fica em `htmlcov/index.html`
- Pode opener coverage no navegador para análise detalhada

## Se Testes Falharem

Você recebe:
```
❌ Tests failed: 2 failed, 86 passed

Use @claude analyze-test-failures para detalhes
```

Então você pede análise e fix automático.
