---
name: grimperium-ci-fix
description: Analisa CI Error Summary Report e corrige todos os erros automaticamente
tools: [bash, file]
context: fork
user-invocable: true
allowed-tools:
  - Bash(poetry *)
  - Bash(black *)
  - Bash(ruff *)
  - Bash(mypy *)
  - Bash(git *)
---

# 🔧 Skill: Fix Grimperium CI Errors

**Propósito:** Automatizar correção de erros de CI/CD

**Quando usar:** Quando você recebe um CI Error Summary Report do GitHub Actions

## O que esta skill faz

1. ✅ Parse do CI Error Report
2. ✅ Identifica tipo de erro (Lint, Type, Tests)
3. ✅ Executa fixes automáticos
4. ✅ Valida que erros foram corrigidos
5. ✅ Commit + Push para GitHub

## Tipos de Erros Suportados

### Lint Errors (Black Format)

Arquivos desformatados são rewritados automaticamente:

```bash
poetry run black src/ tests/
```

### Type Errors (Mypy)

Erros de tipo em return types são corrigidos:

```bash
poetry run mypy src/grimperium --ignore-missing-imports
```

Depois, arquivo afetado é editado manualmente se necessário.

### Test Errors (Pytest)

Testes quebrados são analisados:

```bash
poetry run pytest tests/ -v --tb=short
```

## Como Usar

### Opção 1: Colar CI Error Summary

```
@claude /grimperium-ci-fix

Aqui está o CI Error Report:

# CI/CD Error Summary Report

Generated: 2026-01-07 09:38:30 UTC
Commit: `738992bf4f69`
Branch: `main`
Run: #15

***

## Overall Status

❌ FAILURES DETECTED

[... cola o resto do relatório ...]
```

### Opção 2: Ler arquivo diretamente

```
@claude /grimperium-ci-fix CI_ERROR_SUMMARY.md
```

## Processo Passo a Passo

```
1️⃣ Parse Error Report
   ├─ Lint errors (Black)? → black src/ tests/
   ├─ Type errors (Mypy)? → Listar e corrigir
   └─ Test errors? → Analisar traceback

2️⃣ Execute Fixes
   ├─ poetry run black src/ tests/
   ├─ poetry run mypy src/ --ignore-missing-imports
   └─ poetry run pytest tests/ -v

3️⃣ Validate
   ├─ Ruff: ruff check src/ tests/
   ├─ Mypy: mypy src/grimperium
   └─ Pytest: 88 passed, 0 errors

4️⃣ Commit + Push
   ├─ git add -A
   ├─ git commit -m "Fix: Resolve CI errors - format + type + tests"
   └─ git push origin main
```

## Output

- ✅ Lista de todos os erros corrigidos
- ✅ Comando de cada fix executado
- ✅ Resultado final (pass/fail)
- ✅ Commit hash se bem-sucedido

## Notas

- Se houver erros que exigem decisão humana, eu aviso
- Type errors podem exigir edição manual de signature de função
- Testes pode requerer mudança de lógica (não apenas formatting)

## Next Steps Quando Falhar

Se a skill não conseguir corrigir:

```bash
# Você executa manualmente
poetry run black src/ tests/
poetry run mypy src/grimperium --ignore-missing-imports

# Então pede ajuda
@claude analyze my errors [copy error messages]
```
