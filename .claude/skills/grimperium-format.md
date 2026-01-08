---
name: grimperium-format
description: Formata código com Black e valida com Ruff em um comando
tools: [bash]
context: fork
user-invocable: true
allowed-tools:
  - Bash(black *)
  - Bash(ruff *)
---

# 🎨 Skill: Format + Lint Code (One Command)

**Propósito:** Executar Black + Ruff para limpar código em um comando

**Quando usar:** Antes de commitar, quer ter certeza que code está limpo

## O que esta skill faz

1. ✅ Executa `black src/ tests/` (reformat)
2. ✅ Executa `ruff check src/ tests/` (lint)
3. ✅ Exibe erros restantes (se houver)
4. ✅ Tudo em um comando

## Como Usar

### Formato Padrão

```
@claude /grimperium-format
```

### Formato + Checar (sem reformat)

```
@claude /grimperium-format --check-only
```

### Formato + Arquivo Específico

```
@claude /grimperium-format src/grimperium/core/metrics.py
```

## Output

```
Processing files...

1️⃣ BLACK FORMATTING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
reformatted src/grimperium/core/metrics.py
reformatted src/grimperium/data/loader.py
reformatted tests/unit/test_models.py

All done! 3 files reformatted.

2️⃣ RUFF LINTING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
✅ No lint errors found!

📊 SUMMARY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Black:    3 files reformatted ✅
Ruff:     0 errors ✅
Status:   READY TO COMMIT ✅
```

## Workflow Recomendado

```
Antes de fazer commit:

1. @claude /grimperium-format
   └─ Formata + linta código

2. Verifica output
   └─ Se houver erros, corrige manualmente

3. git add + commit
   └─ Código está clean ✅
```

## Se Houver Erros

```
2️⃣ RUFF LINTING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
src/grimperium/core/metrics.py:45:1: F401 'numpy' imported but unused
src/grimperium/data/loader.py:12:5: E501 line too long (120 > 88)

❌ 2 errors found

📊 SUMMARY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Black:    3 files reformatted ✅
Ruff:     2 errors ❌
Status:   FIX REQUIRED ❌

Next: Fix errors manually or use @claude to fix these ruff errors
```

## Notas

- Black sempre reformata automaticamente
- Ruff só valida (não muda código)
- Se Ruff encontrar erros, você precisa corrigir manualmente ou pedir ajuda
- Roda em background (context: fork) para não bloquear

## Integração com Git

```bash
# Workflow completo
@claude /grimperium-format    # Formata + linta
git add -A                    # Adiciona mudanças
git commit -m "feat: ..."     # Commita
git push                      # Publica
```

## Commands Disponíveis

| Command | O Que Faz |
|---------|-----------|
| `/grimperium-format` | Formata + linta tudo |
| `/grimperium-format --check-only` | Só valida (não muda) |
| `/grimperium-format src/` | Formata só src/ |
| `/grimperium-format tests/` | Formata só tests/ |
