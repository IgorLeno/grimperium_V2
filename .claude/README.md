# Claude Code Configuration & Procedures

Este diretório contém configurações e procedimentos para desenvolvimento com Claude Code.

## Arquivos

### `settings.json`
Configurações do Claude Code com wildcard permissions para comandos frequentes (poetry, pytest, black, ruff, mypy, git).

### `settings.local.json`
Configurações locais do Claude Code (gitignored).

### `skills/`
Skills customizadas para automação de workflows do Grimperium:
- `grimperium-ci.md` - Fix automático de erros CI/CD
- `grimperium-tests.md` - Testes em background (fork)
- `grimperium-format.md` - Formatação + linting em um comando

### `DEVELOPMENT_PROCEDURES.md` ⭐
**Procedimentos obrigatórios de desenvolvimento.**

**LEIA ANTES DE CADA BATCH!**

Contém:
- Fluxo de trabalho padrão
- **Checklist de CHANGELOG (OBRIGATÓRIO)** 📝
- Guidelines de testes
- Formato de commits
- Padrões de documentação

## ⚠️ Novo Hábito Estabelecido em BATCH 3

**SEMPRE atualizar CHANGELOG.md ao final de cada batch!**

Sem exceções. Este é agora um passo obrigatório do workflow.

Ver `DEVELOPMENT_PROCEDURES.md` seção "CHANGELOG Update" para detalhes.

## Quick Start

Antes de iniciar um novo batch:

```bash
# 1. Ler procedimentos
cat .claude/DEVELOPMENT_PROCEDURES.md

# 2. Durante o trabalho - usar TodoWrite
# 3. Ao final - ATUALIZAR CHANGELOG
vim CHANGELOG.md

# 4. Commit
git commit -m "feat(scope): description"
```

## 🎯 Skills Disponíveis (Claude Code 2.1.0)

### 1. grimperium-ci-fix

**O que faz:** Analisa CI Error Report e corrige automaticamente

**Como usar:**
```
@claude /grimperium-ci-fix

[cola CI_ERROR_SUMMARY.md aqui]
```

**Benefício:** Automático fix de lint + type + test errors

---

### 2. grimperium-tests

**O que faz:** Roda testes completos em background

**Como usar:**
```
@claude /grimperium-tests
@claude /grimperium-tests --with-html
```

**Benefício:** Testes rodam paralelo, você continua codificando (context: fork)

---

### 3. grimperium-format

**O que faz:** Formata com Black + linta com Ruff

**Como usar:**
```
@claude /grimperium-format
@claude /grimperium-format src/grimperium/core/
```

**Benefício:** Code sempre limpo antes de commitar

## 🚀 Workflow Recomendado

### Durante Desenvolvimento

```
1. Você está codificando
2. @claude /grimperium-format (quick check)
3. @claude /grimperium-tests (em background)
   - Você continua codificando
   - Testes rodam paralelo
   - Notificação quando terminar
```

### Antes de Commitar

```
1. @claude /grimperium-format
   └─ Formata + linta tudo
2. @claude /grimperium-tests --with-html
   └─ Roda testes com coverage
3. Visualiza htmlcov/index.html
   └─ Verifica coverage
4. git add + commit se tudo OK
```

### Se CI Falhar

```
1. GitHub Actions reporta erro
2. Download CI_ERROR_SUMMARY.md
3. @claude /grimperium-ci-fix < CI_ERROR_SUMMARY.md
4. Claude Code corrige tudo automaticamente
5. git push (novo commit)
6. CI passa ✅
```

## 📊 Features Claude Code 2.1.0

Essas skills aproveitam:

- ✅ **Automatic Skill Hot-reload** — Skills aparecem sem reiniciar
- ✅ **Agent Forking (context: fork)** — Testes rodam em background
- ✅ **Wildcard Bash Permissions** — Menos prompts de segurança

## Estrutura de um Batch

1. **Plan Mode** → Escrever plano
2. **Implementação** → Seguir plano com TodoWrite
3. **Testes** → Verificar que tudo passa (use `/grimperium-tests`)
4. **📝 CHANGELOG** → **OBRIGATÓRIO!**
5. **Commit** → Finalizar batch (use `/grimperium-format` antes)

---

**Última atualização:** 2026-01-07 (Claude Code 2.1.0 Skills Integration)
