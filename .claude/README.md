# Claude Code Configuration & Procedures

Este diretório contém configurações e procedimentos para desenvolvimento com Claude Code.

## Arquivos

### `settings.local.json`
Configurações locais do Claude Code (gitignored).

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

## Estrutura de um Batch

1. **Plan Mode** → Escrever plano
2. **Implementação** → Seguir plano com TodoWrite
3. **Testes** → Verificar que tudo passa
4. **📝 CHANGELOG** → **OBRIGATÓRIO!**
5. **Commit** → Finalizar batch

---

**Última atualização:** 2026-01-07 (BATCH 3)
