# Plano: Auditoria e Otimização do Workflow Claude Code

## Contexto

O projeto Grimperium tem uma estrutura `.claude/` funcional mas minimalista: 4 skills (ci, tests, format, docs), `settings.json` com apenas `plansDirectory`, e `settings.local.json` com permissões. O doc de auditoria identificou 6 melhorias: hooks automáticos, skills especializadas, agente científico, settings expandido e `.claudeignore`.

A integração dos hooks usará os eventos reais do Claude Code (`SessionStart`, `Stop`, `PreToolUse`) — não os nomes inventados no doc (`prePlan`, `postImplement`).

---

## Fase 1: Hooks Shell Scripts

Criar diretório `.claude/hooks/` com 3 scripts executáveis:

### `.claude/hooks/pre-plan.sh`
- Verifica existência de `pyproject.toml` (confirma estar na raiz)
- Detecta se `$VIRTUAL_ENV` está ativo
- Exibe branch atual e contagem de arquivos modificados
- **Exit 0 sempre** (informativo, não bloqueante)

### `.claude/hooks/post-implement.sh`
- Executa `ruff check src/ tests/ --quiet`
- Executa `black --check src/ tests/ --quiet`
- Executa `mypy src/grimperium --ignore-missing-imports`
- Executa `pytest tests/unit/ -q --tb=no`
- Ruff falha = exit 1 (bloqueante); outros = warning apenas
- Script pesado → **não conectar ao Stop hook automaticamente**

### `.claude/hooks/pre-commit-reminder.sh`
- Verifica `git diff --cached --name-only` inclui `CHANGELOG.md`
- Se não incluído: mostra aviso + aguarda Enter (interativo)
- Ideal para uso como **git hook** manual

```bash
chmod +x .claude/hooks/*.sh
```

---

## Fase 2: Atualizar `.claude/settings.json`

Adicionar apenas a seção `hooks` com eventos válidos do Claude Code.
**Manter** campos existentes: `name`, `description`, `plansDirectory`.
**NÃO adicionar** campos inválidos do doc (skills, quality, context, workflow, conventions).

```json
{
  "$schema": "https://json.schemastore.org/claude-code-settings.json",
  "name": "Grimperium Project Settings",
  "description": "Configurações alinhadas com awesome-claude-code para o projeto Grimperium.",
  "plansDirectory": "./docs/plans",
  "hooks": {
    "SessionStart": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "bash .claude/hooks/pre-plan.sh 2>/dev/null || true"
          }
        ]
      }
    ]
  }
}
```

**Racional:** Apenas `pre-plan.sh` faz sentido em `SessionStart` (leve, informativo). `post-implement.sh` seria muito pesado para rodar no `Stop` (testes demoram). `pre-commit-reminder.sh` funciona melhor como git hook manual.

---

## Fase 3: Nova Skill — `grimperium-mopac.md`

**Arquivo:** `.claude/skills/grimperium-mopac.md`

Conteúdo especializado em parsing de arquivos MOPAC:
- Formatos `.aux`, `.out`, `.mop`
- Notação Fortran científica (`D` → `E`: `+0.577997D+02` = 577.997)
- Extração de descritores: `HEAT_OF_FORMATION`, `DIPOLE:DEBYE`, `IONIZATION_POTENTIAL`, `EIGENVALUES`
- Cálculo HOMO/LUMO via índice `NUM_ELECTRONS`
- Padrões de código com exemplos reais
- Arquivos-chave: `src/grimperium/crest_pm7/mopac_descriptors.py`, `mopac_optimizer.py`

---

## Fase 4: Nova Skill — `grimperium-csv.md`

**Arquivo:** `.claude/skills/grimperium-csv.md`

Conteúdo especializado em gestão do schema CSV:
- Schema atual: 54 colunas (Phase C), categorizado por grupos
- Regras de evolução: como adicionar/remover/renomear colunas
- Operações `write_row` e `update_row` com exemplos
- Histórico de breaking changes (Phase C removeu delta_1/2/3, adicionou 10 descritores MOPAC)
- Arquivo-chave: `src/grimperium/crest_pm7/batch/csv_manager.py`

---

## Fase 5: Agente Científico — `scientific-reviewer.md`

**Arquivo:** `.claude/agents/scientific-reviewer.md`
**Criar diretório:** `.claude/agents/`

Agente focado em revisão de código científico:
- Estabilidade numérica (divisão por zero, comparações float)
- Correção científica (conversões de unidade, indexação CREST)
- Performance (evitar O(n²), preferir numpy vetorizado)
- Checklist com exemplos BAD/GOOD
- Gatilho: `@scientific-review`

---

## Fase 6: Criar `.claudeignore`

**Arquivo:** `.claudeignore` (raiz do projeto)

Exclui do contexto Claude:
- `__pycache__/`, `*.pyc`, `*.egg-info/`, `dist/`, `build/`
- `.venv/`, `.tox/`, `.pytest_cache/`, `htmlcov/`
- `data/raw/**` (exceto `data/raw/sample_10.csv`)
- `tmp/**`, `results/large_batches/**`
- Arquivos MOPAC: `*.aux`, `*.out`, `*.mop`, `*.arc`
- `.DS_Store`, `Thumbs.db`

---

## Arquivos Criados/Modificados

| Ação | Arquivo |
|------|---------|
| **CREATE** | `.claude/hooks/pre-plan.sh` |
| **CREATE** | `.claude/hooks/post-implement.sh` |
| **CREATE** | `.claude/hooks/pre-commit-reminder.sh` |
| **CREATE** | `.claude/skills/grimperium-mopac.md` |
| **CREATE** | `.claude/skills/grimperium-csv.md` |
| **CREATE** | `.claude/agents/scientific-reviewer.md` |
| **CREATE** | `.claudeignore` |
| **MODIFY** | `.claude/settings.json` |

---

## Verificação

```bash
# 1. Testar script pre-plan
bash .claude/hooks/pre-plan.sh

# 2. Testar script post-implement (dentro do venv)
bash .claude/hooks/post-implement.sh

# 3. Validar JSON do settings.json
python3 -c "import json; json.load(open('.claude/settings.json')); print('OK')"

# 4. Verificar scripts executáveis
ls -la .claude/hooks/*.sh

# 5. Verificar que .claudeignore existe
cat .claudeignore
```

---

## Não Implementado (Intencionalmente)

- Campos inválidos no `settings.json` (`skills`, `quality`, `context`, `workflow`, `conventions`) — esses não existem no schema real do Claude Code
- Hook `pre-commit-reminder.sh` como `PreToolUse` automático — seria triggado para **todo** comando Bash, não apenas git commits
- Hook `post-implement.sh` como `Stop` automático — rodar testes a cada Stop seria muito lento e intrusivo
