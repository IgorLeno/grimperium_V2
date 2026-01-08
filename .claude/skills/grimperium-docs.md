---
name: grimperium-docs
description: Gera/atualiza documentação Sphinx, READMEs, changelog e relatório técnico
tools: [bash, file]
context: fork
user-invocable: true
allowed-tools:
  - Bash(poetry *)
  - Bash(sphinx *)
  - Bash(git *)
  - Bash(rm *)
  - Bash(find *)
---

# 📚 Skill: Generate Grimperium Documentation

**Propósito:** Automatizar geração de documentação completa do projeto

**Quando usar:**
- Antes de commitar features novas
- Antes de criar releases
- Quando quiser gerar relatório técnico
- Sempre que módulos mudarem

## O que esta skill faz

1. ✅ **Sphinx Documentation** — Gera docs HTML
2. ✅ **Module READMEs** — README.md por módulo
3. ✅ **CHANGELOG.md** — Atualiza com commits recentes
4. ✅ **Relatório Técnico** — Consolida métricas do projeto
5. ✅ **GitHub Pages** — Deploy automático (opcional)

## Como Usar

### Modo Completo (Recomendado)

```
@claude /grimperium-docs
```

Executa tudo:
```
1. poetry install -E docs
2. sphinx-build -b html docs/ docs/html/
3. Gera READMEs por módulo
4. Atualiza CHANGELOG.md
5. Gera relatório técnico
6. Opcional: Deploy para GitHub Pages
```

### Modos Específicos

```
@claude /grimperium-docs --sphinx-only      # Só Sphinx
@claude /grimperium-docs --module-readmes   # Só READMEs por módulo
@claude /grimperium-docs --changelog        # Só changelog
@claude /grimperium-docs --technical-report # Só relatório técnico
@claude /grimperium-docs --github-pages     # Deploy docs
```

## Workflow Recomendado

```
Desenvolvimento → Feature completa → @claude /grimperium-docs
                                                  ↓
                    Documentação gerada automaticamente
                                                  ↓
                    git add docs/ READMEs CHANGELOG.md
                                                  ↓
                    git commit -m "docs: update documentation"
                    git push
```

## Estrutura Gerada

```
docs/html/                    # Sphinx HTML docs
├── index.html               # Página inicial
├── modules.html             # Lista de módulos
├── grimperium.core.html     # Documentação core
├── grimperium.data.html     # Documentação data
└── grimperium.models.html   # Documentação models

READMEs/                      # READMEs por módulo
├── README_core.md
├── README_data.md
├── README_models.md
└── README.md                 # Principal (consolidado)

CHANGELOG.md                  # Atualizado
TECHNICAL_REPORT.md           # Relatório consolidado
```

## Processos Automatizados

### 1. Sphinx Documentation

```
poetry install -E docs
sphinx-apidoc -o docs/source src/grimperium
sphinx-build -b html docs/ docs/html/
```

**Resultado:** `docs/html/index.html` pronto para GitHub Pages

### 2. Module READMEs

Para cada módulo (`core`, `data`, `models`):

```
# README_core.md
## Grimperium Core

### DeltaLearner
```python
from grimperium.core.delta_learning import DeltaLearner
```

### Metrics
```python
from grimperium.core.metrics import compute_rmse
```

**Coverage:** XX%
**Status:** ✅ Production Ready
```

### 3. CHANGELOG.md Update

```
## [Unreleased]

### Added
- [ ] Sua nova feature aqui

### Changed
- [ ] Mudanças aqui

### Fixed
- [ ] Fixes aqui
```

Auto-populado com commits desde último release.

### 4. Technical Report

```
# Grimperium Technical Report

## Project Status

| Metric | Value |
|--------|-------|
| Total Modules | 14 |
| Test Coverage | 95% |
| Test Status | 88 passed, 11 xfailed |
| Python Versions | 3.10-3.12 |
| Dependencies | 12 production, 6 dev |

## Architecture Overview

```
grimperium/
├── core/          # DeltaLearning logic
├── data/          # Loaders + fusion
├── models/        # ML models
└── utils/         # Helpers
```

## Next Steps

- [ ] Batch 4: Model configuration system
- [ ] Batch 5: Hyperparameter optimization
- [ ] Batch 6: Scale to 52k molecules
```

### 5. GitHub Pages (Opcional)

```
git subtree push --prefix docs/html origin gh-pages
```

Docs ficam em: `https://github.com/IgorLeno/grimperium_V2/docs`

## Output Detalhado

```
1️⃣ Sphinx Documentation
   ├─ Generated docs/html/index.html ✅
   ├─ 14 modules documented
   └─ Ready for GitHub Pages

2️⃣ Module READMEs
   ├─ README_core.md created
   ├─ README_data.md created
   ├─ README_models.md created
   └─ README.md (main) updated

3️⃣ CHANGELOG.md
   ├─ [Unreleased] section added
   └─ Auto-populated from git log

4️⃣ Technical Report
   ├─ TECHNICAL_REPORT.md generated
   └─ Project metrics included

5️⃣ GitHub Pages (opcional)
   └─ Deployed to gh-pages ✅

📊 SUMMARY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Files generated: 18
Modules documented: 14
Coverage: 95%
Status: ✅ COMPLETE ✅
```

## Pré-requisitos

```
# Instalar dependências de docs
poetry add --group docs sphinx sphinx-rtd-theme sphinx-autodoc-typehints

# Configurar Sphinx
sphinx-quickstart docs/
```

**Ação esperada:** Instalar dependências se necessário.

## Integração com Git Hooks

Para automatizar completamente:

```bash
# .git/hooks/pre-commit
@claude /grimperium-docs --sphinx-only

# .git/hooks/post-merge
@claude /grimperium-docs
```

## Performance

```
Sphinx: 30-60s (primeira vez)
Module READMEs: 5s
CHANGELOG: 2s
Technical Report: 3s
Total: ~1-2 minutos (background)
```

## Notas Importantes

- Roda em **background** (`context: fork`) — não bloqueia desenvolvimento
- Se Sphinx não estiver configurado, pergunta se quer configurar
- GitHub Pages requer `gh-pages` branch
- Technical Report é gerado a partir de análise do código atual

## Se Algo Falhar

```
❌ Sphinx não configurado → "Quer configurar agora?"
❌ Dependências faltando → "Instalar poetry install -E docs?"
❌ GitHub Pages sem gh-pages → "Criar branch gh-pages?"
```

Skill é resiliente e guia o usuário.

## Próximas Skills Sugeridas

1. `/grimperium-deploy` — Deploy automático
2. `/grimperium-data-analyze` — Análise dataset 52k
3. `/grimperium-benchmark` — Benchmarks vs baselines
