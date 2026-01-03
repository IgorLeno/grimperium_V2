🚀 START HERE - Grimperium v0.2.0 Setup

⚡ TL;DR (30 SEGUNDOS)

Você tem tudo pronto para começar o Grimperium. Siga 3 passos:

1️⃣ Copie o Prompt
Arquivo: PROMPT_001_SCAFFOLDING_INICIAL.md (o segundo arquivo para download)
Ação: Selecionar tudo (Ctrl+A) → Copiar (Ctrl+C)

2️⃣ Cole no Cursor
Cursor: Pressione Cmd+K (Mac) ou Ctrl+K (Windows)
Digite: @plan
Cole o prompt e pressione Enter
Tempo: ~25 minutos

3️⃣ Valide o Resultado
```bash
poetry install
pytest tests/ -v
ruff check .
# Se passar: ✅ Sucesso!
```

---

📋 Quick Navigation

| Preciso... | Arquivo | Tempo |
|-----------|---------|-------|
| **Começar agora** | PROMPT_001_SCAFFOLDING_INICIAL.md | 25 min |
| Entender o que vai acontecer | RESUMO_EXECUTIVO_SCAFFOLDING.md | 5 min |
| Passo-a-passo durante execução | INSTRUCOES_CLAUDE_CODE.md | usar conforme precisa |
| Ver todas as decisões | decisions_final_consolidated.md | 10 min |
| Entender o dataset e delta | dataset_context_and_delta_strategy.md | 10 min |

---

🎯 O Que Será Criado

grimperium/
├── src/grimperium/          ← Core ML ensemble framework
├── tests/                   ← Full test suite
├── docs/                    ← Architecture + guides
├── pyproject.toml           ← Poetry config
├── tox.ini                  ← Multi-Python testing
├── .github/workflows/ci.yml ← GitHub Actions CI
└── README.md                ← High-level overview

Total: ~25 arquivos, ~1500 linhas, pronto para implementação

---

✅ Pré-Requisitos (Checklist)

Antes de começar, confirme:

- [ ] Cursor está aberto
- [ ] Repositório grimperium está zerado (vazio ou só .git)
- [ ] Terminal aberto na raiz do projeto
- [ ] Você tem Python 3.9+ instalado
- [ ] Você pode executar: poetry, python, git

Se tudo OK → Vamos começar!
