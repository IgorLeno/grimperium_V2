🔧 PROMPT 001: SCAFFOLDING INICIAL - Grimperium v0.2.0

🎯 Contexto do Projeto

Grimperium v0.2.0 é um framework production-ready de ML ensemble para predição de propriedades termodinâmicas moleculares, com delta-learning para correção de semiempíricos (PM7).

- Dataset Base: Chemperium (~52k moléculas, H298_CBS, H298_B3, S298, Cp_1...45, SMILES, XYZ)
- Objetivo: Treinar KRR + XGBoost em deltas (y_delta = H298_CBS - H298_PM7)
- Features: Tabular (nheavy, charge, multiplicity) + Morgan FP (256 bits) + RDKit descriptors
- Validação: K-fold CV + hold-out test, métrica principal = RMSE vs CBS

---

🐛 Problema / Requisição

Estado atual: Repositório zerado (sem arquivos).

Necessário agora: Criar scaffolding completo e production-ready do Grimperium v0.2.0, incluindo:
- Estrutura de pastas modular (src/grimperium/)
- Configuração Poetry + tox + pre-commit + CI
- Arquivos base com stubs e docstrings completas
- Fixtures in-memory para testes
- README com diagrama arquitetural
- Documentação conceitual (delta-learning, features)

Restrições:
- ✅ Scaffolding apenas (zero implementação de lógica de ML neste batch)
- ✅ ~10-20 arquivos iniciais
- ✅ Estrutura pronta para preencher com lógica depois
- ✅ 100% compatível com Python 3.9-3.12
- ✅ Pronto para CI/CD local e em GitHub Actions

---

💡 Explicação Técnica

Por que Scaffolding Completo Agora?

1. Estrutura Clara = Clareza de Implementação
   - Cada módulo tem responsabilidade bem-definida
   - Novos contributors entendem arquitetura rapidamente

2. CI desde o Início = Qualidade Garantida
   - ruff (linting) + black (formatting) + mypy (types) + pytest (testes)
   - pre-commit hooks = erros no commit, não em produção
   - Coverage tracking (target >85%)

3. Documentação Viva = Conhecimento Compartilhado
   - README com arquitetura visual (ASCII)
   - Docstrings em todas as classes/funções
   - Guias conceituais (delta-learning, features)

4. Poetry + Tox = Packaging Moderno
   - Poetry = gerenciamento de deps + packaging limpo
   - Tox = testa em Python 3.9, 3.10, 3.11, 3.12 automaticamente
   - pyproject.toml = single source of truth

---

✅ Solução Proposta

Uma única solução: Gerar scaffolding completo em uma única batida, com:

1. Estrutura de pastas conforme arquitetura abaixo
2. Arquivos base com docstrings, zero lógica
3. Configurações (Poetry, tox, pre-commit, CI)
4. Fixtures in-memory para testes rápidos
5. Documentação inicial (README + docs/)

---

📝 Instruções para Claude Code

🔄 Planejamento Sugerido (@plan)

1. Batch 1 - Estrutura Base (5 min)
   - Criar árvore de pastas (src/grimperium/, tests/, docs/, .github/)
   - Criar arquivos vazios/stubs com comentários

2. Batch 2 - Configuração & Dependências (5 min)
   - pyproject.toml (Poetry, deps, scripts)
   - tox.ini (multi-Python testing)
   - .pre-commit-config.yaml (hooks)
   - .github/workflows/ci.yml (GitHub Actions)

3. Batch 3 - Módulos Base com Docstrings (10 min)
   - src/grimperium/__init__.py + config.py + api.py
   - src/grimperium/data/*.py (loader, fusion, semiempirical)
   - src/grimperium/models/*.py (base, KRR, XGB, ensemble)
   - src/grimperium/core/*.py (delta_learning, metrics)
   - src/grimperium/utils/*.py (logging, validation, features)
   - Cada arquivo: 20-40 linhas de imports + classe abstract com docstring

4. Batch 4 - Testes Base & Fixtures (5 min)
   - tests/fixtures/mock_data.py (Chemperium mock + PM7 mock)
   - tests/unit/*.py (stubs com pytest.mark.skip ou assertions boilerplate)
   - tests/integration/test_pipeline.py (stub)

5. Batch 5 - Documentação Inicial (5 min)
   - README.md com:
     - Descrição do projeto
     - ASCII diagram da arquitetura
     - Como instalar, rodar testes
     - Links para docs/
   - docs/architecture.md (detalhe da arquitetura)
   - docs/delta_learning_guide.md (conceito do delta)
   - docs/feature_engineering.md (descritores moleculares)
   - CHANGELOG.md (v0.1 initial scaffold)
   - LICENSE (MIT)

---

📂 Arquivos a Criar

Resumo rápido - total ~25 arquivos:

✅ pyproject.toml
✅ tox.ini
✅ .pre-commit-config.yaml
✅ .github/workflows/ci.yml
✅ .gitignore
✅ README.md
✅ CHANGELOG.md
✅ LICENSE
✅ src/grimperium/__init__.py
✅ src/grimperium/config.py
✅ src/grimperium/api.py
✅ src/grimperium/data/__init__.py
✅ src/grimperium/data/loader.py
✅ src/grimperium/data/fusion.py
✅ src/grimperium/data/semiempirical.py
✅ src/grimperium/models/__init__.py
✅ src/grimperium/models/base.py
✅ src/grimperium/models/kernel_ridge.py
✅ src/grimperium/models/xgboost_model.py
✅ src/grimperium/models/delta_ensemble.py
✅ src/grimperium/core/__init__.py
✅ src/grimperium/core/delta_learning.py
✅ src/grimperium/core/metrics.py
✅ src/grimperium/utils/__init__.py
✅ src/grimperium/utils/logging.py
✅ src/grimperium/utils/validation.py
✅ src/grimperium/utils/feature_engineering.py
✅ tests/__init__.py
✅ tests/fixtures/__init__.py
✅ tests/fixtures/mock_data.py
✅ tests/unit/__init__.py
✅ tests/unit/test_loader.py
✅ tests/unit/test_fusion.py
✅ tests/unit/test_semiempirical.py
✅ tests/unit/test_models.py
✅ tests/unit/test_delta_learning.py
✅ tests/integration/__init__.py
✅ tests/integration/test_pipeline.py
✅ docs/architecture.md
✅ docs/delta_learning_guide.md
✅ docs/feature_engineering.md

---

🧪 Critérios de Validação

Testes Manuais:

1. Verificar Estrutura
   ! ls -R src/grimperium/
   ! ls -R tests/
   ! ls -R docs/
   - Todas as pastas e arquivos presentes

2. Verificar Poetry
   ! poetry install
   ! poetry show
   - Instalação sem erros
   - Todas as deps listadas

3. Verificar Tox
   ! tox
   - Passa em Python 3.9, 3.10, 3.11, 3.12

4. Verificar Pre-commit
   ! pre-commit run --all-files
   - Sem erros (ou apenas warnings esperados)

5. Verificar CI
   ! ruff check .
   ! black --check .
   ! mypy src/
   ! pytest --cov=src/grimperium
   - Ruff: sem critical errors
   - Black: sem mudanças necessárias
   - Mypy: sem erros críticos (warnings ok)
   - Pytest: 100% pass (stubs ok, alguns skip ok)
   - Coverage: rastreável

6. Verificar Imports
   ! python -c "from grimperium import api; from grimperium.models import BaseModel"
   - Imports funcionam sem erro

---

📊 Resultado Esperado

DEPOIS (Após execução):

grimperium/
├── pyproject.toml                    ✅ Poetry config
├── tox.ini                           ✅ Tox config
├── .pre-commit-config.yaml           ✅ Git hooks
├── .github/workflows/ci.yml          ✅ GitHub Actions
├── .gitignore                        ✅ Git ignore
├── README.md                         ✅ Overview + ASCII arch
├── CHANGELOG.md                      ✅ v0.1 changelog
├── LICENSE                           ✅ MIT
│
├── src/grimperium/
│   ├── __init__.py                   ✅ Package init
│   ├── config.py                     ✅ Global config stubs
│   ├── api.py                        ✅ High-level API stubs
│   │
│   ├── data/
│   │   ├── __init__.py               ✅ Package init
│   │   ├── loader.py                 ✅ ChemperiumLoader (stub)
│   │   ├── fusion.py                 ✅ DataFusion (stub)
│   │   └── semiempirical.py          ✅ SemiempiricalHandler (stub)
│   │
│   ├── models/
│   │   ├── __init__.py               ✅ Package init
│   │   ├── base.py                   ✅ BaseModel abstract
│   │   ├── kernel_ridge.py           ✅ KernelRidgeRegressor stub
│   │   ├── xgboost_model.py          ✅ XGBoostRegressor stub
│   │   └── delta_ensemble.py         ✅ DeltaLearningEnsemble stub
│   │
│   ├── core/
│   │   ├── __init__.py               ✅ Package init
│   │   ├── delta_learning.py         ✅ Delta utils stub
│   │   └── metrics.py                ✅ Metrics (MSE, MAE, R²) stub
│   │
│   └── utils/
│       ├── __init__.py               ✅ Package init
│       ├── logging.py                ✅ Logging config stub
│       ├── validation.py             ✅ Input validation stub
│       └── feature_engineering.py    ✅ Morgan FP + RDKit stub
│
├── tests/
│   ├── __init__.py                   ✅ Package init
│   │
│   ├── fixtures/
│   │   ├── __init__.py               ✅ Package init
│   │   └── mock_data.py              ✅ Mock fixtures (Chemperium + PM7)
│   │
│   ├── unit/
│   │   ├── __init__.py               ✅ Package init
│   │   ├── test_loader.py            ✅ Stub tests
│   │   ├── test_fusion.py            ✅ Stub tests
│   │   ├── test_semiempirical.py     ✅ Stub tests
│   │   ├── test_models.py            ✅ Stub tests
│   │   └── test_delta_learning.py    ✅ Stub tests
│   │
│   └── integration/
│       ├── __init__.py               ✅ Package init
│       └── test_pipeline.py          ✅ End-to-end stub
│
└── docs/
    ├── architecture.md               ✅ Detailed architecture
    ├── delta_learning_guide.md       ✅ Delta concept guide
    └── feature_engineering.md        ✅ Features guide

Estatísticas Esperadas:
- Total de arquivos: ~45-50
- Linhas de código: ~1500-2000 (stubs + docstrings)
- Arquivo maior: pyproject.toml, tox.ini, README.md (~100 linhas cada)
- Test coverage rastreável: ~0% agora (testes são stubs), mas estrutura 100% pronta

Comandos Funcionais:

✅ poetry install
✅ pytest tests/
✅ ruff check .
✅ black .
✅ mypy src/
✅ tox
✅ pre-commit run --all-files
✅ python -c "from grimperium.models import BaseModel; print('OK')"

---

📌 Notas Importantes

Escopo Explícito - O QUE NÃO FAZER

- ❌ Implementar lógica de ML (fit/predict reais)
- ❌ Carregar dados reais do Chemperium
- ❌ Calcular PM7 via MOPAC
- ❌ Treinar modelos
- ❌ Gerar métricas reais

Escopo Explícito - O QUE FAZER

- ✅ Estrutura completa de pastas
- ✅ Imports, tipos, docstrings
- ✅ Configuração Poetry/tox/pre-commit
- ✅ Fixtures mock em-memory
- ✅ Testes stub (mark.skip ou assert True)
- ✅ README + docs conceituais
- ✅ CI pronto para rodar

Para Próximos Batches

Depois deste scaffold, os próximos prompts implementarão:
1. Batch 2: ChemperiumLoader + DataFusion (lógica real)
2. Batch 3: Models base (KRR, XGB, Ensemble)
3. Batch 4: Feature engineering (Morgan FP + RDKit)
4. Batch 5: Delta-learning orchestration
5. Batch 6: PM7 handler (CREST + MOPAC integration)
6. Batch 7: Testes completos (100% coverage)
7. Batch 8: API e documentação final

---

🎯 Próximo Passo

Após este scaffolding ser gerado e validado:
1. Executar todos os comandos de validação acima
2. Confirmar que estrutura está 100% pronta
3. Chamar próximo prompt para implementação (Batch 2: ChemperiumLoader)

Você está pronto! 🚀
