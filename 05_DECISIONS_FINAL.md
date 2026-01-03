📊 DECISÕES FINAIS CONSOLIDADAS - Grimperium v0.2.0

✅ Decisões Confirmadas

1. Semiempírico: PM7
- Escolha: PM7 (Stewart, 2013) via MOPAC
- Justificativa: PM7 é atualmente o melhor semiempírico de propósito geral, especialmente para:
  - ✅ Cálculo de entalpia de formação (propriedade alvo)
  - ✅ Geometrias de equilíbrio (importante para XYZ-dependent features)
  - ✅ Multiplicidade de spin (radicais bem descritos)
  - ✅ Suporta uma gama ampla de elementos (C, H, O, N, S, halogênios)
  
- Alternativas consideradas:
  - ⚠️ PM6-D3H+: mais novo, mas menos consolidado para entalpia; bom para complexos metal
  - ⚠️ OM3: específico para espectroscopia UV-Vis, não otimizado para entalpia
  
- Conclusão: PM7 é a escolha mais sólida e estabelecida para esta aplicação

---

2. Pipeline de Cálculo PM7
- Workflow: CREST → MOPAC/PM7
- Etapas:
  1. CREST: conformational search (xTB rápido) para cada SMILES
  2. MOPAC/PM7: otimização de geometria + cálculo de entalpia de formação
  3. Extração de H298_PM7 para cada molécula
  4. Merge com Chemperium dataset

- Desafios:
  - 52.837 moléculas × CREST + MOPAC = computacionalmente custoso (planejado em batches)
  - Possível integração futura com HPC ou Colab GPU
  - Para v0.1: Implementar stubs, deixar lógica de orchestração pronta

---

3. Estratégia Delta: OPÇÃO A (Delta Simples)
- Definição: delta_PM7 = H298_CBS - H298_PM7
- Treinamento: Modelo ML aprender y_delta diretamente
- Predição: H298_CBS ≈ H298_PM7 + delta_model.predict(X)
- Vantagens:
  - ✅ Simples, direto, interpretável
  - ✅ Fácil de integrar em pipeline
  - ✅ Fundação para evoluir para Opção B (ensemble multi-delta) depois
  - ✅ Baseline forte: semiempírico já carrega 70-80% da informação, delta completa <20%

---

4. Features: OPÇÃO D (Híbrida)
- Componentes:
  1. Tabular Básico:
     - nheavy (número de átomos pesados)
     - charge (carga total)
     - multiplicity (multiplicidade de spin)
  
  2. Fingerprints Moleculares:
     - Morgan Fingerprints (RDKit)
     - Tamanho: 256 ou 512 bits (balanço entre performance e interpretabilidade)
  
  3. Descritores RDKit:
     - Descriptors.MolWt() (peso molecular)
     - Descriptors.TPSA() (polar surface area)
     - Descriptors.LogP() (lipofilia)
     - Opcionalmente: NumRotatableBonds, NumHBD, NumHBA, etc.

- Vantagens:
  - ✅ Não requer XYZ (evita overhead geométrico)
  - ✅ Rápido de computar via RDKit
  - ✅ Combina informação tabular + estrutural (SMILES)
  - ✅ Dimensionalidade moderada (~280-550 features total)
  - ✅ Interpretável (features conhecidas, não "black box")

---

5. Métrica de Validação: OPÇÃO A (vs CBS)
- Métrica Principal: 
  - MSE(H298_PM7 + delta_ML vs H298_CBS) ou RMSE
  - MAE(H298_PM7 + delta_ML vs H298_CBS)
  - R² (Coefficient of Determination)

- Meta de Performance:
  - RMSE(delta_corrected) < RMSE(H298_PM7 puro)
  - Idealmente: RMSE(delta_corrected) ≈ 0.5-2 kcal/mol (dependendo da precisão CBS)

- Validação:
  - K-fold Cross-Validation (k=5 ou 10) durante desenvolvimento
  - Hold-out Test Set (20% dados) para avaliação final
  - Análise de resíduos: distribuição de erros, outliers, viés

- Comparativo (Bônus):
  - Mostrar MAE(H298_PM7 vs CBS) vs MAE(SEMIEMP+delta vs CBS)
  - Mostrar que delta_PM7 é comparável ou superior a delta_B3LYP
  - Comunicar ganho percentual

---

📐 Arquitetura Modificada para PM7 + Delta

src/grimperium/
├── data/
│   ├── __init__.py
│   ├── loader.py              # ChemperiumLoader (CBS, B3, Cp, etc.)
│   ├── semiempirical.py       # ✨ NOVO: SemiempiricalHandler (PM7 calc)
│   └── fusion.py              # DataFusion (merge + delta computation)
├── models/
│   ├── __init__.py
│   ├── base.py                # BaseModel (abstract)
│   ├── kernel_ridge.py        # KernelRidgeRegressor
│   ├── xgboost_model.py       # XGBoostRegressor
│   └── delta_ensemble.py      # DeltaLearningEnsemble
├── core/
│   ├── __init__.py
│   ├── delta_learning.py      # Delta-learning logic & utils
│   └── metrics.py             # Métricas (MSE, MAE, R², etc.)
├── utils/
│   ├── __init__.py
│   ├── logging.py             # Logging configurável
│   ├── validation.py          # Validação de entrada
│   └── feature_engineering.py # ✨ NOVO: Morgan FP + RDKit descriptors
├── config.py                  # Configuração global
└── api.py                     # Orquestração high-level

---

🎯 Timeline Revisado (v0.1 → v0.2)

v0.1 (Agora: Dec 2024)
- ✅ Scaffolding completo (arquitetura, CI, docs base)
- ✅ ChemperiumLoader + DataFusion (stubs + testes)
- ✅ BaseModel + KRR + XGB (stubs)
- ✅ Delta-learning core (conceito + interfaces)
- ✅ Feature engineering (tabular + Morgan FP + RDKit)
- ⏳ PM7 Handler (stub, orchestração design)
- ✅ Fixtures in-memory (mock data)

v0.2 (Próximo: Jan-Feb 2025)
- ✅ Implementar PM7 calculation pipeline (CREST + MOPAC)
- ✅ Integrar dados PM7 reais ao loader
- ✅ Treinar KRR + XGB nos deltas reais
- ✅ Validar métricas (RMSE, MAE, R²)
- ✅ Comparar vs B3LYP delta
- ✅ Deploy em Colab (integrations/colab.py)
- ✅ Publicação em PyPI

---

💾 Novas Dependências

[tool.poetry.dependencies]
# ML & Data
scikit-learn = "^1.3.0"
xgboost = "^2.0.0"
pandas = "^2.0.0"
numpy = "^1.24.0"

# Química & Features
rdkit = "^2023.9.0"           # Fingerprints, descritores, SMILES parsing

# Semiempírico (Opcional, vem depois)
# mopac = "^22.1"             # MOPAC/PM7 (via API ou subprocess)
# crest = "^2.13"             # xTB-based conformer search

# Dev & QA
pytest = "^7.4.0"
pytest-cov = "^4.1.0"
pytest-xdist = "^3.3.0"
mypy = "^1.5.0"
ruff = "^0.1.0"
black = "^23.9.0"
pdoc = "^14.0.0"

# Pre-commit & automation
pre-commit = "^3.4.0"

---

🔬 Notas Técnicas Importantes

PM7 vs Alternativas (para Entalpia de Formação)

| Semiempírico | Entalpia Formação | Geometria | Radicals | Speed | Notes |
|---|---|---|---|---|---|
| **PM7** ✅ | Excelente | Muito boa | Bom | Rápido | **ESCOLHIDO** - Melhor overall |
| PM6-D3H+ | Bom | Excelente | Muito bom | Lento | Melhor para complexos |
| OM3 | Bom | Boa | Razoável | Moderado | Otimizado para UV-Vis |
| AM1-SCC | Aceitável | Razoável | Razoável | Rápido | Desatualizado |

---

✨ Próximo: GRANDE PROMPT para Claude Code

Com essas decisões consolidadas, vou gerar um prompt único e bem-estruturado que:

1. ✅ Cria scaffolding completo
2. ✅ Inclui stubs para componentes
3. ✅ Fixtures in-memory com deltas mockados
4. ✅ Tests básicos e estrutura
5. ✅ Documentação
6. ✅ CI integrado

Próximo passo: Você valida esse documento e eu gero o prompt final!
