# Grimperium — Corrigindo Química com Machine Learning

**Python:** 3.10+ | **Licença:** MIT | **Última atualização:** 2026-03-11

---

## O que é este projeto?

Grimperium é um framework científico que usa **machine learning** para tornar cálculos de química quântica mais rápidos e precisos.

### O problema que ele resolve

Calcular propriedades termodinâmicas de moléculas (como a entalpia — a energia armazenada numa molécula) com alta precisão é muito caro computacionalmente. O método mais preciso usado aqui, chamado **CBS-QB3**, pode levar horas por molécula. Já o método mais rápido, chamado **PM7** (semiempírico), termina em segundos — mas comete erros sistemáticos de 10–50 kcal/mol.

### A solução: Delta-Learning

Em vez de substituir o método preciso, o Grimperium ensina um modelo de ML a **corrigir os erros do método rápido**:

```
Predição final = PM7 (rápido, ~80% correto) + Correção ML (aprende os ~20% restantes)
```

Essa abordagem chama-se **delta-learning** e tem duas vantagens principais:
- O modelo de ML tem uma tarefa muito mais simples (aprender erros, não valores absolutos)
- Converge com menos dados de treinamento

### Em números

| Dataset | Moléculas | Uso |
|---|---|---|
| `thermo_cbs_chon.csv` | 29.568 moléculas CHON | Referência de alta precisão (CBS-QB3) |
| `thermo_pm7.csv` | Sendo acumulado | Resultados PM7 gerados pelo pipeline |

O dataset cobre moléculas formadas apenas por **C, H, O e N** — escolha deliberada para manter a física homogênea e facilitar o aprendizado do delta.

---

## Onde estamos agora?

```
[✅ Phase A]  Validação do pipeline CREST + MOPAC (PM7)
[✅ Phase B]  Implementação da CLI interativa
[✅ Phase C]  Error handling, retry system, estabilização do pipeline
[🔄 Agora  ]  Acumulando moléculas PM7 → meta: 1.000 moléculas OK
[⏳ Próximo]  Treinar o primeiro modelo ML com dados reais
[⏳ Futuro ]  Escalar para 3k → 5k → 29k moléculas → deploy
```

**Status dos quality gates:**

| Gate | Status |
|---|---|
| `pytest` | ✅ 494/495 passing |
| `mypy --strict` | ✅ Passando |
| `ruff` | ✅ Passando |
| `black` | ✅ Passando |
| `pre-commit` | ✅ Ativo |
| Cobertura de testes | ~82% (meta: ≥ 85%) |

---

## Como o pipeline funciona?

Para cada molécula (representada como SMILES — uma string de texto):

1. **CREST** gera múltiplas conformações 3D (diferentes "formas" que a molécula pode assumir)
2. **MOPAC PM7** otimiza cada conformação e calcula a entalpia
3. O conformer mais estável (menor entalpia) é selecionado
4. Os resultados são salvos no CSV com descritores eletrônicos (HOMO, LUMO, GAP)
5. *(Futuro)* O modelo ML aplica a correção delta → predição final

---

## Instalação

### Pré-requisitos

Você vai precisar de:
- **Python 3.10 ou superior** — [download aqui](https://www.python.org/downloads/)
- **Git** — [download aqui](https://git-scm.com/)
- **CREST** (gerador de conformações) — [instruções oficiais](https://crest-lab.github.io/crest-docs/)
- **MOPAC** (otimizador PM7) — [download gratuito](https://openmopac.net/)

> **Nota:** CREST e MOPAC são programas externos de química computacional. Para usar apenas o módulo de ML (treino e predição), eles não são obrigatórios.

### Passo a passo

**1. Clone o repositório**
```bash
git clone https://github.com/IgorLeno/grimperium_V2.git
cd grimperium_V2
```

**2. Crie um ambiente virtual**
> Um ambiente virtual isola as dependências do projeto para não conflitar com outros projetos Python na sua máquina.

```bash
# Criar o ambiente
python -m venv venv

# Ativar no Linux/macOS
source venv/bin/activate

# Ativar no Windows
venv\Scripts\activate
```

Você saberá que funcionou quando `(venv)` aparecer no início do terminal.

**3. Instale o projeto e suas dependências**
```bash
# Instala o Grimperium em modo de desenvolvimento
pip install -e .

# Instala as dependências de desenvolvimento (testes, linting)
pip install -e ".[dev]"
```

**4. Configure os pre-commit hooks** *(opcional, recomendado para contribuidores)*
```bash
pre-commit install
```

**5. Valide a instalação**
```bash
# Deve retornar "OK"
python -c "from grimperium import GrimperiumAPI; print('OK')"

# Rode os testes
pytest tests/ -v
```

---

## Uso básico

### Iniciar a CLI interativa
```bash
python -m grimperium
```

A CLI oferece menus para:
- Rodar batches do pipeline CREST + PM7
- Monitorar progresso das moléculas
- Visualizar resultados

### Carregar os dados via código
```python
from grimperium.data.loader import ChemperiumLoader

loader = ChemperiumLoader()

# Carregar o dataset de referência (CBS-QB3)
df_cbs = loader.load_thermo_cbs_chon(max_nheavy=50)
print(f"Moléculas de referência: {len(df_cbs)}")

# Carregar os resultados PM7 acumulados
df_pm7 = loader.load_thermo_pm7()
print(f"Moléculas PM7 processadas: {len(df_pm7)}")
```

### Treinar o modelo delta (quando houver dados suficientes)
```python
from grimperium.core.delta_learning import DeltaLearner

learner = DeltaLearner()
learner.load_data("data/thermo_cbs_chon.csv", "data/thermo_pm7.csv")
learner.compute_features()
learner.train()

metrics = learner.evaluate()
print(f"RMSE delta: {metrics['rmse']:.2f} kcal/mol")
print(f"R²: {metrics['r2']:.4f}")
```

---

## Estrutura do projeto

```text
grimperium_V2/
├── src/grimperium/
│   ├── core/                    # DeltaLearner, métricas, orquestração
│   ├── models/                  # KRR, XGBoost, DeltaLearningEnsemble
│   ├── data/                    # Carregamento e fusão de datasets
│   ├── crest_pm7/               # Pipeline CREST + MOPAC
│   │   ├── batch/               # Gerenciamento de batches e retry system
│   │   ├── conformer_generator.py  # Invoca o CREST
│   │   ├── mopac_optimizer.py   # Invoca o MOPAC
│   │   └── mopac_descriptors.py # Parser dos arquivos .out do MOPAC
│   ├── cli/                     # Interface de linha de comando (Rich)
│   └── utils/                   # Logging, validação, helpers
├── tests/                       # Testes unitários e de integração
├── data/
│   ├── thermo_cbs_chon.csv      # 29.568 moléculas CHON — referência CBS-QB3
│   └── thermo_pm7.csv           # Resultados PM7 acumulados pelo pipeline
└── docs/                        # Documentação técnica detalhada
```

---

## Modelos de ML

O ensemble usa dois modelos complementares:

| Modelo | Papel | Por quê? |
|---|---|---|
| **KRR** (Kernel Ridge Regression) | Aprende correções suaves e contínuas | Kernel RBF captura similaridade molecular |
| **XGBoost** | Captura interações não-lineares complexas | Excelente em dados tabulares |
| **Ensemble** | Média ponderada (50/50 padrão) | Reduz variância, melhora robustez |

### Features moleculares (270 dimensões)

- **Tabular (3):** número de átomos pesados, carga, multiplicidade
- **Morgan Fingerprints (256):** padrões de subestrutura circular (ECFP-like, raio 3)
- **RDKit (11):** peso molecular, TPSA, LogP, ligações rotacionáveis, etc.

---

## Testes e qualidade de código

```bash
# Rodar todos os testes
pytest tests/ -v

# Com relatório de cobertura
pytest tests/ --cov=src/grimperium --cov-report=html
# Abrir o relatório
open htmlcov/index.html   # macOS
xdg-open htmlcov/index.html  # Linux

# Verificar tipos (mypy)
mypy src/ --strict

# Linting
ruff check src/

# Formatação
black src/ tests/

# Rodar todos os gates de uma vez
pre-commit run --all-files
```

---

## Roadmap

- **Agora:** Acumular ≥ 1.000 moléculas PM7 com `quality_grade` A ou B
- **Próximo:** Treinar primeiro modelo — meta RMSE ≤ 15 kcal/mol, R² > 0.9
- **Depois:** Escalar para 3k → 5k → 29k moléculas
- **Futuro:** Hiperparametrização (grid/Bayesian), validação cruzada k-fold, deploy via `api.py`

---

## Tecnologias

- **Python** 3.10+, **pandas**, **numpy**, **scikit-learn**, **XGBoost**
- **RDKit** — cheminformatics e geração de features
- **Rich** + **Questionary** — CLI interativa
- **CREST** v3 + **MOPAC** (OpenMOPAC) — cálculos externos de química computacional
- **pytest**, **mypy**, **ruff**, **black**, **pre-commit** — qualidade de código

---

## Documentação técnica

| Documento | Conteúdo |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Arquitetura detalhada do sistema |
| [`docs/DATASETS.md`](docs/DATASETS.md) | Guia completo dos datasets |
| [`docs/CREST_INTEGRATION.md`](docs/CREST_INTEGRATION.md) | Integração e parâmetros do CREST |
| [`docs/MOPAC_INTEGRATION.md`](docs/MOPAC_INTEGRATION.md) | Integração e descritores do MOPAC |

---

## Autor

Igor Leno — São Paulo, BR

---

## Licença

MIT License — veja [LICENSE](LICENSE) para detalhes.
