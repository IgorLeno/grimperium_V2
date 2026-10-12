# Grimperium — Aproximando velocidade e precisão em cálculos termodinâmicos

[![CI](https://github.com/IgorLeno/grimperium_V2/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/IgorLeno/grimperium_V2/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue?logo=python&logoColor=white)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![codecov](https://codecov.io/gh/IgorLeno/grimperium_V2/branch/main/graph/badge.svg)](https://codecov.io/gh/IgorLeno/grimperium_V2)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Checked with mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://pre-commit.com/)
[![Claude Code](https://img.shields.io/badge/Claude%20Code-assisted-blueviolet?logo=anthropic)](https://claude.ai/code)
[![CodeRabbit](https://img.shields.io/coderabbit/prs/github/IgorLeno/grimperium_V2?utm_source=oss&utm_medium=github&utm_campaign=IgorLeno%2Fgrimperium_V2&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)](https://coderabbit.ai)

---

## Por que o Grimperium existe?

Calcular propriedades termodinâmicas de moléculas envolve uma escolha difícil: **métodos computacionais rápidos costumam perder precisão; métodos mais precisos custam muito mais tempo de processamento**. O Grimperium investiga como reduzir essa diferença, sem afirmar que um modelo rápido passará a substituir qualquer cálculo de referência.

A ideia é **aprender com os resultados de um método mais preciso para corrigir as estimativas de outro, mais rápido**. Com dados que relacionam os dois métodos, treinamos um modelo de *machine learning* para prever a diferença entre eles. Na aplicação, o método rápido calcula o valor inicial e o modelo treinado fornece uma correção. Essa abordagem é chamada de **delta-learning**.

```text
Para construir a correção:
  molécula → cálculo rápido ─┐
                             ├─ diferença entre resultados → treinar ML
  molécula → referência precisa ┘

Para usar o modelo treinado:
  nova molécula → cálculo rápido → correção estimada pelo ML → estimativa corrigida
```

O trabalho exige muitos pares de resultados para treinamento e validação. Por isso, o Grimperium também desenvolve **pipelines em Python para organizar moléculas, integrar ferramentas de química computacional, automatizar cálculos e preparar dados comparáveis**. A correção não é uma constante universal: ela depende dos métodos, da propriedade e do conjunto de dados que sustentam o treinamento.

### Visão do projeto: um hub de métodos corrigidos

O objetivo de longo prazo é construir um **hub de propriedades e métodos computacionais**, no qual diferentes combinações de método rápido, referência de maior precisão e modelo de correção possam ser configuradas, treinadas, comparadas e utilizadas. O princípio é o mesmo, mas **cada propriedade e par de métodos exigirá seus próprios dados, compatibilidades e validações**.

**Escopo implementado hoje:** o catálogo de propriedades está concentrado na **entalpia padrão de formação (ΔHf°)**. A combinação de referência CBS-QB3 e cálculo rápido PM7 é o exemplo central do fluxo de delta-learning. O repositório também possui rotinas de cálculos semiempíricos e interfaces de acompanhamento, mas isso **não significa que outras propriedades termodinâmicas já estejam suportadas ou corrigidas por ML**.

A explicação conceitual e a separação entre implementação e visão futura estão em [Visão e escopo](docs/PRODUCT_VISION.md). Detalhes do treinamento estão em [Guia de delta-learning](docs/delta_learning_guide.md).

### Validação e acompanhamento

- Um modelo de correção deve ser comparado ao **método rápido sem correção** e à **referência de maior precisão**, em dados de avaliação apropriados. A redução do erro precisa ser medida; não deve ser presumida.
- Gráficos e dashboards ajudam a visualizar diferenças entre métodos, erros, variáveis relevantes e progresso das execuções. Ver [dashboard local](apps/dashboard/README.md).
- As métricas e contagens registradas em relatórios e marcos históricos não devem ser confundidas com o **estado atual da base** ou com uma garantia geral de precisão.

### Marcos históricos documentados

A documentação anterior registrou um primeiro modelo treinado e um *quality gate* de validação (MAE 2,50 kcal/mol, RMSE 3,54 kcal/mol e R² 0,997 na avaliação descrita naquele marco, de 2026-03-30). Essas métricas pertencem àquele experimento e **não caracterizam automaticamente todos os modelos, moléculas ou protocolos futuros**.

As contagens de moléculas exibidas em versões anteriores do README correspondiam a momentos distintos do processamento; para números correntes, deve-se consultar o dataset/registro de execução e o dashboard, em vez de somar ou reaproveitar snapshots antigos.

### Pipeline atual em termos simples

No percurso de entalpia de formação, as ferramentas do projeto permitem preparar estruturas moleculares, executar cálculos com métodos semiempíricos (incluindo PM7), registrar resultados e construir pares comparáveis à referência de maior precisão. Os dados são utilizados para treinar e avaliar a correção por ML. A correção é aplicada **quando o fluxo e o modelo compatível estão configurados**; nem toda execução de cálculo semiempírico realiza delta-learning automaticamente.

O **Semi-Imperium** é um fluxo separado de cálculos e catalogação de ΔHf° com AM1/PM3/PM7; ele **não aplica correção de ML** ([escopo do Semi-Imperium](docs/SEMI_IMPERIUM.md)).

---

## Instalação

> **Importante:** O Grimperium depende de **CREST** e **MOPAC**, ferramentas de química computacional que rodam nativamente apenas em **Linux**. Se você usa Windows, veja a seção [Windows (via WSL 2)](#-windows-via-wsl-2) logo abaixo.

---

### 🐧 Linux (instalação nativa)

Esta é a forma recomendada para desenvolvimento e produção.

#### Pré-requisitos

- **Python 3.10+** — `sudo apt install python3 python3-pip python3-venv`
- **Git** — `sudo apt install git`
- **Poetry** — gerenciador de dependências do projeto
- **CREST** — [instruções oficiais](https://crest-lab.github.io/crest-docs/)
- **MOPAC** — [download gratuito](https://openmopac.net/)

> Para apenas explorar o código ou rodar os testes, CREST e MOPAC não são obrigatórios.

#### Passo a passo

**1. Instale o Poetry**
```bash
curl -sSL https://install.python-poetry.org | python3 -

# Adicione ao PATH (se necessário)
export PATH="$HOME/.local/bin:$PATH"

# Verifique a instalação
poetry --version
```

**2. Clone o repositório**
```bash
git clone https://github.com/IgorLeno/grimperium_V2.git
cd grimperium_V2
```

**3. Instale as dependências com Poetry**
```bash
# Cria o ambiente virtual e instala tudo automaticamente
poetry install --with dev --all-extras

# Para instalar também as dependências de desenvolvimento (testes, linting)
poetry install --with dev --all-extras
```

**4. Ative o ambiente virtual**
```bash
# Poetry 2.x: ative o ambiente virtual com:
eval $(poetry env activate)    # Linux / macOS / WSL
# No Windows PowerShell nativo (fora do WSL), use:
# Invoke-Expression (poetry env activate)

# Ou simplesmente prefixe cada comando com `poetry run` (recomendado):
# poetry run grimperium
# poetry run pytest tests/ -v
```

Se ativar o ambiente, você saberá que funcionou quando `(grimperium-...)` aparecer no início do terminal.

**5. Configure os pre-commit hooks** *(recomendado para contribuidores)*
```bash
poetry run pre-commit install
```

**6. Valide a instalação**
```bash
# Deve retornar "OK"
poetry run python -c "from grimperium import GrimperiumAPI; print('OK')"

# Rode os testes
poetry run pytest tests/ -v
```

---

### 🪟 Windows (via WSL 2)

No Windows, o caminho recomendado é usar o **WSL 2** (Windows Subsystem for Linux) — ele cria um ambiente Linux real dentro do Windows, sem precisar de máquina virtual separada.

> **Por que WSL 2 e não instalar direto no Windows?**
> CREST e MOPAC são compilados para Linux e não possuem versões nativas para Windows. O WSL 2 resolve isso de forma transparente.

#### Parte 1 — Instalar o Git for Windows *(para clonar pelo terminal Windows, opcional)*

1. Acesse **[git-scm.com/download/win](https://git-scm.com/download/win)** e baixe o instalador
2. Execute o instalador deixando **todas as opções padrão**
3. Feche e reabra o terminal após a instalação
4. Verifique: `git --version`

#### Parte 2 — Instalar o WSL 2 com Ubuntu

1. Abra o **PowerShell como Administrador** (`Win + X` → *Terminal (Admin)*)
2. Execute:
```powershell
wsl --install
```
3. **Reinicie o computador** quando solicitado
4. Após reiniciar, o Ubuntu abrirá automaticamente pedindo para criar um **usuário e senha Linux** — defina-os e anote
5. Verifique que está no WSL 2:
```powershell
wsl --list --verbose
# Deve mostrar VERSION 2 ao lado do Ubuntu
```

#### Parte 3 — Configurar o ambiente dentro do WSL

Abra o terminal Ubuntu (pelo Menu Iniciar ou digitando `ubuntu` no PowerShell) e execute:

```bash
# Atualizar pacotes
sudo apt update && sudo apt upgrade -y

# Instalar dependências base
sudo apt install -y python3 python3-pip python3-venv git curl

# Instalar Poetry
curl -sSL https://install.python-poetry.org | python3 -
export PATH="$HOME/.local/bin:$PATH"
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc

# Verificar
poetry --version
```

#### Parte 4 — Clonar e instalar o Grimperium

Dentro do terminal Ubuntu (WSL):

```bash
# Clone o repositório
git clone https://github.com/IgorLeno/grimperium_V2.git
cd grimperium_V2

# Instale as dependências
poetry install --with dev --all-extras

# Ative o ambiente
# Poetry 2.x: ative o ambiente virtual com:
eval $(poetry env activate)    # Linux / macOS / WSL
# No Windows PowerShell nativo (fora do WSL), use:
# Invoke-Expression (poetry env activate)

# Ou simplesmente prefixe cada comando com `poetry run` (recomendado):
# poetry run grimperium
# poetry run pytest tests/ -v

# Valide
poetry run python -c "from grimperium import GrimperiumAPI; print('OK')"
poetry run pytest tests/ -v
```

#### Parte 5 — Acessar os arquivos pelo Explorer do Windows *(opcional)*

Você pode navegar pelos arquivos do WSL diretamente no Windows Explorer digitando na barra de endereço:
```
\\wsl$\Ubuntu\home\<seu-usuario>\grimperium_V2
```

Ou, de dentro do terminal WSL:
```bash
# Abre a pasta atual no Explorer do Windows
explorer.exe .
```

#### Resumo — Linux vs Windows

| | Linux nativo | Windows (WSL 2) |
|---|---|---|
| CREST + MOPAC | ✅ Funciona direto | ✅ Funciona via WSL |
| Performance | ✅ Máxima | ⚡ Muito boa (WSL 2 é quase nativo) |
| Configuração | Simples | Requer instalação do WSL |
| Recomendado para | Produção / servidores | Desenvolvimento no Windows |

---

## Uso básico

### Iniciar a CLI interativa
```bash
poetry run grimperium
```

A CLI oferece menus para:
- Rodar batches do pipeline CREST + PM7
- Monitorar progresso das moléculas
- Visualizar resultados

### Iniciar o Semi-Imperium

O Semi-Imperium é um shell independente e focado para os fluxos operacionais
de cálculo e catálogo de entalpia padrão de formação (`ΔHf°`) com MOPAC
AM1/PM3/PM7. Cada Hamiltoniano é tratado independentemente; CREST é uma busca
conformacional opcional, não a fonte da `ΔHf°` final:

```bash
poetry run semi-imperium
# ou, com o ambiente Poetry ativo:
python -m semi_imperium
```

Sua interface principal expõe somente **CALCULATE**, **DATABASE** e
**SETTINGS**. O comando `grimperium` continua disponível com a interface
completa do projeto.

O guia [Semi-Imperium](docs/SEMI_IMPERIUM.md) descreve o escopo científico,
geração versus seleção de conformeros, Energy Top-N, CONFPASS experimental,
verificação de mínimo local, limites sobre mínimo global, reutilização,
proveniência e armazenamento isolado.

> **Sobre unidades:** a CLI legada do Grimperium exibe H298 em kcal/mol e
> kJ/mol simultaneamente (1 kcal = 4,184 kJ). O Semi-Imperium exibe e persiste
> sua `ΔHf°` em kcal/mol no store JSON próprio. Tempos internos usam segundos.

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
poetry run pytest tests/ -v

# Com relatório de cobertura
poetry run pytest tests/ --cov=src/grimperium --cov-report=html
# Abrir o relatório
open htmlcov/index.html   # macOS
xdg-open htmlcov/index.html  # Linux

# Verificar tipos (mypy)
poetry run mypy src/ --strict

# Linting
poetry run ruff check src/

# Formatação
poetry run black src/ tests/

# Rodar todos os gates de uma vez
poetry run pre-commit run --all-files
```

---

## Roadmap

- **Concluído:** Pipeline PM7 estável — 2.222 moléculas OK
- **Concluído:** Primeiro modelo ML treinado e aprovado — RMSE 3,54 kcal/mol, R² 0,997
- **Agora:** Escalar para 3.000 → 5.000 moléculas OK para re-treino com maior volume
- **Próximo:** Re-treino com validação cruzada k-fold, análise de moléculas outlier
- **Futuro:** Hiperparametrização (grid/Bayesian), escalar para 29k moléculas → deploy via `api.py`

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
