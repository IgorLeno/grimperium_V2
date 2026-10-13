# Dashboard local

`index.html` é a entrada do Vite. Abrir esse arquivo com `file://` não executa o
aplicativo: o navegador não transforma os módulos TypeScript/React nem resolve
o caminho `/src/main.tsx`.

Para ver dados reais, inicie a API a partir da raiz do repositório:

```bash
poetry run uvicorn grimperium.webapi.app:create_app --factory --host 127.0.0.1 --port 8000
```

Confira `http://127.0.0.1:8000/api/health` (deve retornar `status: ok`). Em
outro terminal, inicie o frontend:

```bash
cd apps/dashboard
npm run dev
```

Abra o endereço HTTP mostrado pelo Vite (normalmente
`http://localhost:5173/`). O Vite encaminha `/api` para a porta 8000. Se a API
ou um CSV estiver indisponível, a tela mostra o erro e oferece recarga.

Em um aplicativo com suporte a preview de servidores, `.claude/launch.json`
define as duas configurações: inicie `grimperium-api` e depois
`grimperium-dashboard`. As portas são fixas porque o proxy do Vite aponta para
`127.0.0.1:8000`.

## Delta Chemistry, hub e módulos

O dashboard é o shell do **Delta Chemistry**: um hub de rotinas de cálculo de
propriedades com métodos rápidos corrigidos por machine learning (Δ-learning)
a partir de dados de métodos mais precisos. A **barra lateral** tem dois grupos:

| Item | O que é |
|------|---------|
| **Delta Chemistry** (hub, tela principal) | Resumos e Sobre — nada específico de uma propriedade |
| **Delta Entalpia** (módulo) | Entalpia de formação no estado padrão (PM7 corrigido com CBS): Resumo, Dados, Treinamento, Delta Learning e Cálculo |

O menu superior mostra só as abas do item ativo: as telas do Delta Entalpia
existem apenas dentro do módulo. O botão **Menu**, que oculta/mostra a barra
lateral, fica numa barra abaixo do cabeçalho; o cabeçalho guarda só a marca e o
menu superior.

As rotinas são construídas e mantidas no back-end Python e rodam também pela
linha de comando. O painel **só aciona e acompanha**; não monta rotinas.

Para adicionar um módulo, crie o componente de telas (modelo:
`src/DeltaEnthalpyModule.tsx`) e inclua uma entrada em `MODULES`
(`src/modules.ts`) com nome, propriedade e abas. O hub é `HUB` no mesmo
arquivo (componente `src/DeltaModule.tsx`).

## Telas do hub Delta Chemistry

| Menu | O que mostra |
|------|--------------|
| **Resumos** | Compilado dos resultados obtidos pelas rotinas. Ainda não há resultados consolidados: a tela mostra o estado vazio e aponta para o Treinamento. |
| **Sobre** (última) | Explica o hub: a ideia do Δ-learning em três passos, como o hub se organiza, módulos disponíveis, escopo atual e glossário. |

## Telas do módulo Delta Entalpia

Menu superior com submenu (definido em `src/modules.ts`):

| Menu | Submenu | O que mostra | Fonte |
|------|---------|--------------|-------|
| **Resumo** | — | Resumo em texto com números ao vivo, pipeline em 5 passos, mapa das três tabelas, fila CREST + PM7, store Semi-Imperium, ideia do Δ-learning e glossário. | `/api/dataset/progress`, perfis das tabelas, comparação |
| **Dados** | **Experimental** | O que é a tabela; distribuições; comparação CBS × experimental (concordância, paridade, erro por tamanho e fonte); tabela completa pesquisável. | `data/thermo_exp_hf.csv`, `data/thermo_exp_hf_vs_cbs.csv` |
| | **CBS** | O que é a referência CBS; entalpia, tamanho, elementos, correção CBS − B3LYP, spin/carga, SMILES repetidos; tabela completa. | `data/thermo_cbs_chon.csv` |
| | **PM7** | Tabela legada `crest_pm7`: andamento, erro Δ = CBS − PM7, paridade, resultado do CREST/MOPAC, tempos; tabela completa. | `data/thermo_pm7.csv` |
| **Treinamento** | — | Rotina Semi-Imperium (protocolo 2): barra de andamento com denominador explícito, iniciar/parar com moléculas em paralelo e limite, falhas recentes com o motivo. Atualiza a cada 5 s. | `/api/training/status`, `/api/training/start`, `/api/training/stop` |
| **Delta Learning** | — | Explica as etapas do treino do modelo Δ e o portão de qualidade. Ainda sem execução no painel. | — |
| **Cálculo** | — | "Em breve": cálculo de uma molécula pelo nome (Grimperium), depende do Delta Learning. | — |

Cada gráfico traz um bloco **“Como ler”** em linguagem simples. Os CSVs CBS
e PM7 ficam fora do Git; sem eles, a aba correspondente explica onde copiar o
arquivo. Tema claro/escuro segue o sistema operacional.

## Treinamento (Delta Entalpia)

O painel aciona e acompanha a rotina; o cálculo roda no back-end e também pelo terminal:

```bash
poetry run python -m semi_imperium.training status
poetry run python -m semi_imperium.training migrate            # ensaio
poetry run python -m semi_imperium.training migrate --apply    # grava (cria cópia .pre-protocol-*.bak)
poetry run python -m semi_imperium.training run --workers 3 --limit 10
poetry run python -m semi_imperium.training stop
```

A migração é feita uma vez: marca as moléculas do protocolo 1 e devolve
`Selected` à fila. Enquanto o CSV não for migrado, o botão **Iniciar
treinamento** fica desabilitado e a tela mostra o comando. A API inicia a rotina
como subprocesso do mesmo comando (`run`); o estado fica em
`data/semi_imperium/training.sqlite3`. `GRIMPERIUM_ROOT` aponta a API para outra
pasta no formato do repositório (usado para ensaiar em uma cópia dos dados).

## API

| Rota | Conteúdo |
|------|----------|
| `GET /api/health` | Saúde da API e se ela pode acionar a rotina (`can_dispatch`) |
| `GET /api/dataset/progress` | Estados do CSV PM7 e store Semi-Imperium |
| `GET /api/experimental/summary` · `/comparison` | Fontes/flags e pares CBS × experimental |
| `GET /api/tables/{exp,cbs,pm7}/profile` | Distribuições: tamanho, composição CHON, histogramas (99% central; extremos contados à parte), categorias, erro por tamanho e nuvem de paridade (até 4.000 pontos) |
| `GET /api/tables/{exp,cbs,pm7}/rows?q=&offset=&limit=&sort=&desc=` | Linhas paginadas (máx. 200) com busca e ordenação |
| `GET /api/training/status` | Andamento (CSV), estado do job, geometrias guardadas e falhas recentes |
| `POST /api/training/start` · `/stop` | Inicia (`workers` 1–8, `limit` opcional) ou pede parada; só origens locais; 409 se o CSV não foi migrado ou já há job vivo |

Fora `POST /api/training/start` e `/stop`, todas as rotas são `GET`; `tests/webapi/test_app.py` garante isso.

## Verificação

```bash
cd apps/dashboard
npm run lint && npm run typecheck && npm run test && npm run build
cd ../.. && poetry run pytest tests/webapi -q
```
