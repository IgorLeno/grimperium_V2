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
| **Delta Chemistry** (hub, tela principal) | Rotinas, Resumos e Sobre — nada específico de uma propriedade |
| **Grimperium** (módulo) | Entalpia de formação no estado padrão (PM7 corrigido com CBS) |

O menu superior mostra só as abas do item ativo: as telas do Grimperium
(Resumo, Dados, Calcular, Treinar) existem apenas dentro do módulo. O botão
**Menu**, que oculta/mostra a barra lateral, fica numa barra abaixo do
cabeçalho; o cabeçalho guarda só a marca e o menu superior.

Para adicionar um módulo, crie o componente de telas (modelo:
`src/GrimperiumModule.tsx`) e inclua uma entrada em `MODULES`
(`src/modules.ts`) com nome, propriedade e abas. O hub é `HUB` no mesmo
arquivo (componente `src/DeltaModule.tsx`).

## Telas do hub Delta Chemistry

| Menu | O que mostra |
|------|--------------|
| **Rotinas** (primeira) | Menu de rotinas e assistente em 5 passos, onde se escolhe a propriedade. É o mesmo componente da aba *Calcular* do Grimperium (`src/RoutinesView.tsx`) e usa o mesmo armazenamento local, então as rotinas aparecem nas duas. |
| **Resumos** | Compilado dos resultados obtidos pelas rotinas. Ainda sem back-end de jobs não há resultados: a tela mostra o estado vazio. |
| **Sobre** (última) | Explica o hub: a ideia do Δ-learning em três passos, como o hub se organiza, módulos disponíveis, escopo atual e glossário. |

## Telas do módulo Grimperium

Menu superior com submenu (definido em `src/modules.ts`):

| Menu | Submenu | O que mostra | Fonte |
|------|---------|--------------|-------|
| **Resumo** | — | Resumo em texto com números ao vivo, pipeline em 5 passos, mapa das três tabelas, fila CREST + PM7, store Semi-Imperium, ideia do Δ-learning e glossário. | `/api/dataset/progress`, perfis das tabelas, comparação |
| **Dados** | **Experimental** | O que é a tabela; distribuições; comparação CBS × experimental (concordância, paridade, erro por tamanho e fonte); tabela completa pesquisável. | `data/thermo_exp_hf.csv`, `data/thermo_exp_hf_vs_cbs.csv` |
| | **CBS** | O que é a referência CBS; entalpia, tamanho, elementos, correção CBS − B3LYP, spin/carga, SMILES repetidos; tabela completa. | `data/thermo_cbs_chon.csv` |
| | **PM7** | Tabela legada `crest_pm7`: andamento, erro Δ = CBS − PM7, paridade, resultado do CREST/MOPAC, tempos; tabela completa. | `data/thermo_pm7.csv` |
| **Calcular** | — | Menu de rotinas e assistente em 5 passos (nome, propriedade, pipeline, configuração, saída) com parâmetros lidos do Semi-Imperium. | `/api/routines/catalog` |
| **Treinar** | — | Explica as etapas do treino do modelo Δ e o portão de qualidade. Ainda sem execução: depende do back-end de jobs. | — |

Cada gráfico traz um bloco **“Como ler”** em linguagem simples. Os CSVs CBS
e PM7 ficam fora do Git; sem eles, a aba correspondente explica onde copiar o
arquivo. Tema claro/escuro segue o sistema operacional.

## Rotinas (painel de controle)

O dashboard **não executa cálculos**: ele monta a especificação da rotina
(`RoutineSpec`, versão 1, em `src/routines.ts`) e, quando o back-end de jobs
existir, apenas a dispara. Neste momento:

- rotinas ficam salvas no armazenamento local do navegador e podem ser
  exportadas em JSON;
- o botão **Disparar** fica desativado (`dispatch_enabled: false` no catálogo);
- o catálogo (`GET /api/routines/catalog`) lê padrões e opções de
  `semi_imperium` (CREST, seleção/CONFPASS, MOPAC, RDKit) e verifica se
  `crest`/`mopac` estão no PATH, sem executar nada.

## API (somente leitura)

| Rota | Conteúdo |
|------|----------|
| `GET /api/health` | Saúde da API |
| `GET /api/dataset/progress` | Estados do CSV PM7 e store Semi-Imperium |
| `GET /api/experimental/summary` · `/comparison` | Fontes/flags e pares CBS × experimental |
| `GET /api/tables/{exp,cbs,pm7}/profile` | Distribuições: tamanho, composição CHON, histogramas (99% central; extremos contados à parte), categorias, erro por tamanho e nuvem de paridade (até 4.000 pontos) |
| `GET /api/tables/{exp,cbs,pm7}/rows?q=&offset=&limit=&sort=&desc=` | Linhas paginadas (máx. 200) com busca e ordenação |
| `GET /api/routines/catalog` | Propriedades, modos de entrada, programas e saídas das rotinas |

Todas as rotas são `GET`; `tests/webapi/test_app.py` garante isso.

## Verificação

```bash
cd apps/dashboard
npm run lint && npm run typecheck && npm run test && npm run build
cd ../.. && poetry run pytest tests/webapi -q
```
