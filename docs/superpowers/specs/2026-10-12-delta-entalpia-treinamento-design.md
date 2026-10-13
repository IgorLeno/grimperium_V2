# Delta Entalpia: rotina de treinamento (Semi-Imperium) acionada pelo dashboard

Data: 2026-10-12
Status: rascunho para aprovação. Nada implementado.
Plano-mãe: `docs/plans/2026-10-03-hub-restructure.md`.

## Intenção

O dashboard passa a ser o **painel de acionamento e apresentação** das rotinas.
As rotinas são construídas e mantidas só no back-end Python e rodam também pela
linha de comando. O usuário não monta rotinas no dashboard; escolhe uma rotina
pronta, aciona, acompanha e lê os resultados.

O módulo de entalpia de formação passa a se chamar **Delta Entalpia** e reúne três
rotinas:

1. **Treinamento (Semi-Imperium):** calcula, com o protocolo canônico, as
   moléculas PM7 que servirão de dados de treino. **Este documento cobre só esta.**
2. **Delta Learning:** aprende `ΔHf_CBS − ΔHf_PM7` sobre as moléculas calculadas
   (estratégia do Kemperium). Fora de escopo aqui.
3. **Cálculo (Grimperium):** recebe o nome de uma molécula, reaproveita o valor se
   ela já estiver no banco, senão roda o protocolo, e aplica a correção do modelo.
   Fora de escopo aqui.

Os programas devem ser modulares: acoplar e desacoplar etapas sem reescrever o
resto. Por isso o dashboard só conhece "rotina, parâmetros, estado, resultados".

## Objetivos deste projeto

- Rotina de treinamento que **completa** o conjunto de ~29,5 mil moléculas de
  `data/thermo_pm7.csv` (29.568 linhas) com o protocolo 2, de forma retomável.
- Protocolo por molécula: SMILES → RDKit → CREST (GFN2-xTB) → seleção de
  conformeros → otimização MOPAC PM7 → ΔHf PM7 → frequências vibracionais que
  confirmam mínimo local → registro do valor e da geometria final.
- As **8.232 moléculas OK** do protocolo 1 (`crest_pm7`) são mantidas como
  resultado **preliminar**, marcadas por uma coluna que indica o protocolo 1.
  Ficam sem geometria nem frequências. Serão recalculadas pelo protocolo 2
  numa etapa futura, quando houver moléculas suficientes no protocolo 2.
- As 4.427 moléculas `Selected` voltam a `Pending` e entram na fila do protocolo 2.
- Guardar **só a geometria final do conformer selecionado**, fora do CSV, numa
  única tabela SQLite indexada por `mol_id`, de leitura rápida e exportável
  para XYZ.
- Acionar, parar e acompanhar a rotina pelo dashboard e pela CLI, com os mesmos
  comandos por baixo.
- Aba de treinamento com o progresso: quantas moléculas estão no protocolo 2,
  quantas preliminares (protocolo 1), pendentes, com falha, em execução.
- Reorganizar o dashboard no módulo **Delta Entalpia**. Remover o construtor de
  rotinas (aba Rotinas e seu catálogo).

**Critério de sucesso:** pelo dashboard, iniciar a rotina e ver ao menos uma
molécula concluir com ΔHf PM7, frequências confirmando mínimo, geometria
gravada no SQLite, o CSV atualizado e o progresso refletido na tela.

## Não-objetivos

- Delta Learning e Grimperium (cálculo com correção do modelo).
- Recalcular as 8.232 do protocolo 1.
- Rodar o lote completo de ~21 mil moléculas como parte desta entrega. Entrega-se
  a rotina; o lote é decisão do usuário (ver pergunta aberta).
- Sincronização recorrente com o Google Planilhas (envio único é projeto à parte).
- Reaproveitar o server/worker congelado (`src/grimperium/server`, `worker`).
- Servidor multiusuário ou multi-máquina.

## Anti-objetivos

- Perder ou sobrescrever os valores das 8.232 moléculas preliminares.
- Gravar coordenadas na planilha/CSV conjunto.
- Deixar o usuário compor rotinas no dashboard.
- Tratar um valor sem confirmação de mínimo como resultado do protocolo 2.

## Estado atual verificado (2026-10-12)

- `data/thermo_pm7.csv`: 29.568 linhas; status OK 8.232, Selected 4.427,
  Pending 16.876, Skip 29, Rerun 3, Running 1.
- As 8.232 OK têm ΔHf PM7 e propriedades MOPAC na planilha, mas **não** têm
  geometria: `data/molecules_pm7/conformer_details/` tem 2.862 JSONs só com
  energias, e o `geometry_file` aponta para `/tmp/crest_pm7/...`, que não existe
  mais. Não há frequências.
- Já existe o pacote `src/semi_imperium/` (~11 mil linhas, CLI `semi-imperium`):
  resolução de moléculas (inclui PubChem), CREST por subprocesso, seleção de
  conformeros (Energy Top-N e porte experimental do CONFPASS, filtro de
  dobramento opcional), otimização MOPAC com verificação de mínimo
  (`mopac/workflow.py`, `force_parser.py`), persistência própria em
  `data/semi_imperium/` (hoje 1 molécula).
- A API do dashboard (`src/grimperium/webapi/app.py`, 91 linhas) só tem rotas GET
  de leitura. O catálogo de rotinas e o construtor ficam em
  `webapi/catalog.py` e `apps/dashboard/src/routines.ts`, `RoutinesView.tsx`.
- O dashboard já tem o hub Delta Chemistry com módulos (`modules.ts`).
- `xtb`, `crest` e `mopac` estão instalados. O CONFPASS original está em
  `~/Estagio/confpass-original` (exige Python 3.8); o porte usado em produção é
  `src/semi_imperium/conformers/confpass_port.py`.

## Restrições

- Máquina única, 16 núcleos, 31 GB; usuário único.
- A autoridade do esquema do CSV é `src/grimperium/crest_pm7/batch/csv_manager.py`
  e `tests/test_csv_schema.py`. A nova coluna de protocolo deve passar por eles.
- Portões de qualidade do `CLAUDE.md` (black, ruff, mypy --strict, pytest;
  lint/typecheck/testes/build do dashboard). `CHANGELOG.md` em `[Unreleased]`.
- Verificação científica deve cruzar loaders, fusão, ML e pipeline, não só o
  módulo editado.

## Decisões do usuário

- Manter as 8.232 como preliminares; recalcular depois.
- `Selected` volta a Pending.
- Guardar só a geometria final do conformer selecionado, em tabela única SQLite
  indexada por `mol_id`, fora do CSV.
- Rotinas só no back-end; dashboard só aciona e apresenta; CLI continua.
- Módulo renomeado para Delta Entalpia, com Treinamento, Delta Learning e Cálculo.

## Escolhas minhas (revisar)

- Acionamento: a API do dashboard inicia a rotina como **subprocesso** do mesmo
  comando da CLI; o estado vive no SQLite e no CSV, não em memória do servidor.
  Parar = sinalizar o subprocesso; retomar = reiniciar (a rotina pula o que já
  está feito).
- Paralelismo: N moléculas em paralelo, configurável, com padrão conservador
  para não saturar os 16 núcleos (CREST e MOPAC já usam threads).
- Ordem da fila: moléculas menores (`nheavy`) primeiro.
- Dashboard: Delta Learning e Cálculo aparecem como abas "em breve" até seus
  projetos existirem.
- A coluna de protocolo aceita ao menos "1" (preliminar) e "2"; o nome exato e
  onde o SQLite fica em `data/semi_imperium/` são do construtor.
- O que o SQLite e a persistência já existente em `data/semi_imperium/` guardam
  não deve duplicar dado: o construtor decide se estende a persistência atual.

## Escolhas minhas após revisão (revisar)

- **Parâmetros do protocolo 2:** os padrões atuais de `src/semi_imperium/`
  (CREST, MOPAC PM7, carga e multiplicidade vindas da linha do CSV). Todos são
  gravados na assinatura do protocolo. Só a estratégia de seleção de conformeros
  espera a decisão do usuário (pergunta aberta); a rotina a recebe como
  parâmetro e a grava na assinatura.
- **Frequência imaginária:** usar a recuperação limitada que o fluxo MOPAC de
  `semi_imperium` já tem para ponto de sela. Se a verificação continuar falhando,
  a molécula vira `Failed`, sem valor de protocolo 2, e não bloqueia a fila.
- **Fila:** entram `Pending`, `Selected` (resetado), `Rerun` e `Running` órfão
  (execução interrompida). `Skip` continua fora. `OK` do protocolo 1 continua
  fora, pois é preliminar.
- **Migração do CSV:** entregue como comando com `--dry-run` e backup
  automático. Só é executada no `data/thermo_pm7.csv` real depois de o usuário
  confirmar no momento da execução. O CSV migrado não entra no PR sem essa
  confirmação.
- **Fonte de verdade:** o CSV continua sendo o conjunto de dados (os loaders de
  ML leem dele) e tem **um único escritor**, o processo coordenador da rotina.
  Os trabalhadores paralelos devolvem resultados a ele. O SQLite guarda
  geometrias e o estado de execução; o dashboard lê os dois.

## Deixado ao construtor

- Nomes de colunas, tabelas, endpoints e comandos.
- Formato do estado de execução e da fila.
- Estratégia de retomada após queda no meio de uma molécula (status `Running`).
- Como migrar o CSV (marcação do protocolo 1 e reset de `Selected`) com backup.

## Decisão do usuário (2026-10-12): seleção de conformeros

O protocolo 2 usa **CONFPASS + filtro de dobramento**, como em
`2026-10-03-conformer-selection-spec.md`, e não o Energy Top-N. A assinatura do
protocolo 2 grava essa estratégia e os parâmetros calibrados do filtro
(`reports/folding_calibration/report.md`). O porte do CONFPASS
(`confpass_port.py`) ainda é experimental. O plano-mãe diz que o lote não roda
até o protocolo estar fechado, então o construtor deve validar o porte contra a
especificação antes de recomendar o lote completo. Quem decide rodá-lo é o usuário.
