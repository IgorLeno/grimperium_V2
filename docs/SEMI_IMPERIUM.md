# Semi-Imperium

Semi-Imperium é um programa delimitado dentro deste repositório para calcular e
catalogar a **entalpia padrão de formação** (`ΔHf°`) por métodos semiempíricos.
Ele possui pacote, comandos, sessão e armazenamento próprios. Não é uma tela
alternativa do fluxo de delta-learning do Grimperium e não aplica correção por
machine learning aos seus resultados.

## Início rápido

Com as dependências instaladas, inicie o programa por qualquer uma destas rotas:

```bash
poetry run semi-imperium
# ou, dentro do ambiente virtual:
python -m semi_imperium
```

A interface superior contém somente **CALCULATE**, **DATABASE** e **SETTINGS**.
O comando legado `poetry run grimperium` continua iniciando o Grimperium e seus
fluxos originais.

## O que o número significa

O resultado de cada Hamiltoniano é a `ΔHf°` calculada pelo MOPAC, em
**kcal/mol**, para a fase e as opções registradas na configuração efetiva do
run. A interface e o store JSON do Semi-Imperium usam kcal/mol; a nota de
conversão para kJ/mol do fluxo legado do Grimperium não se aplica a este store.

AM1, PM3 e PM7 são Hamiltonianos semiempíricos do MOPAC. Eles são solicitações
**independentes**: cada um otimiza os conformeros selecionados na sua própria
superfície de energia, escolhe seu próprio candidato e possui sua própria
assinatura, estado de verificação e `ΔHf°`. O menor candidato segundo AM1 não é
presumido como o menor segundo PM3 ou PM7.

O Semi-Imperium não transforma automaticamente esses valores em referência
experimental ou em resultado de alta precisão. A proveniência deve acompanhar
qualquer número exportado ou comparado.

## CREST: geração não é seleção

Há duas decisões distintas antes do MOPAC:

1. **Geração do conjunto.** Com CREST ligado, o programa faz uma busca
   conformacional e recebe um ensemble. Com CREST desligado, ele gera uma única
   estrutura 3D inicial com RDKit; isso é um ponto de partida, não uma busca
   conformacional equivalente.
2. **Seleção do conjunto finito.** Uma estratégia ordena ou filtra o ensemble e
   limita quais estruturas cada Hamiltoniano realmente otimizará.

O CREST explora conformações e fornece geometrias e energias de busca. Ele não
fornece a `ΔHf°` final do Semi-Imperium e não decide sozinho qual mínimo do
MOPAC será aceito. Ligar ou desligar o CREST altera a evidência científica e,
por isso, participa da assinatura do cálculo.

### Energy Top-N (padrão)

`crest_energy_top_n` ordena o ensemble pela energia reportada pela busca e leva
ao MOPAC no máximo os `N` primeiros. O padrão é **Top-10**. Uma janela de energia
opcional pode reduzir ainda mais o conjunto antes do corte de `N`.

Em uma rota sem CREST há somente uma estrutura e nenhuma energia de busca a
inventar; essa estrutura única é selecionada com essa ausência registrada.

### CONFPASS (experimental)

`confpass_prioritization` é **EXPERIMENTAL** e precisa ser escolhido
explicitamente. O adaptador entrega o ensemble completo ao backend antes do
corte Top-N e preserva ordem de átomos, coordenadas, conectividade e proveniência
na conversão XYZ para SDF.

Limitações atuais:

- a priorização ainda não foi validada contra uma referência científica deste
  projeto;
- é necessário fornecer um backend CONFPASS e uma topologia compatível;
- o shell mostra a opção para inspeção/configuração, mas a composição de
  produção atual não injeta backend nem topologia CONFPASS; escolhê-la falha de
  forma explícita, portanto o caminho só é exercitável hoje por integração de
  API;
- a classe de completude PAS é metadado consultivo: ela pode ser registrada,
  mas nunca constitui evidência que justifique a seleção;
- um resultado marcado como experimental não deve ser apresentado como
  equivalente à estratégia padrão validada apenas porque o workflow terminou.

## Otimização e verificação de mínimo

Para cada Hamiltoniano, o MOPAC otimiza separadamente todo o conjunto finito
selecionado. Uma otimização convergida produz primeiro uma `ΔHf°` **provisória**.
Quando a política exige mínimo, o workflow executa `FORCE` na geometria
otimizada e examina as frequências, tratando modos triviais projetados e a
região numérica de baixa frequência de forma explícita.

Os estados têm significado científico, não apenas operacional:

- `verified / minimum_confirmed`: a geometria examinada foi confirmada como
  mínimo local e pode fornecer a `ΔHf°` final daquele cálculo;
- `saddle / saddle_point`: foi detectado ponto de sela; a energia provisória é
  preservada para diagnóstico, mas não promovida a valor final verificado;
- `unverified`: houve uma geometria/energia utilizável, porém a verificação não
  confirmou um mínimo;
- `failed`: não houve resultado científico utilizável.

A recuperação de sela por deslocamento ao longo de modo normal é limitada por
um orçamento registrado na assinatura. Falha de `FORCE`, frequência ilegível ou
esgotamento desse orçamento não autorizam o programa a inventar um mínimo.

O executor científico implementado nesta versão aceita
`require_minimum`. Embora o modelo de configuração ainda enumere `none` e
`advisory`, essas políticas não estão compostas no workflow executável e são
rejeitadas. Para um cálculo real, selecione `require_minimum` em SETTINGS; não
interprete a presença das outras opções como suporte operacional.

### O que a verificação não prova

Confirmar frequências reais demonstra, dentro do modelo e dos critérios usados,
um **mínimo local** da geometria examinada. Não prova o mínimo global da
molécula. A busca explora um espaço finito, a seleção otimiza apenas um
subconjunto finito e tanto CREST quanto MOPAC possuem aproximações e limites de
amostragem/convergência. Portanto, "menor entre os candidatos verificados" não
é sinônimo de "mínimo global matematicamente provado".

## Assinaturas, reutilização e banco local

O banco do Semi-Imperium é um store JSON próprio, por padrão em
`data/semi_imperium/`:

```text
data/semi_imperium/
├── runs/<run_id>/run.json
└── calculations/<molecule_id>/<signature>/<calculation_id>.json
```

Ele não escreve nos CSVs de batch nem nos datasets do Grimperium. Cada run
registra a configuração efetiva e a proveniência, incluindo identidade do
método/propriedade, versões do Semi-Imperium e do Grimperium, host/origem e
notas quando fornecidas. O modelo aceita versões dos programas externos, mas o
fluxo interativo atual não as preenche automaticamente; um mapa vazio significa
"não capturado", não uma versão inferida. Cada cálculo preserva também a
identidade molecular resolvida, estado, verificação, artefatos relativos e
timestamps.

A reutilização requer a mesma identidade molecular e a mesma assinatura
SHA-256 de configuração científica. Entram na assinatura:

- busca conformacional, inclusive CREST ligado/desligado;
- estratégia e limites da seleção;
- Hamiltoniano e opções científicas do MOPAC;
- política e orçamento de verificação de mínimo.

Threads, timeouts, caminhos e executáveis são detalhes de execução e não entram
na assinatura. Em uma solicitação AM1+PM3+PM7, a reutilização é decidida por
Hamiltoniano: um PM7 compatível não autoriza reutilizar AM1 ou PM3. Recalcular
cria outro run e outro registro; um resultado terminal anterior não é
sobrescrito silenciosamente.

## Entradas e sessões

CALCULATE aceita nome químico ou SMILES e mantém uma tabela com uma ou várias
moléculas. Nome químico precisa ser resolvido para uma estrutura; SMILES passa
por canonicalização e validação local. Entrada inválida, nome sem resolução,
ambiguidade e incompatibilidade de carga/multiplicidade ficam explícitos e não
seguem para os executáveis. A API suporta escolher candidato e recuperar uma
resolução com SMILES manual, mas a tela CALCULATE atual expõe apenas edição ou
remoção da linha; não prometa as recuperações da API como ações da interface.

Cada linha controla independentemente CREST e o conjunto AM1/PM3/PM7. Antes de
escrever ou executar, a revisão mostra resultados compatíveis e oferece as
ações possíveis: reutilizar, calcular somente Hamiltonianos ausentes ou
recalcular tudo.

### Composição externa atual

Readiness comprova que o executável está no `PATH`; não comprova que um adapter
foi ligado ao executor. O workspace padrão desta versão cria o executor sem um
runner CREST. Assim, uma linha com CREST ligado falha explicitamente antes de
executar em vez de cair silenciosamente na rota RDKit. Para integrar CREST é
necessário injetar um runner pela API; o smoke real abaixo faz isso de modo
visível. Com CREST desligado, a rota RDKit+MOPAC é utilizável após configurar
`require_minimum`.

## Testes e evidência científica

Os testes automatizados normais usam backends em memória ou runners que
produzem artefatos controlados. Eles demonstram contratos, roteamento, estados,
persistência e tratamento de falhas; **não são cálculos científicos reais** e
não validam a exatidão química de CREST, MOPAC ou CONFPASS.

O smoke real, pequeno e opt-in, fica em
`tests/integration/test_semi_imperium_real_smoke.py`. Ele só deve ser ativado em
uma máquina que possua CREST e MOPAC:

```bash
SEMI_IMPERIUM_REAL_SMOKE=1 poetry run pytest \
  tests/integration/test_semi_imperium_real_smoke.py -v -s
```

Sem a variável, ou sem qualquer binário, o teste informa `SKIPPED`; mocks não
substituem a ausência das ferramentas. O resultado observado na execução de
desenvolvimento desta versão está registrado em
[`SEMI_IMPERIUM_SMOKE.md`](SEMI_IMPERIUM_SMOKE.md).

