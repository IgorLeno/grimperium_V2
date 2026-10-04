# Calibração do filtro de dobramento extremo

Gerado em 2026-10-04T11:32:25.116734+00:00 por `scripts/calibrate_folding_filter.py`.
Especificação: `docs/plans/2026-10-03-conformer-selection-spec.md`, seção 4.1.

## Entrada

- Arquivo: `data/thermo_cbs_chon_v2.csv`
- SHA-256: `20693a0d7c2b6c66052a263687124321f761888e4773fe2edf72e61a4831ffa8`
- Linhas: 27760

Percepção de topologia (geometria vs SMILES da linha):

| status | linhas |
|---|---|
| ok | 27746 |
| topology_mismatch | 12 |
| disconnected_geometry | 2 |

Calibração usa apenas `ok` (27746 linhas).

## Definição

- `fold_contacts`: pares de átomos pesados com distância topológica ≥ k e distância < f · (r_vdW,i + r_vdW,j).
- Raios de van der Waals: tabela periódica do RDKit (2026.03.3): C 1,70; N 1,60; O 1,55 Å.
- Excluídos: pares doador/aceitador O/N com H do doador a < 2.5 Å do aceitador (ligação H clássica).
- Descartada se `fold_contacts ≥ c_min`. `rg_ratio` fora desta calibração.
- Orçamento: ≤ 1.00% das referências descartadas.
- Coluna `fechada` = multiplicidade 1 e carga 0 (população da validação).

## Proposta

**k = 6, f = 1.00, c_min = 1** — descarta 20 referências (0.07%; camada fechada 0.07%).

Regra de escolha (no código, `propose`): só k ≥ 6 e f ≤ 1.00; dentro do orçamento, maior f, depois menor k, depois menor c_min. k = 4 sinaliza contatos 1,5 forçados por substituição (anisóis com tBu orto) e k = 5 contatos 1,6 em grupos conjugados planos (oximas, enaminas); só k ≥ 6 sinaliza pontas de cadeia se encontrando. f > 1 é proximidade, não contato extremo.

Taxa de descarte por `nheavy` com a proposta:

| nheavy | linhas | descartadas | taxa |
|---|---|---|---|
| 1 | 1 | 0 | 0.00% |
| 2 | 5 | 0 | 0.00% |
| 3 | 16 | 0 | 0.00% |
| 4 | 31 | 0 | 0.00% |
| 5 | 183 | 0 | 0.00% |
| 6 | 631 | 0 | 0.00% |
| 7 | 1586 | 0 | 0.00% |
| 8 | 2159 | 0 | 0.00% |
| 9 | 5390 | 2 | 0.04% |
| 10 | 8911 | 8 | 0.09% |
| 11 | 7395 | 7 | 0.09% |
| 12 | 391 | 1 | 0.26% |
| 13 | 301 | 1 | 0.33% |
| 14 | 511 | 0 | 0.00% |
| 15 | 122 | 1 | 0.82% |
| 16 | 80 | 0 | 0.00% |
| 17 | 7 | 0 | 0.00% |
| 18 | 14 | 0 | 0.00% |
| 19 | 5 | 0 | 0.00% |
| 20 | 6 | 0 | 0.00% |
| 22 | 1 | 0 | 0.00% |

## Fronteira (maior f dentro do orçamento por k, c_min)

| k | c_min | f | taxa | fechada |
|---|---|---|---|---|
| 4 | 1 | 0.85 | 0.44% | 0.42% |
| 4 | 2 | 0.90 | 0.36% | 0.37% |
| 4 | 3 | 0.95 | 0.46% | 0.48% |
| 4 | 4 | 1.00 | 0.29% | 0.29% |
| 5 | 1 | 0.95 | 0.27% | 0.26% |
| 5 | 2 | 1.05 | 0.37% | 0.37% |
| 5 | 3 | 1.10 | 0.28% | 0.28% |
| 5 | 4 | 1.10 | 0.10% | 0.10% |
| 6 | 1 | 1.10 | 0.87% | 0.89% |
| 6 | 2 | 1.10 | 0.16% | 0.17% |
| 6 | 3 | 1.10 | 0.06% | 0.06% |
| 6 | 4 | 1.10 | 0.01% | 0.01% |
| 7 | 1 | 1.10 | 0.26% | 0.27% |
| 7 | 2 | 1.10 | 0.04% | 0.04% |
| 7 | 3 | 1.10 | 0.01% | 0.01% |
| 7 | 4 | 1.10 | 0.01% | 0.01% |

## Grade completa (taxa de descarte, todas as linhas `ok`)

| k | f | c_min=1 | c_min=2 | c_min=3 | c_min=4 |
|---|---|---|---|---|---|
| 4 | 0.70 | 0.00% | 0.00% | 0.00% | 0.00% |
| 4 | 0.75 | 0.00% | 0.00% | 0.00% | 0.00% |
| 4 | 0.80 | 0.00% | 0.00% | 0.00% | 0.00% |
| 4 | 0.85 | 0.44% | 0.00% | 0.00% | 0.00% |
| 4 | 0.90 | 3.20% | 0.36% | 0.17% | 0.00% |
| 4 | 0.95 | 9.59% | 1.84% | 0.46% | 0.06% |
| 4 | 1.00 | 22.59% | 5.58% | 1.21% | 0.29% |
| 4 | 1.05 | 37.63% | 13.77% | 3.80% | 1.14% |
| 4 | 1.10 | 46.43% | 21.60% | 8.52% | 3.44% |
| 5 | 0.70 | 0.00% | 0.00% | 0.00% | 0.00% |
| 5 | 0.75 | 0.00% | 0.00% | 0.00% | 0.00% |
| 5 | 0.80 | 0.00% | 0.00% | 0.00% | 0.00% |
| 5 | 0.85 | 0.00% | 0.00% | 0.00% | 0.00% |
| 5 | 0.90 | 0.08% | 0.00% | 0.00% | 0.00% |
| 5 | 0.95 | 0.27% | 0.00% | 0.00% | 0.00% |
| 5 | 1.00 | 1.14% | 0.08% | 0.01% | 0.00% |
| 5 | 1.05 | 3.58% | 0.37% | 0.09% | 0.03% |
| 5 | 1.10 | 6.66% | 1.30% | 0.28% | 0.10% |
| 6 | 0.70 | 0.00% | 0.00% | 0.00% | 0.00% |
| 6 | 0.75 | 0.00% | 0.00% | 0.00% | 0.00% |
| 6 | 0.80 | 0.00% | 0.00% | 0.00% | 0.00% |
| 6 | 0.85 | 0.00% | 0.00% | 0.00% | 0.00% |
| 6 | 0.90 | 0.00% | 0.00% | 0.00% | 0.00% |
| 6 | 0.95 | 0.00% | 0.00% | 0.00% | 0.00% |
| 6 | 1.00 | 0.07% | 0.01% | 0.00% | 0.00% |
| 6 | 1.05 | 0.35% | 0.04% | 0.02% | 0.01% |
| 6 | 1.10 | 0.87% | 0.16% | 0.06% | 0.01% |
| 7 | 0.70 | 0.00% | 0.00% | 0.00% | 0.00% |
| 7 | 0.75 | 0.00% | 0.00% | 0.00% | 0.00% |
| 7 | 0.80 | 0.00% | 0.00% | 0.00% | 0.00% |
| 7 | 0.85 | 0.00% | 0.00% | 0.00% | 0.00% |
| 7 | 0.90 | 0.00% | 0.00% | 0.00% | 0.00% |
| 7 | 0.95 | 0.00% | 0.00% | 0.00% | 0.00% |
| 7 | 1.00 | 0.04% | 0.00% | 0.00% | 0.00% |
| 7 | 1.05 | 0.09% | 0.02% | 0.01% | 0.00% |
| 7 | 1.10 | 0.26% | 0.04% | 0.01% | 0.01% |

## Exemplos sinalizados pela proposta (top 25)

Ordem: mais contatos, depois menor razão vdW. Lista completa em `flagged.csv`.
Índices de átomo são 0-based na ordem do bloco `xyz`.

| mol_id | SMILES | nheavy | contatos | par | ligações | d (Å) | d/Σr_vdW |
|---|---|---|---|---|---|---|---|
| cbs_30466 | `O=C(O)CCCCCCC(=O)O` | 12 | 2 | 2-9 | 8 | 3.142 | 0.967 |
| cbs_14159 | `NCCC(=O)NCCC(=O)O` | 11 | 2 | 2-10 | 6 | 3.156 | 0.971 |
| cbs_06656 | `CC(=O)OCCOC(C)=O` | 10 | 2 | 1-8 | 6 | 3.211 | 0.988 |
| cbs_02480 | `O=C(O)CCOCCC(=O)O` | 11 | 2 | 2-10 | 8 | 3.069 | 0.990 |
| cbs_04291 | `[H]/N=C(\N)N/C(=C\[N+](=O)[O-])NC` | 11 | 1 | 5-10 | 6 | 2.939 | 0.933 |
| cbs_20663 | `COCCCN(C)CCO` | 10 | 1 | 0-9 | 8 | 3.186 | 0.980 |
| cbs_17052 | `C=CC[C@@H](C)/C=C(/C)CO` | 10 | 1 | 5-8 | 6 | 3.188 | 0.981 |
| cbs_11399 | `C/C(=C/CCC#N)CO` | 9 | 1 | 5-8 | 6 | 3.192 | 0.982 |
| cbs_30457 | `OCCOCCOCCOCCO` | 13 | 1 | 5-12 | 7 | 3.214 | 0.989 |
| cbs_20577 | `COCCNC(C)(C)CO` | 10 | 1 | 0-9 | 7 | 3.220 | 0.991 |
| cbs_12084 | `CO[C@@H](C)CO[C@@H](C)CO` | 10 | 1 | 0-7 | 7 | 3.225 | 0.992 |
| cbs_20287 | `CN(C)CCN(C)CCO` | 10 | 1 | 0-7 | 7 | 3.227 | 0.993 |
| cbs_11437 | `COCCCOCCO` | 9 | 1 | 0-8 | 8 | 3.234 | 0.995 |
| cbs_20120 | `OCCNCCNCCO` | 10 | 1 | 1-9 | 8 | 3.232 | 0.995 |
| cbs_11127 | `CC(=O)NCCOCC(=O)O` | 11 | 1 | 3-10 | 6 | 3.136 | 0.996 |
| cbs_11289 | `O=C[C@@H](O)[C@@H](O)[C@@H](O)CO` | 10 | 1 | 0-6 | 6 | 3.092 | 0.997 |
| cbs_05742 | `CC(=O)C(C)(C)CCC(=O)O` | 11 | 1 | 2-8 | 6 | 3.244 | 0.998 |
| cbs_07414 | `COC[C@H](C)OCC(=O)OC` | 11 | 1 | 1-7 | 6 | 3.095 | 0.999 |
| cbs_30329 | `CC(=O)OCC(COC(C)=O)OC(C)=O` | 15 | 1 | 9-12 | 6 | 3.245 | 0.999 |
| cbs_24293 | `CN(C)/[CH]=[C](/C#N)[CH][N+]([CH3])C` | 11 | 1 | 5-9 | 6 | 3.399 | 1.000 |

## Limites

- Só referências CBS-QB3 (uma geometria por linha); não mede quantos conformeros CREST dobrados o filtro pegaria.
- `r_max` (`rg_ratio`) exige ensembles reais; calibrado na validação.
- Topologia percebida da geometria; em produção vem da molécula resolvida.
