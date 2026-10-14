# Plano: espelho das tabelas no Google Sheets

Status: proposta (2026-10-12). Nada implementado ainda.

## Decisão

A planilha online é um **espelho somente leitura**. Os CSVs em `data/`
continuam sendo a fonte da verdade; edições feitas no Sheets nunca voltam
para o pipeline. Isso evita divergência entre o que foi calculado e o que
está publicado.

## Escopo proposto

| Aba da planilha | Origem | Observação |
|-----------------|--------|------------|
| `Experimental` | `data/experimental.csv` | **Excluir ou isolar linhas NIST**: os termos do NIST SRD 69 são “all rights reserved” (ver `NIST_TERMS` em `scripts/build_experimental_hf.py`). Publicar só fontes com licença compatível, ou manter a planilha privada. |
| `CBS` | `data/cbs_referencia.csv` | ThermoCBS é CC-BY 4.0: manter a citação no cabeçalho. ~27,8 mil linhas cabem no limite de 10 milhões de células. Sem a coluna `xyz` (pesada). |
| `PM7` | `data/pm7_fila_pipeline.csv` | Rotulada como legado `crest_pm7`. |
| `Comparação` | `data/experimental_vs_cbs.csv` | Mesma regra de licença da aba Experimental. |
| `Leia-me` | gerada | Data da publicação, sha256 de cada CSV, citações e aviso “somente leitura”. |

## Mecânica

1. Script `scripts/publish_sheets.py` (CLI, nunca chamado pelo dashboard
   automaticamente): lê os CSVs, substitui cada aba inteira e grava a aba
   `Leia-me` com hashes e data.
2. Autenticação por **conta de serviço** do Google Cloud; a chave fica fora do
   Git (variável de ambiente apontando para o JSON). A planilha é compartilhada
   com o e-mail da conta de serviço como editora e com leitores como “leitor”.
3. Dependência opcional nova (`gspread` ou `google-api-python-client`) em um
   extra `sheets`, para não pesar a instalação principal.
4. Depois do back-end de jobs, a publicação pode virar uma rotina disparada
   pela aba Rotinas.

## Decisões pendentes (do usuário)

- Planilha pública (link) ou restrita a e-mails?
- Valores NIST: omitir, ou publicar em planilha privada separada?
- Frequência: manual após cada build de dataset, ou agendada?
- Qual conta Google será dona da planilha?
