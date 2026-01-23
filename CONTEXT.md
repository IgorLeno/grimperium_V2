# Contexto do Projeto Grimperium

## Visao geral
- Framework de ML para predicao de propriedades termodinamicas com delta-learning.
- Stack: Python 3.10+, Rich, Questionary, Pandas, Pytest, Ruff, Black, MyPy.
- Fase atual: Phase C (CLI Interactive Application), BATCH 12 (correcoes criticas).

## Fluxo principal (CREST -> PM7 -> CSV)
1. Selecao de batch via CLI (batch_view/databases_view).
2. Processamento por molecula (crest_pm7/batch/execution_manager.py).
3. CREST gera conformers -> MOPAC/PM7 otimiza cada conformer.
4. PM7Result -> CSV (batch/csv_manager.py).
5. Enriquecimento do CSV (csv_enhancements.py) com:
   - abs_diff e abs_diff_% (|H298_cbs - H298_pm7|)
   - delta_1/2/3 (|H298_cbs - hof| para top 3 conformers)
   - conformer_selected (indice 1-based do menor delta)
   - Settings do batch (CREST/MOPAC)

## Comandos principais
- Testes: pytest tests/ -v --cov=src/grimperium
- Lint: ruff check src/ tests/
- Format: black src/ tests/
- Type check: mypy src/ --strict

## Pastas criticas
- src/grimperium/crest_pm7/ -> pipeline CREST/PM7
- src/grimperium/crest_pm7/batch/ -> batch, CSV, execucao
- src/grimperium/cli/ -> views e fluxo da interface
- tests/ -> unit + integration (fixtures em tests/conftest.py)

## Estado atual (jan/2026)
- Batch 12 em andamento, foco em correcoes do CSV e fluxo de batch.
- Arquivo de dados principal: data/thermo_pm7.csv.
- Problemas reportados:
  - abs_diff/abs_diff_% nao preenchidos (header com espaco: "abs_diff ").
  - delta_1 retornando 0 indevidamente (semantica correta e |H298_cbs - hof|).
  - num_conformers_selected menor que 3 mesmo com crest_conformers_generated >= 3.
  - Muitas linhas "Processing" sem info relevante durante execucao.

## Expectativas funcionais
- Conformers: se crest_conformers_generated >= 3, selecionar 3; so <3 se CREST gerar <3.
- Deltas:
  - delta_1/2/3 = |H298_cbs - hof| para top 3 conformers (ordenados por energia).
  - conformer_selected = indice (1-based) do menor delta.
- CSV: colunas normalizadas (sem espacos) e valores salvos apos cada molecula.

## Arquivos com ajustes recentes
- src/grimperium/crest_pm7/batch/csv_manager.py
- src/grimperium/crest_pm7/csv_enhancements.py
- src/grimperium/crest_pm7/conformer_selector.py
- src/grimperium/crest_pm7/batch/execution_manager.py
