# Etapa 2: fila de treinamento no SQLite, semeada do CBS

Decisões aprovadas: um protocolo só (protocolo 2); legado fora do caminho; CBS CSV
somente leitura; estado das moléculas no SQLite `data/semi_imperium/training.sqlite3`.

## Critério de seleção (um por SMILES)

`data/cbs_referencia.csv` tem 27.760 linhas e 27.335 SMILES únicos. Para cada SMILES
fica a linha de `conformer_rank_by_h298` mínimo (o conformer de referência de menor
H298); empate por `mol_id`. O CSV é só lido, nunca escrito.

## Itens

- [x] `queue.py` (tabela `training_queue` em `db.py`): semeadura idempotente, fila,
      `mark_running`, `record_success`, `record_failure`, `read_progress`
- [x] Testes novos de fila + adaptação de `test_runner`, `test_cli`
- [x] `runner.py`, `__main__.py`, `protocol.py` (sem `LEGACY_PROTOCOL`, sem `migrate`)
- [x] Remover `ledger.py`, `migration.py` e seus testes
- [x] webapi: sem `legacy`, `pm7_csv`, tabela `pm7`; progresso contra a fila
- [x] dashboard: sem `Pm7View`, sem "PM7 (legado)", sem card "PM7 concluído"
- [x] CHANGELOG `[Unreleased]` + docs próximos
- [x] Quality gates (black, ruff, mypy, pytest, vitest/build)
