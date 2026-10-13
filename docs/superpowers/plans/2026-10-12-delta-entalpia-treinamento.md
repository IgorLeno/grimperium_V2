# Delta Entalpia: rotina de treinamento (Semi-Imperium) — Plano de implementação

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rotina de treinamento retomável (protocolo 2 do Semi-Imperium) que completa as 29.568 moléculas de `data/thermo_pm7.csv`, acionável pela CLI e pelo dashboard, com a geometria final de cada molécula em SQLite.

**Architecture:** Um pacote novo `src/semi_imperium/training/` reúne protocolo, migração do CSV, ledger do CSV (único escritor), tabela SQLite de geometrias, estado do job e um coordenador que despacha moléculas para um `ProcessPoolExecutor`; cada processo roda o `ScientificCalculationExecutor` que já existe. A CLI é `python -m semi_imperium.training`. A API do dashboard ganha rotas `POST /api/training/start|stop` e `GET /api/training/status` que lançam e acompanham esse mesmo subprocesso. O dashboard vira o módulo **Delta Entalpia** (Treinamento funcional; Delta Learning e Cálculo como "em breve") e perde o construtor de rotinas.

**Tech Stack:** Python 3.14 (poetry), pandas, sqlite3 (stdlib), FastAPI/pydantic, pytest; React 19 + TypeScript + Vite + vitest. Programas externos já instalados: `xtb`, `crest`, `mopac`.

**Spec:** `docs/superpowers/specs/2026-10-12-delta-entalpia-treinamento-design.md` (aprovado em 2026-10-12). Ler antes de começar.

## Global Constraints

- Sempre `poetry run ...` (nenhum venv ativo por padrão).
- Portões de qualidade antes do fim de cada tarefa que mexe em Python: `poetry run black src/ tests/`, `poetry run ruff check src/ tests/`, `poetry run mypy src/ --strict` (rodar antes `poetry run python scripts/patch_rdkit_stubs.py`), `poetry run pytest <alvo> -v`. No dashboard: `npm run lint`, `npm run typecheck`, `npm test`, `npm run build` em `apps/dashboard`.
- A autoridade do esquema do CSV é `src/grimperium/crest_pm7/batch/csv_manager.py` e `tests/test_csv_schema.py`. O CSV só é escrito por `BatchCSVManager` (escrita atômica com `.bak`).
- `data/thermo_pm7.csv` é ignorado pelo git (`*.csv`); nunca fazer `git add -f` nele. O arquivo `training.sqlite3` (e `-wal`/`-shm`) nunca é versionado.
- As 8.232 moléculas `OK` do protocolo 1 **nunca** são recalculadas nem perdem valor. Elas só ganham `protocol = "1"`.
- Só o conformer selecionado final é guardado (SQLite). Artefatos MOPAC/CREST de uma molécula bem-sucedida são apagados; os de uma falha são mantidos para diagnóstico.
- Nada de `git push`. Commits locais na branch `feat/delta-entalpia-treinamento`, uma preocupação por commit, mensagem em português no estilo do repositório, terminando com a linha `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
- Não executar o lote completo (~21 mil moléculas) em nenhuma tarefa. Operações sobre o `data/` real (migração, uma molécula) só na Tarefa 12, depois de o usuário confirmar.
- O servidor/worker congelado (`src/grimperium/server`, `src/grimperium/worker`) não é usado nem alterado.
- Unidades: entalpias em kcal/mol; `H298_pm7` vem de `energy_hof_kcal_mol` do MOPAC (PM7).
- Valores de `protocol`: `"1"` (preliminar, protocolo `crest_pm7`) e `"2"` (canônico). Strings, nunca números.

## Review Focus

1. **CSV ainda não migrado** (sem coluna `protocol`, ou com linhas `OK` sem marca): o ledger recusa trabalhar, `POST /api/training/start` responde 409 com o comando de migração, `GET /api/training/status` responde normalmente com `migrated=false`.
2. **Processo morto no meio de uma molécula:** linhas ficam `Running` e o job fica `running` sem batimento. O job é reportado como `stale`, um novo `run` pode começar e as linhas `Running` voltam para a fila.
3. **Geometria que não é mínimo** (sela ou verificação falha): nunca vira `OK`/protocolo 2; vira `Rerun` e, ao atingir o limite, `Skip`, e o motivo aparece em `recent_failures`.
4. **Radicais e espécies carregadas** (653 com multiplicidade 2; 13 com carga): `charge` e `multiplicity` do CSV (lidos como `1.0`/`2.0`) chegam como inteiros ao `MolecularIdentity`.
5. **POST de outra origem** (página qualquer tentando disparar processo em `127.0.0.1:8000`): rejeitado com 403; `workers` e `limit` fora dos limites: 422.

## File Structure

Criar:
- `src/semi_imperium/training/__init__.py` — exporta a API pública do pacote.
- `src/semi_imperium/training/protocol.py` — definição do protocolo 2 (settings, configuração, assinatura).
- `src/semi_imperium/training/db.py` — conexão SQLite (WAL) e esquema das duas tabelas de apoio.
- `src/semi_imperium/training/geometry_store.py` — `GeometryRecord`, `GeometryStore`.
- `src/semi_imperium/training/job_state.py` — `JobState`, `JobSnapshot`, falhas recentes.
- `src/semi_imperium/training/migration.py` — `migrate()` do CSV.
- `src/semi_imperium/training/ledger.py` — `TrainingLedger`, `read_progress()`.
- `src/semi_imperium/training/worker.py` — `WorkerTask`, `WorkerResult`, `run_task()`.
- `src/semi_imperium/training/runner.py` — `RunSettings`, `run_training()`.
- `src/semi_imperium/training/__main__.py` — CLI (`status`, `migrate`, `run`, `stop`).
- `src/grimperium/webapi/training.py` — status/start/stop da API.
- `tests/training/__init__.py`, `tests/training/conftest.py`, e um arquivo de teste por módulo acima.
- `tests/webapi/test_training_api.py`.
- `apps/dashboard/src/TrainingView.tsx`, `ComingSoonView.tsx`, `DeltaEnthalpyModule.tsx` (renomeia `GrimperiumModule.tsx`), `TrainingView.test.tsx`.

Modificar:
- `src/grimperium/crest_pm7/batch/csv_manager.py` (coluna `protocol`), `tests/test_csv_schema.py`, `tests/unit/test_csv_schema.py`, `tests/test_pr6d_schema_reconciliation.py`.
- `src/semi_imperium/workflows/calculation.py` e `src/semi_imperium/workflows/execution.py` (`ExecutionOutcome.geometry`).
- `src/grimperium/webapi/app.py`, `models.py`, `config.py`; remover `catalog.py`.
- `apps/dashboard/src/{api.ts,useApi.ts,modules.ts,DeltaModule.tsx,ResultsView.tsx,SummaryView.tsx,App.tsx,App.test.tsx,style.css}`; remover `RoutinesView.tsx`, `routines.ts`, `stats.test.ts` (se só testa `routines.ts`).
- `.gitignore`, `CHANGELOG.md`, `apps/dashboard/README.md`.

---

### Task 0: Branch, linha de base e diagnóstico da suíte `webapi`

**Files:**
- Modify: `.gitignore`

- [ ] **Step 1: Criar a branch e versionar o design e o plano**

```bash
cd /home/plasma-test/Projetos/grimperium
git switch -c feat/delta-entalpia-treinamento
git add docs/superpowers/specs/2026-10-12-delta-entalpia-treinamento-design.md docs/superpowers/plans/2026-10-12-delta-entalpia-treinamento.md
git commit -m "docs: design e plano da rotina de treinamento Delta Entalpia" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 2: Ignorar o SQLite do treinamento**

Acrescentar ao fim de `.gitignore`:

```gitignore

# Delta Entalpia: estado e geometrias do treinamento (SQLite + WAL)
data/semi_imperium/*.sqlite3*
```

Verificar: `git check-ignore -v data/semi_imperium/training.sqlite3 data/semi_imperium/training.sqlite3-wal` deve imprimir as duas linhas.

- [ ] **Step 3: Rodar a linha de base com limite de tempo**

Na sessão anterior, `tests/webapi` parou no primeiro teste. Reproduzir antes de mudar qualquer coisa:

```bash
timeout 180 poetry run pytest tests/webapi -x -v 2>&1 | tail -30
timeout 300 poetry run pytest tests/test_csv_schema.py tests/unit/test_csv_schema.py tests/test_pr6d_schema_reconciliation.py tests/semi_imperium -q 2>&1 | tail -15
```

Se `tests/webapi` travar ou falhar, **parar** e diagnosticar com `superpowers:systematic-debugging` (suspeita inicial: `store_progress` lendo o store real ou o `TestClient`). Registrar o resultado exato (passou/travou/falhou, qual teste) no topo do `CHANGELOG` rascunho da Tarefa 12. Não seguir para a Tarefa 1 com a suíte `webapi` travada.

- [ ] **Step 4: Commit**

```bash
git add .gitignore
git commit -m "chore: ignora o SQLite do treinamento Delta Entalpia" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 1: Definição do protocolo 2

**Files:**
- Create: `src/semi_imperium/training/__init__.py`, `src/semi_imperium/training/protocol.py`
- Create: `tests/training/__init__.py`, `tests/training/test_protocol.py`

**Interfaces:**
- Produces:
  - `LEGACY_PROTOCOL: str = "1"`, `CANONICAL_PROTOCOL: str = "2"`, `HAMILTONIAN: str = "PM7"`
  - `canonical_settings(base: SemiImperiumSettings | None = None) -> SemiImperiumSettings`
  - `canonical_configuration(settings: SemiImperiumSettings) -> EffectiveConfiguration`
  - `canonical_signature(settings: SemiImperiumSettings) -> str`

- [ ] **Step 1: Escrever os testes**

`tests/training/__init__.py` vazio. `tests/training/test_protocol.py`:

```python
"""Protocol 2 is CONFPASS + folding filter + PM7 with minimum required."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from semi_imperium.domain import ConformerSelectionStrategy, VerificationPolicy
from semi_imperium.settings import RuntimeSettings, SemiImperiumSettings
from semi_imperium.training.protocol import (
    CANONICAL_PROTOCOL,
    HAMILTONIAN,
    LEGACY_PROTOCOL,
    canonical_configuration,
    canonical_settings,
    canonical_signature,
)


def test_protocol_labels_are_strings() -> None:
    assert (LEGACY_PROTOCOL, CANONICAL_PROTOCOL, HAMILTONIAN) == ("1", "2", "PM7")


def test_canonical_settings_select_with_confpass_and_folding_filter() -> None:
    settings = canonical_settings()
    selection = settings.conformer_selection
    assert selection.resolved_strategy is ConformerSelectionStrategy.CONFPASS_PRIORITIZATION
    assert selection.folding_filter.enabled is True
    assert settings.conformer_search.enabled is True
    assert settings.verification.policy is VerificationPolicy.REQUIRE_MINIMUM


def test_canonical_configuration_is_pm7_with_crest() -> None:
    configuration = canonical_configuration(canonical_settings())
    assert configuration.semiempirical.hamiltonian == "PM7"
    assert configuration.conformer_search.enabled is True


def test_signature_ignores_runtime_knobs() -> None:
    base = canonical_settings()
    other_runtime = replace(
        base,
        runtime=RuntimeSettings(crest_threads=9, work_dir=Path("/tmp/elsewhere")),
    )
    assert canonical_signature(base) == canonical_signature(other_runtime)


def test_signature_changes_with_the_selection_strategy() -> None:
    energy_top_n = replace(
        canonical_settings(),
        conformer_selection=replace(
            canonical_settings().conformer_selection,
            strategy=ConformerSelectionStrategy.CREST_ENERGY_TOP_N.value,
        ),
    )
    assert canonical_signature(energy_top_n) != canonical_signature(canonical_settings())


def test_base_settings_are_not_mutated() -> None:
    base = SemiImperiumSettings()
    canonical_settings(base)
    assert base.conformer_selection.resolved_strategy is ConformerSelectionStrategy.CREST_ENERGY_TOP_N
    assert base.verification.policy is VerificationPolicy.NONE
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/training/test_protocol.py -v`
Expected: FAIL (`ModuleNotFoundError: semi_imperium.training`).

- [ ] **Step 3: Implementar**

`src/semi_imperium/training/__init__.py`:

```python
"""Training routine of the Delta Entalpia module (Semi-Imperium protocol 2)."""
```

`src/semi_imperium/training/protocol.py`:

```python
"""Protocol 2: the canonical Semi-Imperium protocol behind the training data.

Protocol 1 is the legacy ``crest_pm7`` pipeline that produced the first 8,232
rows of ``thermo_pm7.csv``. Protocol 2 selects conformers with CONFPASS plus
the folding filter (specification 2026-10-03) and demands a frequency check
that confirms a minimum. Both choices are part of the calculation signature,
so a number computed under another configuration is never mistaken for it.
"""

from __future__ import annotations

from dataclasses import replace

from semi_imperium.domain import (
    ConformerSelectionStrategy,
    EffectiveConfiguration,
    VerificationPolicy,
)
from semi_imperium.settings import SemiImperiumSettings

LEGACY_PROTOCOL = "1"
CANONICAL_PROTOCOL = "2"
HAMILTONIAN = "PM7"


def canonical_settings(
    base: SemiImperiumSettings | None = None,
) -> SemiImperiumSettings:
    """Return ``base`` (or the defaults) adjusted to protocol 2."""
    settings = base or SemiImperiumSettings()
    selection = settings.conformer_selection
    return replace(
        settings,
        conformer_search=replace(settings.conformer_search, enabled=True),
        conformer_selection=replace(
            selection,
            strategy=ConformerSelectionStrategy.CONFPASS_PRIORITIZATION.value,
            folding_filter=replace(selection.folding_filter, enabled=True),
        ),
        verification=replace(
            settings.verification, policy=VerificationPolicy.REQUIRE_MINIMUM
        ),
    )


def canonical_configuration(settings: SemiImperiumSettings) -> EffectiveConfiguration:
    """Resolve the one configuration every protocol-2 calculation runs under."""
    return settings.configuration_for(HAMILTONIAN, crest_enabled=True)


def canonical_signature(settings: SemiImperiumSettings) -> str:
    """Digest that identifies protocol 2 under ``settings``."""
    return canonical_configuration(settings).signature().digest
```

- [ ] **Step 4: Rodar e ver passar**

Run: `poetry run pytest tests/training/test_protocol.py -v`
Expected: PASS (6 testes).

- [ ] **Step 5: Commit**

```bash
git add src/semi_imperium/training tests/training
git commit -m "feat(training): define o protocolo 2 (CONFPASS + filtro, PM7, mínimo exigido)" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Geometria final no resultado da execução

**Files:**
- Modify: `src/semi_imperium/workflows/calculation.py` (dataclass `ExecutionOutcome`)
- Modify: `src/semi_imperium/workflows/execution.py` (`outcome_from_hamiltonian_result`)
- Test: `tests/integration/test_semi_imperium_end_to_end.py::test_conformer_boundaries_and_hamiltonian_surfaces_remain_independent` (já projeta AM1/PM3/PM7 com `outcome_from_hamiltonian_result`, no dicionário `projected`)

**Interfaces:**
- Produces: `ExecutionOutcome.geometry: ConformerGeometry | None = None` (campo novo, último, com padrão `None`, para não quebrar chamadas existentes). Em `outcome_from_hamiltonian_result`, `geometry = chosen.optimized_geometry` nos ramos que devolvem `result=data` (VERIFIED, SADDLE, UNVERIFIED).

- [ ] **Step 1: Escrever a asserção (no teste de integração que já projeta os três Hamiltonianos)**

Em `tests/integration/test_semi_imperium_end_to_end.py`, logo depois da asserção que compara `{name: item.state ...}` com `{"AM1": VERIFIED, "PM3": UNVERIFIED, "PM7": SADDLE}` (linha ~553), acrescentar:

```python
    for name in ("AM1", "PM3", "PM7"):
        hamiltonian_result = result.minima.for_hamiltonian(name)
        chosen = {a.attempt_id: a for a in hamiltonian_result.attempts}[
            hamiltonian_result.verified_attempt_id
            or hamiltonian_result.provisional_lowest_attempt_id
        ]
        assert projected[name].geometry is not None
        assert projected[name].geometry == chosen.optimized_geometry
```

Se o `MixedMopacBackend` do arquivo não preencher `provisional_lowest_attempt_id` para PM3/PM7, restringir o laço aos nomes em que o id existe (a asserção de estado logo acima já prova quais produzem resultado). Para o ramo sem geometria, acrescentar em `tests/unit/semi_imperium/test_calculate_workflow.py` (onde `ExecutionOutcome(` já é construído) uma asserção simples de que `ExecutionOutcome(state=CalculationState.FAILED).geometry is None`.

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/integration/test_semi_imperium_end_to_end.py tests/unit/semi_imperium/test_calculate_workflow.py -v -k "boundaries or outcome"`
Expected: FAIL (`AttributeError: 'ExecutionOutcome' object has no attribute 'geometry'`).

- [ ] **Step 3: Implementar**

Em `workflows/calculation.py`, dentro de `if TYPE_CHECKING:` (criar o bloco se não existir, junto dos imports) importar `from semi_imperium.conformers import ConformerGeometry` e, em `ExecutionOutcome`, acrescentar depois de `error_message`:

```python
    geometry: ConformerGeometry | None = None
    """Final optimized geometry of the chosen conformer; ``None`` on failure."""
```

(O módulo já usa `from __future__ import annotations`; se não usar, adicioná-lo. A importação fica só sob `TYPE_CHECKING` para evitar ciclo com `semi_imperium.conformers`.)

Em `workflows/execution.py`, nos três `return ExecutionOutcome(... result=data ...)` de `outcome_from_hamiltonian_result`, acrescentar `geometry=chosen.optimized_geometry,`.

- [ ] **Step 4: Rodar a suíte do pacote**

Run: `poetry run pytest tests/semi_imperium -q`
Expected: PASS (nenhum teste existente quebra; o campo novo tem padrão).

- [ ] **Step 5: Commit**

```bash
git add src/semi_imperium/workflows tests
git commit -m "feat(semi-imperium): ExecutionOutcome carrega a geometria final do conformer" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Tabela SQLite de geometrias

**Files:**
- Create: `src/semi_imperium/training/db.py`, `src/semi_imperium/training/geometry_store.py`
- Test: `tests/training/test_geometry_store.py`

**Interfaces:**
- Produces:
  - `db.connect(path: Path) -> sqlite3.Connection` — cria diretório/arquivo e as tabelas `final_geometry`, `training_job`, `training_failure`; WAL; `busy_timeout` 30 s; autocommit (`isolation_level=None`).
  - `GeometryRecord(mol_id, smiles, hamiltonian, protocol, signature, calculation_id, conformer_index, elements, coordinates, hof_kcal_mol, lowest_frequency_cm1, minimum_confirmed, created_at="")` com `.to_xyz() -> str`.
  - `GeometryStore(path)` com `save(record)`, `get(mol_id) -> GeometryRecord | None`, `count() -> int`, `mol_ids() -> set[str]`. Leituras num arquivo inexistente devolvem vazio **sem criar o arquivo**.

- [ ] **Step 1: Escrever os testes**

`tests/training/test_geometry_store.py`:

```python
from __future__ import annotations

from pathlib import Path

from semi_imperium.training.geometry_store import GeometryRecord, GeometryStore


def _record(mol_id: str = "mol_00001", hof: float = -12.5) -> GeometryRecord:
    return GeometryRecord(
        mol_id=mol_id,
        smiles="CO",
        hamiltonian="PM7",
        protocol="2",
        signature="abc123",
        calculation_id="calc-1",
        conformer_index=3,
        elements=("C", "O", "H"),
        coordinates=((0.0, 0.0, 0.0), (1.4, 0.0, 0.0), (-0.5, 0.9, 0.0)),
        hof_kcal_mol=hof,
        lowest_frequency_cm1=212.4,
        minimum_confirmed=True,
    )


def test_reading_a_missing_database_does_not_create_it(tmp_path: Path) -> None:
    path = tmp_path / "training.sqlite3"
    store = GeometryStore(path)
    assert store.count() == 0
    assert store.get("mol_00001") is None
    assert store.mol_ids() == set()
    assert not path.exists()


def test_save_then_get_round_trips(tmp_path: Path) -> None:
    store = GeometryStore(tmp_path / "training.sqlite3")
    store.save(_record())
    loaded = store.get("mol_00001")
    assert loaded is not None
    assert loaded.elements == ("C", "O", "H")
    assert loaded.coordinates[1] == (1.4, 0.0, 0.0)
    assert loaded.hof_kcal_mol == -12.5
    assert loaded.minimum_confirmed is True
    assert loaded.created_at  # filled on save


def test_saving_the_same_molecule_replaces_the_row(tmp_path: Path) -> None:
    store = GeometryStore(tmp_path / "training.sqlite3")
    store.save(_record(hof=-1.0))
    store.save(_record(hof=-2.0))
    assert store.count() == 1
    loaded = store.get("mol_00001")
    assert loaded is not None and loaded.hof_kcal_mol == -2.0


def test_xyz_export_keeps_atom_order_and_count(tmp_path: Path) -> None:
    lines = _record().to_xyz().splitlines()
    assert lines[0] == "3"
    assert "mol_00001" in lines[1]
    assert [line.split()[0] for line in lines[2:]] == ["C", "O", "H"]


def test_mol_ids_lists_every_saved_molecule(tmp_path: Path) -> None:
    store = GeometryStore(tmp_path / "training.sqlite3")
    store.save(_record("mol_00001"))
    store.save(_record("mol_00002"))
    assert store.mol_ids() == {"mol_00001", "mol_00002"}
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/training/test_geometry_store.py -v`
Expected: FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Implementar `db.py`**

```python
"""SQLite plumbing shared by the training stores (one file, WAL mode).

Several processes touch this file at once (the coordinator, the dashboard API
and the CLI), so every call opens a short-lived connection instead of sharing
one, and the schema is created idempotently on connect.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

DB_FILENAME = "training.sqlite3"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS final_geometry (
    mol_id TEXT PRIMARY KEY,
    smiles TEXT NOT NULL,
    hamiltonian TEXT NOT NULL,
    protocol TEXT NOT NULL,
    signature TEXT NOT NULL,
    calculation_id TEXT NOT NULL,
    conformer_index INTEGER,
    n_atoms INTEGER NOT NULL,
    elements TEXT NOT NULL,
    coordinates TEXT NOT NULL,
    hof_kcal_mol REAL,
    lowest_frequency_cm1 REAL,
    minimum_confirmed INTEGER NOT NULL,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS training_job (
    id INTEGER PRIMARY KEY CHECK (id = 1),
    state TEXT NOT NULL,
    pid INTEGER,
    started_at TEXT,
    heartbeat_at TEXT,
    finished_at TEXT,
    workers INTEGER NOT NULL DEFAULT 0,
    limit_count INTEGER,
    done INTEGER NOT NULL DEFAULT 0,
    failed INTEGER NOT NULL DEFAULT 0,
    message TEXT NOT NULL DEFAULT '',
    stop_requested INTEGER NOT NULL DEFAULT 0
);
CREATE TABLE IF NOT EXISTS training_failure (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mol_id TEXT NOT NULL,
    reason TEXT NOT NULL,
    at TEXT NOT NULL
);
"""


def database_path(store_root: Path) -> Path:
    """Location of the training database inside a Semi-Imperium store."""
    return store_root / DB_FILENAME


def connect(path: Path) -> sqlite3.Connection:
    """Open ``path`` (creating it and its schema); autocommit, WAL, 30 s wait."""
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path, timeout=30.0, isolation_level=None)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA busy_timeout=30000")
    connection.executescript(_SCHEMA)
    return connection
```

- [ ] **Step 4: Implementar `geometry_store.py`**

```python
"""Final geometry of every protocol-2 molecule, one row per ``mol_id``.

Only the optimized geometry of the selected conformer is kept; the CSV never
holds coordinates. ``to_xyz`` exports one molecule for other programs.
"""

from __future__ import annotations

import json
from contextlib import closing
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from semi_imperium.training.db import connect


@dataclass(frozen=True)
class GeometryRecord:
    """One molecule's final geometry plus the numbers that justify it."""

    mol_id: str
    smiles: str
    hamiltonian: str
    protocol: str
    signature: str
    calculation_id: str
    conformer_index: int | None
    elements: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    hof_kcal_mol: float | None
    lowest_frequency_cm1: float | None
    minimum_confirmed: bool
    created_at: str = ""

    def to_xyz(self) -> str:
        """Render the geometry as one XYZ block, preserving atom order."""
        comment = (
            f"{self.mol_id} {self.hamiltonian} protocol={self.protocol} "
            f"hof_kcal_mol={self.hof_kcal_mol}"
        )
        lines = [str(len(self.elements)), comment]
        for element, (x, y, z) in zip(self.elements, self.coordinates, strict=True):
            lines.append(f"{element:<2} {x:16.10f} {y:16.10f} {z:16.10f}")
        return "\n".join(lines) + "\n"


class GeometryStore:
    """SQLite-backed table of final geometries (see :mod:`...training.db`)."""

    def __init__(self, path: Path | str) -> None:
        self._path = Path(path)

    def save(self, record: GeometryRecord) -> None:
        """Insert or replace the row of ``record.mol_id``."""
        created_at = record.created_at or datetime.now(timezone.utc).isoformat()
        with closing(connect(self._path)) as connection:
            connection.execute(
                "INSERT OR REPLACE INTO final_geometry (mol_id, smiles, hamiltonian,"
                " protocol, signature, calculation_id, conformer_index, n_atoms,"
                " elements, coordinates, hof_kcal_mol, lowest_frequency_cm1,"
                " minimum_confirmed, created_at)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    record.mol_id,
                    record.smiles,
                    record.hamiltonian,
                    record.protocol,
                    record.signature,
                    record.calculation_id,
                    record.conformer_index,
                    len(record.elements),
                    json.dumps(list(record.elements)),
                    json.dumps([list(point) for point in record.coordinates]),
                    record.hof_kcal_mol,
                    record.lowest_frequency_cm1,
                    int(record.minimum_confirmed),
                    created_at,
                ),
            )

    def get(self, mol_id: str) -> GeometryRecord | None:
        """Return the stored geometry of ``mol_id``, or ``None``."""
        if not self._path.exists():
            return None
        with closing(connect(self._path)) as connection:
            row = connection.execute(
                "SELECT * FROM final_geometry WHERE mol_id = ?", (mol_id,)
            ).fetchone()
        if row is None:
            return None
        return GeometryRecord(
            mol_id=row["mol_id"],
            smiles=row["smiles"],
            hamiltonian=row["hamiltonian"],
            protocol=row["protocol"],
            signature=row["signature"],
            calculation_id=row["calculation_id"],
            conformer_index=row["conformer_index"],
            elements=tuple(json.loads(row["elements"])),
            coordinates=tuple(
                (float(p[0]), float(p[1]), float(p[2]))
                for p in json.loads(row["coordinates"])
            ),
            hof_kcal_mol=row["hof_kcal_mol"],
            lowest_frequency_cm1=row["lowest_frequency_cm1"],
            minimum_confirmed=bool(row["minimum_confirmed"]),
            created_at=row["created_at"],
        )

    def count(self) -> int:
        """Number of stored geometries (0 when the database does not exist)."""
        if not self._path.exists():
            return 0
        with closing(connect(self._path)) as connection:
            return int(
                connection.execute("SELECT COUNT(*) FROM final_geometry").fetchone()[0]
            )

    def mol_ids(self) -> set[str]:
        """Every ``mol_id`` that has a stored geometry."""
        if not self._path.exists():
            return set()
        with closing(connect(self._path)) as connection:
            return {
                str(row["mol_id"])
                for row in connection.execute("SELECT mol_id FROM final_geometry")
            }
```

- [ ] **Step 5: Rodar e ver passar**

Run: `poetry run pytest tests/training/test_geometry_store.py -v`
Expected: PASS (5 testes).

- [ ] **Step 6: Commit**

```bash
git add src/semi_imperium/training tests/training
git commit -m "feat(training): tabela SQLite com a geometria final de cada molécula" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Coluna `protocol` no esquema do CSV

**Files:**
- Modify: `src/grimperium/crest_pm7/batch/csv_manager.py` (grupos de colunas, `get_schema`, `_CSV_DTYPE`, docstrings "61")
- Modify: `tests/test_csv_schema.py`, `tests/unit/test_csv_schema.py`, `tests/test_pr6d_schema_reconciliation.py`

**Interfaces:**
- Produces: `BatchCSVManager.PROTOCOL_COLUMNS = ["protocol"]`; `get_schema()` devolve 62 colunas com `"protocol"` por último; `load_csv` lê `protocol` como `str`.

- [ ] **Step 1: Atualizar os testes (falham antes da mudança)**

Em `tests/test_csv_schema.py`: acrescentar `"protocol",` ao fim de `EXPECTED_COLUMNS` (depois de `"delta_correction"`), com o comentário `# Protocol provenance (62)`, e trocar a docstring/contagens "61" por "62". Em `tests/unit/test_csv_schema.py` linha ~37: `assert len(schema) == 62, f"Expected 62 columns, got {len(schema)}"`, e a mesma troca onde houver lista literal. Em `tests/test_pr6d_schema_reconciliation.py`: trocar `== 61` por `== 62` (linhas ~31 e ~65) e renomear os dois testes para `..._with_62_columns` / `..._62_column_contract`.

Acrescentar em `tests/test_csv_schema.py` um teste de leitura em texto:

```python
def test_protocol_column_round_trips_as_text(tmp_path: Path) -> None:
    from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager

    csv_path = tmp_path / "thermo_pm7.csv"
    manager = BatchCSVManager(csv_path)
    manager.create_empty_csv()
    frame = pd.DataFrame(
        [
            {"mol_id": "a", "status": "OK", "smiles": "C", "nheavy": 1, "protocol": "1"},
            {"mol_id": "b", "status": "Pending", "smiles": "CC", "nheavy": 2, "protocol": ""},
        ]
    )
    frame.to_csv(csv_path, index=False)
    loaded = BatchCSVManager(csv_path).load_csv()
    assert loaded["protocol"].fillna("").tolist() == ["1", ""]
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/test_csv_schema.py tests/unit/test_csv_schema.py tests/test_pr6d_schema_reconciliation.py -q`
Expected: FAIL (esquema ainda com 61 colunas; `protocol` lido como float).

- [ ] **Step 3: Implementar**

Em `csv_manager.py`: após `PREDICTION_OUTPUT_COLUMNS`:

```python
    PROTOCOL_COLUMNS = ["protocol"]
    """Which protocol produced the row: "1" legacy crest_pm7, "2" canonical."""
```

Em `get_schema()` acrescentar `+ self.PROTOCOL_COLUMNS` ao fim da soma e atualizar docstring e o comentário `# ── Column group definitions (62 columns total) ──`. Em `_CSV_DTYPE` acrescentar `"protocol": str,`.

- [ ] **Step 4: Rodar e corrigir o que mais enumerar colunas**

Run: `poetry run pytest tests/test_csv_schema.py tests/unit tests/test_pr6d_schema_reconciliation.py tests/test_batch_output_contract.py tests/test_csv_manager_retry.py tests/worker -q`
Expected: PASS. Se algum teste exigir que toda coluna do esquema pertença a `BATCH_STATE_COLUMNS` ou `SCIENTIFIC_RESULT_COLUMNS` (`output_contracts.py`), acrescentar `"protocol"` a `SCIENTIFIC_RESULT_COLUMNS` e repetir.

- [ ] **Step 5: Commit**

```bash
git add src/grimperium/crest_pm7/batch tests
git commit -m "feat(csv): coluna protocol no esquema do thermo_pm7 (62 colunas)" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Migração do CSV (protocolo 1 e reset de `Selected`)

**Files:**
- Create: `src/semi_imperium/training/migration.py`
- Test: `tests/training/test_migration.py`, `tests/training/conftest.py`

**Interfaces:**
- Consumes: `BatchCSVManager` (Tarefa 4), `LEGACY_PROTOCOL` (Tarefa 1).
- Produces: `PROTOCOL_COLUMN = "protocol"`; `MigrationReport(rows, legacy_marked, selected_reset, already_marked, applied, backup_path)`; `migrate(csv_path: Path, *, apply: bool = False) -> MigrationReport` (padrão é ensaio: não escreve nada).

- [ ] **Step 1: Fixture compartilhada e testes**

`tests/training/conftest.py`:

```python
"""Shared builders for the training tests."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager

CsvFactory = Callable[[list[dict[str, Any]]], Path]


@pytest.fixture
def make_csv(tmp_path: Path) -> CsvFactory:
    """Write a ``thermo_pm7.csv``-shaped file with the given rows."""

    def _make(rows: list[dict[str, Any]]) -> Path:
        csv_path = tmp_path / "thermo_pm7.csv"
        schema = BatchCSVManager(csv_path=None).get_schema()
        frame = pd.DataFrame(rows).reindex(columns=schema)
        frame.to_csv(csv_path, index=False)
        return csv_path

    return _make


def row(mol_id: str, status: str, **extra: Any) -> dict[str, Any]:
    """One CSV row with the identity columns the routine needs."""
    base: dict[str, Any] = {
        "mol_id": mol_id,
        "status": status,
        "smiles": "C",
        "multiplicity": 1,
        "charge": 0,
        "nheavy": 1,
        "H298_cbs": -17.8,
        "reruns": 0,
    }
    base.update(extra)
    return base
```

`tests/training/test_migration.py`:

```python
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from semi_imperium.training.migration import migrate
from tests.training.conftest import CsvFactory, row


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(make_csv: CsvFactory) -> Path:
    return make_csv(
        [
            row("a", "OK", H298_pm7=-3.83),
            row("b", "OK", H298_pm7=-5.5),
            row("c", "Selected"),
            row("d", "Pending"),
            row("e", "Skip"),
        ]
    )


def test_dry_run_reports_and_writes_nothing(make_csv: CsvFactory) -> None:
    path = _rows(make_csv)
    before = _digest(path)
    report = migrate(path)
    assert (report.rows, report.legacy_marked, report.selected_reset) == (5, 2, 1)
    assert report.applied is False and report.backup_path is None
    assert _digest(path) == before
    assert not list(path.parent.glob("*.pre-protocol-*"))


def test_apply_marks_legacy_rows_resets_selected_and_backs_up(
    make_csv: CsvFactory,
) -> None:
    path = _rows(make_csv)
    original = pd.read_csv(path, dtype=str, keep_default_na=False)
    report = migrate(path, apply=True)
    assert report.applied is True
    assert report.backup_path is not None and report.backup_path.is_file()
    after = pd.read_csv(path, dtype=str, keep_default_na=False)
    by_id = after.set_index("mol_id")
    assert by_id.loc["a", "protocol"] == "1" and by_id.loc["b", "protocol"] == "1"
    assert by_id.loc["c", "status"] == "Pending"
    assert by_id.loc["d", "protocol"] == "" and by_id.loc["e", "status"] == "Skip"
    # The preliminary values themselves are untouched.
    assert by_id.loc["a", "H298_pm7"] == original.set_index("mol_id").loc["a", "H298_pm7"]
    assert by_id.loc["b", "H298_pm7"] == original.set_index("mol_id").loc["b", "H298_pm7"]


def test_apply_is_idempotent(make_csv: CsvFactory) -> None:
    path = _rows(make_csv)
    migrate(path, apply=True)
    once = _digest(path)
    report = migrate(path, apply=True)
    assert (report.legacy_marked, report.selected_reset) == (0, 0)
    assert report.already_marked == 2
    assert _digest(path) == once


def test_rows_already_in_protocol_2_keep_their_mark(make_csv: CsvFactory) -> None:
    path = make_csv([row("a", "OK", protocol="2"), row("b", "OK")])
    migrate(path, apply=True)
    by_id = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("mol_id")
    assert by_id.loc["a", "protocol"] == "2"
    assert by_id.loc["b", "protocol"] == "1"


def test_missing_file_is_an_error_not_an_empty_csv(tmp_path: Path) -> None:
    import pytest

    with pytest.raises(FileNotFoundError):
        migrate(tmp_path / "absent.csv")
    assert not (tmp_path / "absent.csv").exists()
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/training/test_migration.py -v`
Expected: FAIL (`ModuleNotFoundError: semi_imperium.training.migration`).

- [ ] **Step 3: Implementar**

`src/semi_imperium/training/migration.py`:

```python
"""One-time CSV migration to the two-protocol layout.

* every ``OK`` row without a protocol is preliminary protocol-1 work and
  gets ``protocol = "1"`` (its values stay exactly as they are);
* every ``Selected`` row goes back to ``Pending`` so protocol 2 recomputes it.

The default is a dry run. Applying first copies the CSV next to itself, then
saves through ``BatchCSVManager`` (atomic write).
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from semi_imperium.training.protocol import LEGACY_PROTOCOL

PROTOCOL_COLUMN = "protocol"


@dataclass(frozen=True)
class MigrationReport:
    """What the migration found and, when applied, did."""

    rows: int
    legacy_marked: int
    selected_reset: int
    already_marked: int
    applied: bool
    backup_path: Path | None


def migrate(csv_path: Path, *, apply: bool = False) -> MigrationReport:
    """Report (and with ``apply=True`` perform) the migration of ``csv_path``."""
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    manager = BatchCSVManager(csv_path)
    frame = manager.load_csv()
    if PROTOCOL_COLUMN not in frame.columns:
        frame[PROTOCOL_COLUMN] = ""
    frame[PROTOCOL_COLUMN] = frame[PROTOCOL_COLUMN].fillna("").astype(str)

    ok = frame["status"] == MoleculeStatus.OK.value
    unmarked = ok & (frame[PROTOCOL_COLUMN] == "")
    selected = frame["status"] == MoleculeStatus.SELECTED.value
    report = MigrationReport(
        rows=len(frame),
        legacy_marked=int(unmarked.sum()),
        selected_reset=int(selected.sum()),
        already_marked=int((ok & (frame[PROTOCOL_COLUMN] != "")).sum()),
        applied=False,
        backup_path=None,
    )
    if not apply:
        return report
    if report.legacy_marked == 0 and report.selected_reset == 0:
        return MigrationReport(**{**report.__dict__, "applied": True})

    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    backup = csv_path.with_name(f"{csv_path.name}.pre-protocol-{stamp}.bak")
    shutil.copy2(csv_path, backup)
    frame.loc[unmarked, PROTOCOL_COLUMN] = LEGACY_PROTOCOL
    frame.loc[selected, "status"] = MoleculeStatus.PENDING.value
    manager.save_csv()
    return MigrationReport(**{**report.__dict__, "applied": True, "backup_path": backup})
```

- [ ] **Step 4: Rodar e ver passar**

Run: `poetry run pytest tests/training/test_migration.py -v`
Expected: PASS (5 testes). Se `test_apply_is_idempotent` falhar por diferença de bytes causada só pela reserialização do pandas na primeira escrita, comparar o `DataFrame` lido em vez do hash (a primeira migração pode reformatar floats; a segunda não escreve nada porque `legacy_marked == 0 and selected_reset == 0`).

- [ ] **Step 5: Commit**

```bash
git add src/semi_imperium/training tests/training
git commit -m "feat(training): migração do CSV (marca protocolo 1 e devolve Selected à fila)" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Ledger do CSV (fila, progresso e escrita de resultados)

**Files:**
- Create: `src/semi_imperium/training/ledger.py`
- Test: `tests/training/test_ledger.py`

**Interfaces:**
- Consumes: `BatchCSVManager`, `PROTOCOL_COLUMN`, `CANONICAL_PROTOCOL`.
- Produces:
  - `QueuedMolecule(mol_id, smiles, charge, multiplicity, nheavy, h298_cbs)`
  - `TrainingProgress(total, protocol2, protocol1, ok_unmarked, pending, selected, rerun, running, skipped, migrated)`
  - `read_progress(csv_path: Path) -> TrainingProgress` (leitura pura; `FileNotFoundError` se faltar o arquivo)
  - `NotMigratedError(RuntimeError)`
  - `TrainingLedger(csv_path, *, max_reruns=3)` com `load()`, `queue(limit: int | None = None) -> list[QueuedMolecule]`, `mark_running(mol_id)`, `record_success(mol_id, *, hof_kcal_mol, conformers_evaluated, elapsed_s, batch_id)`, `record_failure(mol_id, reason) -> str` (devolve o novo status).

- [ ] **Step 1: Escrever os testes**

`tests/training/test_ledger.py`:

```python
from __future__ import annotations

import pandas as pd
import pytest

from semi_imperium.training.ledger import (
    NotMigratedError,
    TrainingLedger,
    read_progress,
)
from tests.training.conftest import CsvFactory, row


def _migrated(make_csv: CsvFactory) -> object:
    return make_csv(
        [
            row("big", "Pending", nheavy=9, protocol=""),
            row("small", "Pending", nheavy=2, protocol=""),
            row("mid", "Rerun", nheavy=5, protocol="", reruns=1),
            row("orphan", "Running", nheavy=3, protocol=""),
            row("sel", "Selected", nheavy=4, protocol=""),
            row("old", "OK", nheavy=1, protocol="1", H298_pm7=-3.83),
            row("done", "OK", nheavy=1, protocol="2", H298_pm7=-1.0),
            row("skip", "Skip", nheavy=1, protocol=""),
        ]
    )


def test_progress_counts_each_bucket(make_csv: CsvFactory) -> None:
    progress = read_progress(_migrated(make_csv))  # type: ignore[arg-type]
    assert progress.total == 8
    assert (progress.protocol2, progress.protocol1, progress.ok_unmarked) == (1, 1, 0)
    assert (progress.pending, progress.selected, progress.rerun) == (2, 1, 1)
    assert (progress.running, progress.skipped) == (1, 1)
    assert progress.migrated is True


def test_progress_reports_unmigrated_csv(make_csv: CsvFactory) -> None:
    path = make_csv([row("a", "OK"), row("b", "Pending")])
    progress = read_progress(path)
    assert progress.migrated is False
    assert progress.ok_unmarked == 1


def test_ledger_refuses_an_unmigrated_csv(make_csv: CsvFactory) -> None:
    ledger = TrainingLedger(make_csv([row("a", "OK")]))
    with pytest.raises(NotMigratedError):
        ledger.load()


def test_queue_is_smallest_molecules_first_and_skips_finished_rows(
    make_csv: CsvFactory,
) -> None:
    ledger = TrainingLedger(_migrated(make_csv))  # type: ignore[arg-type]
    ledger.load()
    assert [m.mol_id for m in ledger.queue()] == ["small", "orphan", "sel", "mid", "big"]
    assert [m.mol_id for m in ledger.queue(limit=2)] == ["small", "orphan"]


def test_queue_casts_charge_and_multiplicity_to_int(make_csv: CsvFactory) -> None:
    path = make_csv(
        [row("rad", "Pending", multiplicity=2.0, charge=-1.0, protocol="")]
    )
    ledger = TrainingLedger(path)
    ledger.load()
    (molecule,) = ledger.queue()
    assert (molecule.charge, molecule.multiplicity) == (-1, 2)
    assert isinstance(molecule.charge, int) and isinstance(molecule.multiplicity, int)


def test_record_success_writes_value_delta_and_protocol(make_csv: CsvFactory) -> None:
    path = make_csv([row("a", "Pending", H298_cbs=-17.8, protocol="")])
    ledger = TrainingLedger(path)
    ledger.load()
    ledger.mark_running("a")
    ledger.record_success(
        "a", hof_kcal_mol=-16.5, conformers_evaluated=3, elapsed_s=120.0, batch_id="p2-x"
    )
    saved = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["a", "status"] == "OK"
    assert saved.loc["a", "protocol"] == "2"
    assert float(saved.loc["a", "H298_pm7"]) == pytest.approx(-16.5)
    assert float(saved.loc["a", "target_delta_kcalmol"]) == pytest.approx(-1.3)
    assert float(saved.loc["a", "total_time"]) == pytest.approx(2.0)
    assert saved.loc["a", "crest_status"] == "SUCCESS"
    assert saved.loc["a", "mopac_status"] == "OK"
    assert saved.loc["a", "batch_id"] == "p2-x"


def test_failures_become_rerun_until_the_limit_then_skip(make_csv: CsvFactory) -> None:
    path = make_csv([row("a", "Pending", protocol="")])
    ledger = TrainingLedger(path, max_reruns=2)
    ledger.load()
    assert ledger.record_failure("a", "saddle point") == "Rerun"
    assert ledger.record_failure("a", "saddle point") == "Skip"
    saved = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["a", "status"] == "Skip"
    assert saved.loc["a", "reruns"] in ("2", "2.0")
    assert saved.loc["a", "protocol"] == ""


def test_legacy_values_are_never_touched_by_a_run(make_csv: CsvFactory) -> None:
    path = make_csv(
        [row("old", "OK", protocol="1", H298_pm7=-3.83), row("new", "Pending", protocol="")]
    )
    ledger = TrainingLedger(path)
    ledger.load()
    ledger.mark_running("new")
    ledger.record_success(
        "new", hof_kcal_mol=-1.0, conformers_evaluated=1, elapsed_s=6.0, batch_id="p2"
    )
    saved = pd.read_csv(path, dtype=str, keep_default_na=False).set_index("mol_id")
    assert (saved.loc["old", "status"], saved.loc["old", "protocol"]) == ("OK", "1")
    assert float(saved.loc["old", "H298_pm7"]) == pytest.approx(-3.83)
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/training/test_ledger.py -v`
Expected: FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Implementar**

`src/semi_imperium/training/ledger.py`:

```python
"""The CSV side of the training routine.

``thermo_pm7.csv`` stays the dataset (the ML loaders read it), so the routine
writes its results there, through ``BatchCSVManager`` and from a single
process: the coordinator. Workers never touch the file.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from grimperium.crest_pm7.batch.enums import MoleculeStatus
from semi_imperium.training.migration import PROTOCOL_COLUMN
from semi_imperium.training.protocol import CANONICAL_PROTOCOL

#: Rows the routine still has to compute. ``Running`` is an orphan of an
#: interrupted run (the coordinator waits for in-flight work before exiting).
QUEUE_STATUSES = (
    MoleculeStatus.PENDING.value,
    MoleculeStatus.SELECTED.value,
    MoleculeStatus.RERUN.value,
    MoleculeStatus.RUNNING.value,
)


class NotMigratedError(RuntimeError):
    """The CSV has no ``protocol`` column or has unmarked ``OK`` rows."""


@dataclass(frozen=True)
class QueuedMolecule:
    """One molecule waiting for protocol 2."""

    mol_id: str
    smiles: str
    charge: int
    multiplicity: int
    nheavy: int
    h298_cbs: float | None


@dataclass(frozen=True)
class TrainingProgress:
    """Row counts per stage, with the denominator (``total``) explicit."""

    total: int
    protocol2: int
    protocol1: int
    ok_unmarked: int
    pending: int
    selected: int
    rerun: int
    running: int
    skipped: int
    migrated: bool


def read_progress(csv_path: Path) -> TrainingProgress:
    """Count rows by stage without loading the full CSV or touching it."""
    if not csv_path.is_file():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    header = pd.read_csv(csv_path, nrows=0).columns
    has_protocol = PROTOCOL_COLUMN in header
    columns = ["status", *([PROTOCOL_COLUMN] if has_protocol else [])]
    frame = pd.read_csv(csv_path, usecols=columns, dtype=str, keep_default_na=False)
    status = frame["status"].str.strip()
    if has_protocol:
        protocol = frame[PROTOCOL_COLUMN].str.strip().str.replace(r"\.0$", "", regex=True)
    else:
        protocol = pd.Series("", index=frame.index)
    ok = status == MoleculeStatus.OK.value
    unmarked = int((ok & ~protocol.isin(["1", "2"])).sum())
    return TrainingProgress(
        total=len(frame),
        protocol2=int((ok & (protocol == "2")).sum()),
        protocol1=int((ok & (protocol == "1")).sum()),
        ok_unmarked=unmarked,
        pending=int((status == MoleculeStatus.PENDING.value).sum()),
        selected=int((status == MoleculeStatus.SELECTED.value).sum()),
        rerun=int((status == MoleculeStatus.RERUN.value).sum()),
        running=int((status == MoleculeStatus.RUNNING.value).sum()),
        skipped=int((status == MoleculeStatus.SKIP.value).sum()),
        migrated=has_protocol and unmarked == 0,
    )


class TrainingLedger:
    """Queue, state transitions and result columns of the training CSV."""

    def __init__(self, csv_path: Path, *, max_reruns: int = 3) -> None:
        self._manager = BatchCSVManager(csv_path, max_reruns=max_reruns)

    def load(self) -> None:
        """Load the CSV and insist that it was migrated."""
        frame = self._manager.load_csv()
        if PROTOCOL_COLUMN not in frame.columns:
            raise NotMigratedError(
                "thermo_pm7.csv has no 'protocol' column; run "
                "'python -m semi_imperium.training migrate --apply' first"
            )
        frame[PROTOCOL_COLUMN] = frame[PROTOCOL_COLUMN].fillna("").astype(str)
        unmarked = (frame["status"] == MoleculeStatus.OK.value) & (
            frame[PROTOCOL_COLUMN] == ""
        )
        if unmarked.any():
            raise NotMigratedError(
                f"{int(unmarked.sum())} 'OK' rows have no protocol mark; run "
                "'python -m semi_imperium.training migrate --apply' first"
            )

    def _frame(self) -> pd.DataFrame:
        if self._manager.df is None:
            raise RuntimeError("TrainingLedger.load() must be called first")
        return self._manager.df

    def queue(self, limit: int | None = None) -> list[QueuedMolecule]:
        """Queued molecules, smallest first (ties by ``mol_id``)."""
        frame = self._frame()
        rows = frame.loc[frame["status"].isin(QUEUE_STATUSES)].copy()
        rows["_size"] = pd.to_numeric(rows["nheavy"], errors="coerce").fillna(1e9)
        rows = rows.sort_values(["_size", "mol_id"], kind="stable")
        queue = [
            QueuedMolecule(
                mol_id=str(item["mol_id"]),
                smiles=str(item["smiles"]),
                charge=_as_int(item.get("charge"), 0),
                multiplicity=_as_int(item.get("multiplicity"), 1),
                nheavy=_as_int(item.get("nheavy"), 0),
                h298_cbs=_as_float(item.get("H298_cbs")),
            )
            for _, item in rows.iterrows()
        ]
        return queue if limit is None else queue[:limit]

    def mark_running(self, mol_id: str) -> None:
        """Flag ``mol_id`` as in flight (persisted immediately)."""
        self._manager.apply_operational_status(mol_id, MoleculeStatus.RUNNING.value)

    def record_success(
        self,
        mol_id: str,
        *,
        hof_kcal_mol: float,
        conformers_evaluated: int,
        elapsed_s: float,
        batch_id: str,
    ) -> None:
        """Write a protocol-2 result: value, delta to CBS, provenance columns."""
        cbs = self._cell(mol_id, "H298_cbs")
        delta = None if cbs is None else cbs - hof_kcal_mol
        pct = (
            None
            if cbs is None or cbs == 0
            else round(abs(cbs - hof_kcal_mol) / abs(cbs) * 100, 2)
        )
        update: dict[str, Any] = {
            PROTOCOL_COLUMN: CANONICAL_PROTOCOL,
            "H298_pm7": hof_kcal_mol,
            "target_delta_kcalmol": delta,
            "abs_diff_%": pct,
            "crest_status": "SUCCESS",
            "mopac_status": "OK",
            "num_conformers_selected": conformers_evaluated,
            "total_time": round(elapsed_s / 60.0, 2),
            "timestamp": datetime.now().strftime("%m/%d-%H:%M"),
            "batch_id": batch_id,
        }
        self._manager.apply_operational_status(
            mol_id, MoleculeStatus.OK.value, result_update=update
        )

    def record_failure(self, mol_id: str, reason: str) -> str:
        """Count a failed attempt: ``Rerun`` below the limit, else ``Skip``."""
        del reason  # the reason is persisted by JobState.record_failure
        reruns = _as_int(self._cell(mol_id, "reruns"), 0) + 1
        status = (
            MoleculeStatus.SKIP.value
            if reruns >= self._manager.max_reruns
            else MoleculeStatus.RERUN.value
        )
        self._manager.apply_operational_status(mol_id, status, reruns=reruns)
        return status

    def _cell(self, mol_id: str, column: str) -> float | None:
        frame = self._frame()
        match = frame.loc[frame["mol_id"] == mol_id, column]
        return None if match.empty else _as_float(match.iloc[0])


def _as_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def _as_int(value: Any, default: int) -> int:
    number = _as_float(value)
    return default if number is None else int(number)
```

- [ ] **Step 4: Rodar e ver passar**

Run: `poetry run pytest tests/training/test_ledger.py -v`
Expected: PASS (8 testes). Se `apply_operational_status` recusar atribuir texto a uma coluna `float64` (por exemplo `timestamp`), `_CSV_DTYPE` já força `str` para ela; se aparecer o mesmo problema em `batch_id`/`crest_status`, confirmar que estão em `_CSV_DTYPE` (estão) e que o CSV de teste foi relido com `BatchCSVManager.load_csv`.

- [ ] **Step 5: Commit**

```bash
git add src/semi_imperium/training tests/training
git commit -m "feat(training): ledger do CSV (fila, progresso e resultados do protocolo 2)" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Estado do job e falhas recentes

**Files:**
- Create: `src/semi_imperium/training/job_state.py`
- Test: `tests/training/test_job_state.py`

**Interfaces:**
- Produces:
  - `STALE_AFTER_SECONDS = 180.0`
  - `JobAlreadyRunning(RuntimeError)`
  - `JobSnapshot(state, pid, started_at, heartbeat_at, finished_at, workers, limit, done, failed, message, stop_requested)`; `state ∈ {"idle","running","stopping","finished","failed","stale"}`; propriedade `is_active` (`running` ou `stopping`).
  - `JobState(path: Path)` com `begin(*, pid, workers, limit) -> None`, `heartbeat(*, done, failed) -> None`, `request_stop() -> bool` (False se não houver job ativo), `stop_requested() -> bool`, `mark_stopping() -> None`, `finish(state, message) -> None`, `snapshot() -> JobSnapshot`, `record_failure(mol_id, reason) -> None`, `recent_failures(limit=20) -> list[tuple[str, str, str]]` (`mol_id, reason, at`, mais recentes primeiro).
  - Leituras (`snapshot`, `recent_failures`, `stop_requested`) num arquivo inexistente não criam o arquivo.

- [ ] **Step 1: Escrever os testes**

```python
from __future__ import annotations

import os
from pathlib import Path

import pytest

from semi_imperium.training.job_state import (
    JobAlreadyRunning,
    JobState,
)


def _state(tmp_path: Path) -> JobState:
    return JobState(tmp_path / "training.sqlite3")


def test_missing_database_is_idle_and_is_not_created(tmp_path: Path) -> None:
    state = _state(tmp_path)
    snapshot = state.snapshot()
    assert snapshot.state == "idle" and not snapshot.is_active
    assert state.recent_failures() == []
    assert state.stop_requested() is False
    assert not (tmp_path / "training.sqlite3").exists()


def test_begin_heartbeat_finish_cycle(tmp_path: Path) -> None:
    state = _state(tmp_path)
    state.begin(pid=os.getpid(), workers=3, limit=10)
    snapshot = state.snapshot()
    assert snapshot.state == "running" and snapshot.is_active
    assert (snapshot.workers, snapshot.limit, snapshot.pid) == (3, 10, os.getpid())
    state.heartbeat(done=4, failed=1)
    assert (state.snapshot().done, state.snapshot().failed) == (4, 1)
    state.finish("finished", "ok")
    assert state.snapshot().state == "finished" and not state.snapshot().is_active


def test_second_begin_is_refused_while_alive(tmp_path: Path) -> None:
    state = _state(tmp_path)
    state.begin(pid=os.getpid(), workers=1, limit=None)
    with pytest.raises(JobAlreadyRunning):
        state.begin(pid=os.getpid(), workers=1, limit=None)


def test_dead_pid_is_reported_stale_and_can_be_replaced(tmp_path: Path) -> None:
    state = _state(tmp_path)
    state.begin(pid=2_000_000_000, workers=1, limit=None)  # no such process
    assert state.snapshot().state == "stale"
    assert state.snapshot().is_active is False
    state.begin(pid=os.getpid(), workers=2, limit=None)
    assert state.snapshot().state == "running"


def test_old_heartbeat_with_live_pid_is_stale(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import semi_imperium.training.job_state as module

    state = _state(tmp_path)
    state.begin(pid=os.getpid(), workers=1, limit=None)
    monkeypatch.setattr(module, "STALE_AFTER_SECONDS", -1.0)
    assert state.snapshot().state == "stale"


def test_request_stop_only_when_active(tmp_path: Path) -> None:
    state = _state(tmp_path)
    assert state.request_stop() is False
    state.begin(pid=os.getpid(), workers=1, limit=None)
    assert state.request_stop() is True
    assert state.stop_requested() is True
    state.mark_stopping()
    assert state.snapshot().state == "stopping"
    assert state.snapshot().is_active is True


def test_failures_are_kept_newest_first(tmp_path: Path) -> None:
    state = _state(tmp_path)
    state.record_failure("a", "saddle point")
    state.record_failure("b", "CREST timeout")
    failures = state.recent_failures(limit=1)
    assert [(mol, reason) for mol, reason, _ in failures] == [("b", "CREST timeout")]
    assert len(state.recent_failures()) == 2
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/training/test_job_state.py -v`
Expected: FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Implementar**

```python
"""Single-row job state in SQLite: who runs, how far, and whether to stop.

The coordinator writes a heartbeat; readers (dashboard API, CLI) decide whether
a job that claims to run is really alive by checking the pid and the
heartbeat age. A job that is not alive is reported as ``stale`` and does not
block a new one.
"""

from __future__ import annotations

import os
import sqlite3
from contextlib import closing
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from semi_imperium.training.db import connect

STALE_AFTER_SECONDS = 180.0
_ACTIVE = ("running", "stopping")


class JobAlreadyRunning(RuntimeError):
    """Another training job is alive."""


@dataclass(frozen=True)
class JobSnapshot:
    """What readers see of the job, with staleness already decided."""

    state: str
    pid: int | None
    started_at: str | None
    heartbeat_at: str | None
    finished_at: str | None
    workers: int
    limit: int | None
    done: int
    failed: int
    message: str
    stop_requested: bool

    @property
    def is_active(self) -> bool:
        """Whether a live process is running or winding down."""
        return self.state in _ACTIVE


_IDLE = JobSnapshot("idle", None, None, None, None, 0, None, 0, 0, "", False)


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _pid_alive(pid: int | None) -> bool:
    if pid is None:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OverflowError:
        return False
    return True


def _snapshot(row: sqlite3.Row) -> JobSnapshot:
    state = str(row["state"])
    if state in _ACTIVE:
        beat = row["heartbeat_at"]
        age = (
            (_now() - datetime.fromisoformat(beat)).total_seconds()
            if beat
            else float("inf")
        )
        if not _pid_alive(row["pid"]) or age > STALE_AFTER_SECONDS:
            state = "stale"
    return JobSnapshot(
        state=state,
        pid=row["pid"],
        started_at=row["started_at"],
        heartbeat_at=row["heartbeat_at"],
        finished_at=row["finished_at"],
        workers=int(row["workers"]),
        limit=row["limit_count"],
        done=int(row["done"]),
        failed=int(row["failed"]),
        message=str(row["message"]),
        stop_requested=bool(row["stop_requested"]),
    )


class JobState:
    """Reads and writes the ``training_job`` and ``training_failure`` tables."""

    def __init__(self, path: Path | str) -> None:
        self._path = Path(path)

    def snapshot(self) -> JobSnapshot:
        """Current job state (``idle`` if nothing ever ran)."""
        if not self._path.exists():
            return _IDLE
        with closing(connect(self._path)) as connection:
            row = connection.execute("SELECT * FROM training_job WHERE id = 1").fetchone()
        return _IDLE if row is None else _snapshot(row)

    def begin(self, *, pid: int, workers: int, limit: int | None) -> None:
        """Claim the job slot, atomically refusing when a live job exists."""
        now = _now().isoformat()
        with closing(connect(self._path)) as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                row = connection.execute(
                    "SELECT * FROM training_job WHERE id = 1"
                ).fetchone()
                if row is not None and _snapshot(row).is_active:
                    raise JobAlreadyRunning(
                        f"training job already running (pid {row['pid']})"
                    )
                connection.execute(
                    "INSERT OR REPLACE INTO training_job (id, state, pid, started_at,"
                    " heartbeat_at, finished_at, workers, limit_count, done, failed,"
                    " message, stop_requested)"
                    " VALUES (1, 'running', ?, ?, ?, NULL, ?, ?, 0, 0, '', 0)",
                    (pid, now, now, workers, limit),
                )
                connection.execute("COMMIT")
            except BaseException:
                connection.execute("ROLLBACK")
                raise

    def heartbeat(self, *, done: int, failed: int) -> None:
        """Record progress and prove the coordinator is alive."""
        with closing(connect(self._path)) as connection:
            connection.execute(
                "UPDATE training_job SET heartbeat_at = ?, done = ?, failed = ?"
                " WHERE id = 1",
                (_now().isoformat(), done, failed),
            )

    def request_stop(self) -> bool:
        """Ask the running job to stop dispatching; ``False`` if none is active."""
        if not self.snapshot().is_active:
            return False
        with closing(connect(self._path)) as connection:
            connection.execute("UPDATE training_job SET stop_requested = 1 WHERE id = 1")
        return True

    def stop_requested(self) -> bool:
        """Whether a stop was requested for the current job."""
        return self.snapshot().stop_requested

    def mark_stopping(self) -> None:
        """Show that the coordinator saw the stop and is draining in-flight work."""
        with closing(connect(self._path)) as connection:
            connection.execute(
                "UPDATE training_job SET state = 'stopping' WHERE id = 1"
                " AND state = 'running'"
            )

    def finish(self, state: str, message: str) -> None:
        """Close the job as ``finished`` or ``failed``."""
        with closing(connect(self._path)) as connection:
            connection.execute(
                "UPDATE training_job SET state = ?, message = ?, finished_at = ?"
                " WHERE id = 1",
                (state, message, _now().isoformat()),
            )

    def record_failure(self, mol_id: str, reason: str) -> None:
        """Keep why one molecule failed (the CSV only keeps the status)."""
        with closing(connect(self._path)) as connection:
            connection.execute(
                "INSERT INTO training_failure (mol_id, reason, at) VALUES (?, ?, ?)",
                (mol_id, reason, _now().isoformat()),
            )

    def recent_failures(self, limit: int = 20) -> list[tuple[str, str, str]]:
        """``(mol_id, reason, at)`` of the latest failures, newest first."""
        if not self._path.exists():
            return []
        with closing(connect(self._path)) as connection:
            rows = connection.execute(
                "SELECT mol_id, reason, at FROM training_failure"
                " ORDER BY id DESC LIMIT ?",
                (limit,),
            ).fetchall()
        return [(str(r["mol_id"]), str(r["reason"]), str(r["at"])) for r in rows]
```

- [ ] **Step 4: Rodar e ver passar**

Run: `poetry run pytest tests/training/test_job_state.py -v`
Expected: PASS (7 testes).

- [ ] **Step 5: Commit**

```bash
git add src/semi_imperium/training tests/training
git commit -m "feat(training): estado do job (batimento, parada, falhas recentes) em SQLite" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Worker (uma molécula) e coordenador

**Files:**
- Create: `src/semi_imperium/training/worker.py`, `src/semi_imperium/training/runner.py`
- Test: `tests/training/test_worker.py`, `tests/training/test_runner.py`

**Interfaces:**
- Consumes: tudo das Tarefas 1–7, `ScientificCalculationExecutor`, `SubprocessCrestRunner`, `ExecutionOutcome.geometry`.
- Produces (`worker.py`):
  - `WorkerTask(mol_id, smiles, charge, multiplicity, run_id, store_root, work_dir, crest_executable, xtb_executable, mopac_executable, crest_threads, crest_timeout_s, mopac_timeout_s, keep_artifacts)` — dataclass frozen, só tipos primitivos (picklável).
  - `WorkerResult(mol_id, ok, state, hof_kcal_mol, lowest_frequency_cm1, conformers_evaluated, conformer_index, elements, coordinates, calculation_id, signature, smiles, elapsed_s, error)`.
  - `run_task(task: WorkerTask, *, executor: CalculationExecutor | None = None) -> WorkerResult` — função de módulo (picklável); `executor` só é injetado em testes.
- Produces (`runner.py`):
  - `RunSettings(csv_path, store_root, work_dir, workers=3, limit=None, crest_executable="crest", xtb_executable="xtb", mopac_executable="mopac", crest_threads=4, crest_timeout_s=3600.0, mopac_timeout_s=900.0, max_reruns=3, keep_artifacts=False)`
  - `RunSummary(attempted, ok, failed, stopped)`
  - `run_training(settings, *, pool_factory=..., task_fn=run_task, install_signals=True) -> RunSummary`

- [ ] **Step 1: Testes do worker**

```python
from __future__ import annotations

from pathlib import Path

from semi_imperium.conformers import ConformerGeometry
from semi_imperium.domain import (
    CalculationResultData,
    CalculationState,
    VerificationOutcome,
)
from semi_imperium.training.worker import WorkerTask, run_task
from semi_imperium.workflows.calculation import ExecutionOutcome, ExecutionRequest

GEOMETRY = ConformerGeometry(
    elements=("C", "H", "H", "H", "H"),
    coordinates=(
        (0.0, 0.0, 0.0),
        (0.63, 0.63, 0.63),
        (-0.63, -0.63, 0.63),
        (-0.63, 0.63, -0.63),
        (0.63, -0.63, -0.63),
    ),
)


class FakeExecutor:
    def __init__(self, outcome: ExecutionOutcome) -> None:
        self.outcome = outcome
        self.requests: list[ExecutionRequest] = []

    def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
        self.requests.append(request)
        return self.outcome


def _task(tmp_path: Path, **overrides: object) -> WorkerTask:
    values: dict[str, object] = {
        "mol_id": "mol_00001",
        "smiles": "C",
        "charge": 0,
        "multiplicity": 1,
        "run_id": "run-test",
        "store_root": str(tmp_path / "store"),
        "work_dir": str(tmp_path / "work"),
        "crest_executable": "crest",
        "xtb_executable": "xtb",
        "mopac_executable": "mopac",
        "crest_threads": 1,
        "crest_timeout_s": 60.0,
        "mopac_timeout_s": 60.0,
        "keep_artifacts": False,
    }
    values.update(overrides)
    return WorkerTask(**values)  # type: ignore[arg-type]


def _verified() -> ExecutionOutcome:
    return ExecutionOutcome(
        state=CalculationState.VERIFIED,
        verification=VerificationOutcome.MINIMUM_CONFIRMED,
        result=CalculationResultData(
            energy_hof_kcal_mol=-17.9,
            conformer_index=2,
            conformers_evaluated=3,
            lowest_frequency_cm1=1306.0,
        ),
        geometry=GEOMETRY,
    )


def test_verified_minimum_is_ok_and_carries_the_geometry(tmp_path: Path) -> None:
    executor = FakeExecutor(_verified())
    result = run_task(_task(tmp_path), executor=executor)
    assert result.ok is True and result.state == "verified"
    assert result.hof_kcal_mol == -17.9 and result.conformers_evaluated == 3
    assert result.elements == GEOMETRY.elements
    assert result.coordinates == GEOMETRY.coordinates
    request = executor.requests[0]
    assert request.hamiltonian == "PM7"
    assert request.configuration.conformer_selection.folding_filter.enabled is True
    assert result.signature == request.configuration.signature().digest


def test_charge_and_multiplicity_reach_the_identity(tmp_path: Path) -> None:
    executor = FakeExecutor(_verified())
    run_task(_task(tmp_path, smiles="[CH3]", multiplicity=2), executor=executor)
    identity = executor.requests[0].identity
    assert (identity.charge, identity.multiplicity) == (0, 2)


def test_a_saddle_is_never_accepted(tmp_path: Path) -> None:
    outcome = ExecutionOutcome(
        state=CalculationState.SADDLE,
        verification=VerificationOutcome.SADDLE_POINT,
        result=CalculationResultData(energy_hof_kcal_mol=-10.0, lowest_frequency_cm1=-420.0),
        geometry=GEOMETRY,
    )
    result = run_task(_task(tmp_path), executor=FakeExecutor(outcome))
    assert result.ok is False and result.state == "saddle"
    assert "saddle" in (result.error or "")


def test_executor_exception_becomes_a_failed_result(tmp_path: Path) -> None:
    class Boom:
        def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
            raise RuntimeError("crest exploded")

    result = run_task(_task(tmp_path), executor=Boom())
    assert result.ok is False and result.state == "failed"
    assert "crest exploded" in (result.error or "")


def test_success_removes_artifacts_but_failure_keeps_them(tmp_path: Path) -> None:
    store = tmp_path / "store"

    class Leaves:
        def __init__(self, outcome: ExecutionOutcome) -> None:
            self.outcome = outcome

        def execute(self, request: ExecutionRequest) -> ExecutionOutcome:
            folder = store / "artifacts" / request.calculation_id
            folder.mkdir(parents=True)
            (folder / "trace.out").write_text("x")
            return self.outcome

    ok = run_task(_task(tmp_path), executor=Leaves(_verified()))
    assert not (store / "artifacts" / ok.calculation_id).exists()

    failed_outcome = ExecutionOutcome(
        state=CalculationState.FAILED,
        verification=VerificationOutcome.FAILED,
        error_message="no geometry",
    )
    bad = run_task(_task(tmp_path, run_id="run-2"), executor=Leaves(failed_outcome))
    assert (store / "artifacts" / bad.calculation_id / "trace.out").exists()


def test_the_calculation_record_is_persisted_in_the_store(tmp_path: Path) -> None:
    from semi_imperium.persistence import SemiImperiumStore

    result = run_task(_task(tmp_path), executor=FakeExecutor(_verified()))
    records = SemiImperiumStore(tmp_path / "store").list_calculations()
    assert [r.calculation_id for r in records] == [result.calculation_id]
    assert records[0].state is CalculationState.VERIFIED
    assert records[0].result is not None and records[0].result.artifact_paths == ()
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/training/test_worker.py -v`
Expected: FAIL (`ModuleNotFoundError: semi_imperium.training.worker`).

- [ ] **Step 3: Implementar `worker.py`**

```python
"""One protocol-2 molecule, run in a worker process.

``run_task`` is a plain module-level function over primitive arguments so a
``ProcessPoolExecutor`` can pickle it. The worker persists the Semi-Imperium
calculation record (the full audit trail) but never touches the CSV or the
geometry table: those belong to the coordinator, the only writer.
"""

from __future__ import annotations

import shutil
import time
from dataclasses import dataclass, replace
from pathlib import Path

from semi_imperium.domain import (
    CalculationRecord,
    CalculationState,
    MolecularIdentity,
    VerificationOutcome,
)
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.settings import RuntimeSettings, SemiImperiumSettings
from semi_imperium.training.protocol import (
    HAMILTONIAN,
    canonical_configuration,
    canonical_settings,
)
from semi_imperium.workflows.calculation import (
    CalculationExecutor,
    ExecutionOutcome,
    ExecutionRequest,
    build_provenance,
)
from semi_imperium.workflows.execution import ARTIFACTS_DIRNAME


@dataclass(frozen=True)
class WorkerTask:
    """Everything one worker needs, as picklable primitives."""

    mol_id: str
    smiles: str
    charge: int
    multiplicity: int
    run_id: str
    store_root: str
    work_dir: str
    crest_executable: str
    xtb_executable: str
    mopac_executable: str
    crest_threads: int
    crest_timeout_s: float
    mopac_timeout_s: float
    keep_artifacts: bool


@dataclass(frozen=True)
class WorkerResult:
    """What the coordinator needs to record one finished molecule."""

    mol_id: str
    smiles: str
    ok: bool
    state: str
    hof_kcal_mol: float | None
    lowest_frequency_cm1: float | None
    conformers_evaluated: int
    conformer_index: int | None
    elements: tuple[str, ...]
    coordinates: tuple[tuple[float, float, float], ...]
    calculation_id: str
    signature: str
    elapsed_s: float
    error: str | None


def _settings(task: WorkerTask) -> SemiImperiumSettings:
    runtime = RuntimeSettings(
        crest_executable=task.crest_executable,
        mopac_executable=task.mopac_executable,
        crest_threads=task.crest_threads,
        crest_timeout_seconds=task.crest_timeout_s,
        mopac_timeout_seconds=task.mopac_timeout_s,
        work_dir=Path(task.work_dir),
        store_root=Path(task.store_root),
    )
    return canonical_settings(SemiImperiumSettings(runtime=runtime))


def _default_executor(
    settings: SemiImperiumSettings, task: WorkerTask
) -> CalculationExecutor:
    from semi_imperium.conformers.crest_runner import SubprocessCrestRunner
    from semi_imperium.workflows.execution import ScientificCalculationExecutor

    runner = SubprocessCrestRunner(
        work_root=Path(task.work_dir) / "crest",
        executable=task.crest_executable,
        xtb_executable=task.xtb_executable,
        threads=task.crest_threads,
        timeout_seconds=task.crest_timeout_s,
    )
    return ScientificCalculationExecutor(settings, crest_runner=runner)


def run_task(
    task: WorkerTask, *, executor: CalculationExecutor | None = None
) -> WorkerResult:
    """Compute one molecule under protocol 2 and persist its audit record."""
    started = time.perf_counter()
    settings = _settings(task)
    configuration = canonical_configuration(settings)
    signature = configuration.signature()
    identity = MolecularIdentity.from_smiles(
        task.smiles, charge=task.charge, multiplicity=task.multiplicity
    )
    store = SemiImperiumStore(task.store_root)
    record = CalculationRecord(
        run_id=task.run_id,
        molecule=identity,
        signature=signature,
        provenance=build_provenance(settings),
    ).transition_to(CalculationState.RUNNING)
    store.save_calculation(record)

    runner = executor or _default_executor(settings, task)
    try:
        outcome = runner.execute(
            ExecutionRequest(
                run_id=task.run_id,
                calculation_id=record.calculation_id,
                identity=identity,
                configuration=configuration,
                hamiltonian=HAMILTONIAN,
            )
        )
    except Exception as exc:  # a crashing executor is a failed result
        outcome = ExecutionOutcome(
            state=CalculationState.FAILED,
            verification=VerificationOutcome.FAILED,
            error_message=f"{type(exc).__name__}: {exc}",
        )

    accepted = (
        outcome.state is CalculationState.VERIFIED
        and outcome.result is not None
        and outcome.result.energy_hof_kcal_mol is not None
        and outcome.geometry is not None
    )
    artifacts = Path(task.store_root) / ARTIFACTS_DIRNAME / record.calculation_id
    result_data = outcome.result
    if accepted and not task.keep_artifacts:
        shutil.rmtree(artifacts, ignore_errors=True)
        shutil.rmtree(Path(task.work_dir) / "crest" / identity.molecule_id, ignore_errors=True)
        if result_data is not None:
            result_data = replace(result_data, artifact_paths=())

    record = record.transition_to(
        outcome.state,
        verification=outcome.verification,
        result=result_data,
        error_message=outcome.error_message,
    )
    store.save_calculation(record)

    geometry = outcome.geometry
    error = None
    if not accepted:
        error = f"{outcome.state.value}: " + (
            outcome.error_message or outcome.verification.value
        )
    return WorkerResult(
        mol_id=task.mol_id,
        smiles=task.smiles,
        ok=accepted,
        state=outcome.state.value,
        hof_kcal_mol=None if result_data is None else result_data.energy_hof_kcal_mol,
        lowest_frequency_cm1=(
            None if result_data is None else result_data.lowest_frequency_cm1
        ),
        conformers_evaluated=0 if result_data is None else result_data.conformers_evaluated,
        conformer_index=None if result_data is None else result_data.conformer_index,
        elements=() if geometry is None else geometry.elements,
        coordinates=() if geometry is None else geometry.coordinates,
        calculation_id=record.calculation_id,
        signature=signature.digest,
        elapsed_s=time.perf_counter() - started,
        error=error,
    )
```

(Se `semi_imperium.domain` não exportar `CalculationRecord`/`MolecularIdentity`/`VerificationOutcome`, importar dos módulos `semi_imperium.domain.records` / `.identity` / `.enums`; `tests/webapi/test_app.py` já importa `CalculationRecord` e `MolecularIdentity` de `semi_imperium.domain`.)

- [ ] **Step 4: Rodar e ver passar**

Run: `poetry run pytest tests/training/test_worker.py -v`
Expected: PASS (6 testes).

- [ ] **Step 5: Testes do coordenador**

`tests/training/test_runner.py`:

```python
from __future__ import annotations

from concurrent.futures import Executor, ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import pytest

from semi_imperium.training.geometry_store import GeometryStore
from semi_imperium.training.job_state import JobAlreadyRunning, JobState
from semi_imperium.training.ledger import NotMigratedError
from semi_imperium.training.runner import RunSettings, run_training
from semi_imperium.training.worker import WorkerResult, WorkerTask
from tests.training.conftest import CsvFactory, row


def _threads(workers: int) -> Executor:
    return ThreadPoolExecutor(max_workers=workers)


def _fake(fail: set[str] | None = None):
    fail = fail or set()

    def task_fn(task: WorkerTask) -> WorkerResult:
        bad = task.mol_id in fail
        return WorkerResult(
            mol_id=task.mol_id,
            smiles=task.smiles,
            ok=not bad,
            state="saddle" if bad else "verified",
            hof_kcal_mol=None if bad else -1.5,
            lowest_frequency_cm1=None if bad else 100.0,
            conformers_evaluated=0 if bad else 2,
            conformer_index=None if bad else 1,
            elements=() if bad else ("C",),
            coordinates=() if bad else ((0.0, 0.0, 0.0),),
            calculation_id=f"calc-{task.mol_id}",
            signature="sig",
            elapsed_s=6.0,
            error="saddle: imaginary mode" if bad else None,
        )

    return task_fn


def _settings(tmp_path: Path, csv: Path, **overrides: object) -> RunSettings:
    values: dict[str, object] = {
        "csv_path": csv,
        "store_root": tmp_path / "store",
        "work_dir": tmp_path / "work",
        "workers": 2,
    }
    values.update(overrides)
    return RunSettings(**values)  # type: ignore[arg-type]


def _csv(make_csv: CsvFactory) -> Path:
    return make_csv(
        [
            row("a", "Pending", nheavy=1, protocol=""),
            row("b", "Pending", nheavy=2, protocol=""),
            row("c", "Pending", nheavy=3, protocol=""),
            row("old", "OK", nheavy=1, protocol="1", H298_pm7=-3.83),
        ]
    )


def test_run_records_results_in_csv_geometry_store_and_job(
    tmp_path: Path, make_csv: CsvFactory
) -> None:
    csv = _csv(make_csv)
    summary = run_training(
        _settings(tmp_path, csv), pool_factory=_threads, task_fn=_fake(), install_signals=False
    )
    assert (summary.attempted, summary.ok, summary.failed, summary.stopped) == (3, 3, 0, False)
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert [saved.loc[m, "protocol"] for m in ("a", "b", "c")] == ["2", "2", "2"]
    assert saved.loc["old", "protocol"] == "1"
    store = GeometryStore(tmp_path / "store" / "training.sqlite3")
    assert store.mol_ids() == {"a", "b", "c"}
    job = JobState(tmp_path / "store" / "training.sqlite3").snapshot()
    assert job.state == "finished" and job.done == 3 and job.failed == 0


def test_limit_caps_the_number_of_molecules(tmp_path: Path, make_csv: CsvFactory) -> None:
    csv = _csv(make_csv)
    summary = run_training(
        _settings(tmp_path, csv, limit=1),
        pool_factory=_threads,
        task_fn=_fake(),
        install_signals=False,
    )
    assert summary.attempted == 1
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["a", "status"] == "OK"
    assert saved.loc["b", "status"] == "Pending"


def test_failure_goes_to_rerun_and_its_reason_is_kept(
    tmp_path: Path, make_csv: CsvFactory
) -> None:
    csv = _csv(make_csv)
    summary = run_training(
        _settings(tmp_path, csv),
        pool_factory=_threads,
        task_fn=_fake(fail={"b"}),
        install_signals=False,
    )
    assert (summary.ok, summary.failed) == (2, 1)
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert (saved.loc["b", "status"], saved.loc["b", "protocol"]) == ("Rerun", "")
    store = GeometryStore(tmp_path / "store" / "training.sqlite3")
    assert "b" not in store.mol_ids()
    failures = JobState(tmp_path / "store" / "training.sqlite3").recent_failures()
    assert failures[0][0] == "b" and "saddle" in failures[0][1]


def test_orphan_running_rows_are_requeued(tmp_path: Path, make_csv: CsvFactory) -> None:
    csv = make_csv([row("orphan", "Running", nheavy=1, protocol="")])
    summary = run_training(
        _settings(tmp_path, csv), pool_factory=_threads, task_fn=_fake(), install_signals=False
    )
    assert summary.ok == 1


def test_stop_request_prevents_new_dispatch(tmp_path: Path, make_csv: CsvFactory) -> None:
    csv = _csv(make_csv)
    state = JobState(tmp_path / "store" / "training.sqlite3")

    def stopping_task(task: WorkerTask) -> WorkerResult:
        state.request_stop()  # the first molecule asks the job to stop
        return _fake()(task)

    summary = run_training(
        _settings(tmp_path, csv, workers=1),
        pool_factory=_threads,
        task_fn=stopping_task,
        install_signals=False,
    )
    assert summary.stopped is True and summary.attempted == 1
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["b", "status"] == "Pending"


def test_a_live_job_blocks_a_second_run(tmp_path: Path, make_csv: CsvFactory) -> None:
    import os

    csv = _csv(make_csv)
    JobState(tmp_path / "store" / "training.sqlite3").begin(
        pid=os.getpid(), workers=1, limit=None
    )
    with pytest.raises(JobAlreadyRunning):
        run_training(
            _settings(tmp_path, csv),
            pool_factory=_threads,
            task_fn=_fake(),
            install_signals=False,
        )


def test_unmigrated_csv_is_refused_before_any_work(
    tmp_path: Path, make_csv: CsvFactory
) -> None:
    csv = make_csv([row("a", "OK"), row("b", "Pending")])
    with pytest.raises(NotMigratedError):
        run_training(
            _settings(tmp_path, csv),
            pool_factory=_threads,
            task_fn=_fake(),
            install_signals=False,
        )
    assert JobState(tmp_path / "store" / "training.sqlite3").snapshot().state == "idle"


def test_empty_queue_finishes_cleanly(tmp_path: Path, make_csv: CsvFactory) -> None:
    csv = make_csv([row("old", "OK", protocol="1")])
    summary = run_training(
        _settings(tmp_path, csv), pool_factory=_threads, task_fn=_fake(), install_signals=False
    )
    assert summary.attempted == 0
```

- [ ] **Step 6: Rodar e ver falhar**

Run: `poetry run pytest tests/training/test_runner.py -v`
Expected: FAIL (`ModuleNotFoundError: semi_imperium.training.runner`).

- [ ] **Step 7: Implementar `runner.py`**

```python
"""Coordinator of the training routine: the single writer.

It reads the queue, dispatches one molecule per free worker, and records each
finished molecule: first the geometry (SQLite), then the CSV row, so a row
marked ``OK`` always has its geometry. Stopping is cooperative: no new
molecule is dispatched, in-flight ones finish, untouched rows stay ``Pending``.
"""

from __future__ import annotations

import multiprocessing
import os
import signal
import threading
import uuid
from collections import deque
from collections.abc import Callable
from concurrent.futures import (
    FIRST_COMPLETED,
    Executor,
    Future,
    ProcessPoolExecutor,
    wait,
)
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from semi_imperium.domain import RunRecord, RunState
from semi_imperium.persistence import SemiImperiumStore
from semi_imperium.training.db import database_path
from semi_imperium.training.geometry_store import GeometryRecord, GeometryStore
from semi_imperium.training.job_state import JobState
from semi_imperium.training.ledger import QueuedMolecule, TrainingLedger
from semi_imperium.training.protocol import (
    CANONICAL_PROTOCOL,
    HAMILTONIAN,
    canonical_configuration,
    canonical_settings,
)
from semi_imperium.training.worker import WorkerResult, WorkerTask, run_task
from semi_imperium.workflows.calculation import build_provenance

PoolFactory = Callable[[int], Executor]
TaskFn = Callable[[WorkerTask], WorkerResult]


@dataclass(frozen=True)
class RunSettings:
    """Inputs of one training run (execution knobs, not scientific choices)."""

    csv_path: Path
    store_root: Path
    work_dir: Path
    workers: int = 3
    limit: int | None = None
    crest_executable: str = "crest"
    xtb_executable: str = "xtb"
    mopac_executable: str = "mopac"
    crest_threads: int = 4
    crest_timeout_s: float = 3600.0
    mopac_timeout_s: float = 900.0
    max_reruns: int = 3
    keep_artifacts: bool = False

    def __post_init__(self) -> None:
        if self.workers < 1:
            raise ValueError(f"RunSettings.workers must be >= 1, got {self.workers}")
        if self.limit is not None and self.limit < 1:
            raise ValueError(f"RunSettings.limit must be >= 1, got {self.limit}")


@dataclass(frozen=True)
class RunSummary:
    """Outcome of one run."""

    attempted: int
    ok: int
    failed: int
    stopped: bool


def _spawn_pool(workers: int) -> Executor:
    return ProcessPoolExecutor(
        max_workers=workers, mp_context=multiprocessing.get_context("spawn")
    )


def _run_id() -> str:
    return f"run-{datetime.now().strftime('%Y%m%dT%H%M%S')}-{uuid.uuid4().hex[:6]}"


def run_training(
    settings: RunSettings,
    *,
    pool_factory: PoolFactory = _spawn_pool,
    task_fn: TaskFn = run_task,
    install_signals: bool = True,
) -> RunSummary:
    """Run the queue to exhaustion (or to ``limit``, or until stopped)."""
    ledger = TrainingLedger(settings.csv_path, max_reruns=settings.max_reruns)
    ledger.load()  # raises NotMigratedError before any state is touched
    queue = ledger.queue(settings.limit)

    db_path = database_path(settings.store_root)
    job = JobState(db_path)
    geometries = GeometryStore(db_path)
    job.begin(pid=os.getpid(), workers=settings.workers, limit=settings.limit)

    stop = threading.Event()
    if install_signals:
        signal.signal(signal.SIGTERM, lambda *_: stop.set())
        signal.signal(signal.SIGINT, lambda *_: stop.set())

    run_id = _run_id()
    batch_id = f"protocol{CANONICAL_PROTOCOL}-{run_id}"
    sem_settings = canonical_settings()
    store = SemiImperiumStore(settings.store_root)
    run = RunRecord(
        run_id=run_id,
        configuration=canonical_configuration(sem_settings),
        provenance=build_provenance(sem_settings),
        label="delta-entalpia-treinamento",
    )
    store.save_run(run)
    running = run.transition_to(RunState.RUNNING)
    store.save_run(running)

    ok = failed = 0
    attempted = 0
    stopped = False
    try:
        pending = deque(queue)
        in_flight: dict[Future[WorkerResult], QueuedMolecule] = {}
        with pool_factory(settings.workers) as pool:
            while pending or in_flight:
                halt = stop.is_set() or job.stop_requested()
                if halt and not stopped:
                    stopped = True
                    job.mark_stopping()
                while pending and len(in_flight) < settings.workers and not halt:
                    molecule = pending.popleft()
                    ledger.mark_running(molecule.mol_id)
                    in_flight[pool.submit(task_fn, _task(settings, run_id, molecule))] = molecule
                    attempted += 1
                if not in_flight:
                    break
                done, _ = wait(in_flight, timeout=30.0, return_when=FIRST_COMPLETED)
                for future in done:
                    molecule = in_flight.pop(future)
                    if _record(ledger, geometries, job, molecule, future, batch_id):
                        ok += 1
                    else:
                        failed += 1
                job.heartbeat(done=ok, failed=failed)
        run_state = (
            RunState.COMPLETED
            if failed == 0
            else RunState.FAILED
            if ok == 0 and attempted > 0
            else RunState.PARTIAL
        )
        store.save_run(
            running.transition_to(
                run_state,
                error_message="every molecule failed" if run_state is RunState.FAILED else None,
            )
        )
        job.finish("finished", "parado a pedido" if stopped else "concluído")
    except BaseException as exc:
        job.finish("failed", f"{type(exc).__name__}: {exc}")
        raise
    return RunSummary(attempted=attempted, ok=ok, failed=failed, stopped=stopped)


def _task(settings: RunSettings, run_id: str, molecule: QueuedMolecule) -> WorkerTask:
    return WorkerTask(
        mol_id=molecule.mol_id,
        smiles=molecule.smiles,
        charge=molecule.charge,
        multiplicity=molecule.multiplicity,
        run_id=run_id,
        store_root=str(settings.store_root),
        work_dir=str(settings.work_dir),
        crest_executable=settings.crest_executable,
        xtb_executable=settings.xtb_executable,
        mopac_executable=settings.mopac_executable,
        crest_threads=settings.crest_threads,
        crest_timeout_s=settings.crest_timeout_s,
        mopac_timeout_s=settings.mopac_timeout_s,
        keep_artifacts=settings.keep_artifacts,
    )


def _record(
    ledger: TrainingLedger,
    geometries: GeometryStore,
    job: JobState,
    molecule: QueuedMolecule,
    future: Future[WorkerResult],
    batch_id: str,
) -> bool:
    """Persist one finished molecule; return whether it was accepted."""
    try:
        result = future.result()
    except Exception as exc:  # worker process died or was killed
        reason = f"{type(exc).__name__}: {exc}"
        ledger.record_failure(molecule.mol_id, reason)
        job.record_failure(molecule.mol_id, reason)
        return False
    if not result.ok or result.hof_kcal_mol is None:
        ledger.record_failure(molecule.mol_id, result.error or result.state)
        job.record_failure(molecule.mol_id, result.error or result.state)
        return False
    geometries.save(
        GeometryRecord(
            mol_id=result.mol_id,
            smiles=result.smiles,
            hamiltonian=HAMILTONIAN,
            protocol=CANONICAL_PROTOCOL,
            signature=result.signature,
            calculation_id=result.calculation_id,
            conformer_index=result.conformer_index,
            elements=result.elements,
            coordinates=result.coordinates,
            hof_kcal_mol=result.hof_kcal_mol,
            lowest_frequency_cm1=result.lowest_frequency_cm1,
            minimum_confirmed=True,
        )
    )
    ledger.record_success(
        result.mol_id,
        hof_kcal_mol=result.hof_kcal_mol,
        conformers_evaluated=result.conformers_evaluated,
        elapsed_s=result.elapsed_s,
        batch_id=batch_id,
    )
    return True
```

- [ ] **Step 8: Rodar e ver passar**

Run: `poetry run pytest tests/training -v`
Expected: PASS (todos os testes de `tests/training`). Se `RunState` não aceitar `RUNNING → PARTIAL` ou `RUNNING → COMPLETED`, ler `ALLOWED_RUN_TRANSITIONS` em `semi_imperium/domain/records.py` e ajustar só o estado final escolhido.

- [ ] **Step 9: Portões de qualidade do pacote**

```bash
poetry run black src/semi_imperium/training tests/training
poetry run ruff check src/semi_imperium/training tests/training
poetry run python scripts/patch_rdkit_stubs.py
poetry run mypy src/semi_imperium/training --strict
```
Corrigir o que aparecer (a anotação `# type: ignore[arg-type]` dos testes só vale em `tests/`, que o mypy do projeto não checa).

- [ ] **Step 10: Commit**

```bash
git add src/semi_imperium/training tests/training
git commit -m "feat(training): worker por molécula e coordenador retomável com escritor único" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 9: CLI `python -m semi_imperium.training`

**Files:**
- Create: `src/semi_imperium/training/__main__.py`
- Test: `tests/training/test_cli.py`

**Interfaces:**
- Produces: `main(argv: Sequence[str] | None = None) -> int` com subcomandos:
  - `status` — imprime o progresso (contagens com denominador) e o estado do job; exit 0.
  - `migrate [--apply]` — ensaio por padrão; imprime o relatório; exit 0.
  - `run [--workers N] [--limit N] [--crest-threads N] [--crest-timeout S] [--mopac-timeout S] [--keep-artifacts]` — exit 0 ao concluir, 2 se o CSV não foi migrado, 3 se já há job vivo.
  - `stop` — pede parada; exit 0 se havia job ativo, 1 caso contrário.
  - Globais: `--csv PATH` (padrão `data/thermo_pm7.csv`), `--store PATH` (padrão `data/semi_imperium`), `--work-dir PATH` (padrão `<tmp>/semi_imperium_training`).

- [ ] **Step 1: Escrever os testes**

```python
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from semi_imperium.training.__main__ import main
from tests.training.conftest import CsvFactory, row


def _args(tmp_path: Path, csv: Path, *rest: str) -> list[str]:
    return ["--csv", str(csv), "--store", str(tmp_path / "store"), *rest]


def test_status_prints_progress(
    tmp_path: Path, make_csv: CsvFactory, capsys: pytest.CaptureFixture[str]
) -> None:
    csv = make_csv([row("a", "OK", protocol="1"), row("b", "Pending", protocol="")])
    assert main([*_args(tmp_path, csv), "status"]) == 0
    out = capsys.readouterr().out
    assert "protocolo 1" in out and "1 de 2" in out
    assert "parado" in out or "ocioso" in out


def test_migrate_is_a_dry_run_unless_applied(
    tmp_path: Path, make_csv: CsvFactory, capsys: pytest.CaptureFixture[str]
) -> None:
    csv = make_csv([row("a", "OK"), row("b", "Selected")])
    assert main([*_args(tmp_path, csv), "migrate"]) == 0
    assert "ensaio" in capsys.readouterr().out
    assert "protocol" not in pd.read_csv(csv, nrows=0).columns or (
        pd.read_csv(csv, dtype=str, keep_default_na=False)["protocol"].eq("").all()
    )
    assert main([*_args(tmp_path, csv), "migrate", "--apply"]) == 0
    saved = pd.read_csv(csv, dtype=str, keep_default_na=False).set_index("mol_id")
    assert saved.loc["a", "protocol"] == "1" and saved.loc["b", "status"] == "Pending"


def test_run_refuses_an_unmigrated_csv_with_exit_2(
    tmp_path: Path, make_csv: CsvFactory, capsys: pytest.CaptureFixture[str]
) -> None:
    csv = make_csv([row("a", "OK"), row("b", "Pending")])
    assert main([*_args(tmp_path, csv), "run", "--limit", "1"]) == 2
    assert "migrate --apply" in capsys.readouterr().err


def test_stop_without_a_job_exits_1(tmp_path: Path, make_csv: CsvFactory) -> None:
    csv = make_csv([row("a", "Pending", protocol="")])
    assert main([*_args(tmp_path, csv), "stop"]) == 1
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/training/test_cli.py -v`
Expected: FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Implementar**

```python
"""Command line of the training routine: ``python -m semi_imperium.training``."""

from __future__ import annotations

import argparse
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from semi_imperium.training.db import database_path
from semi_imperium.training.job_state import JobAlreadyRunning, JobState
from semi_imperium.training.ledger import NotMigratedError, read_progress
from semi_imperium.training.migration import migrate
from semi_imperium.training.runner import RunSettings, run_training

# src/semi_imperium/training/__main__.py -> repository root
_ROOT = Path(__file__).resolve().parents[3]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m semi_imperium.training",
        description="Delta Entalpia: rotina de treinamento (Semi-Imperium, protocolo 2).",
    )
    parser.add_argument("--csv", type=Path, default=_ROOT / "data" / "thermo_pm7.csv")
    parser.add_argument("--store", type=Path, default=_ROOT / "data" / "semi_imperium")
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "semi_imperium_training",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("status", help="mostra o progresso e o estado do job")
    mig = commands.add_parser("migrate", help="marca o protocolo 1 e devolve Selected à fila")
    mig.add_argument("--apply", action="store_true", help="grava (padrão: ensaio)")
    run = commands.add_parser("run", help="processa a fila com o protocolo 2")
    run.add_argument("--workers", type=int, default=3)
    run.add_argument("--limit", type=int, default=None)
    run.add_argument("--crest-threads", type=int, default=4)
    run.add_argument("--crest-timeout", type=float, default=3600.0)
    run.add_argument("--mopac-timeout", type=float, default=900.0)
    run.add_argument("--keep-artifacts", action="store_true")
    commands.add_parser("stop", help="pede ao job em execução que pare")
    return parser


def _status(args: argparse.Namespace) -> int:
    progress = read_progress(args.csv)
    total = progress.total

    def line(label: str, count: int) -> str:
        share = 0.0 if total == 0 else 100.0 * count / total
        return f"  {label}: {count} de {total} ({share:.1f}%)"

    print(f"Progresso de {args.csv}")
    print(line("protocolo 2 (canônico)", progress.protocol2))
    print(line("protocolo 1 (preliminar)", progress.protocol1))
    print(line("na fila (Pending)", progress.pending))
    print(line("selecionadas (Selected)", progress.selected))
    print(line("para repetir (Rerun)", progress.rerun))
    print(line("em execução (Running)", progress.running))
    print(line("descartadas (Skip)", progress.skipped))
    if not progress.migrated:
        print("  CSV ainda não migrado: rode 'migrate --apply'.")
    job = JobState(database_path(args.store)).snapshot()
    label = {"idle": "ocioso"}.get(job.state, job.state)
    print(f"Job: {label} (feitas {job.done}, falhas {job.failed}) {job.message}".rstrip())
    return 0


def _migrate(args: argparse.Namespace) -> int:
    report = migrate(args.csv, apply=args.apply)
    mode = "aplicado" if report.applied else "ensaio (nada foi gravado)"
    print(f"Migração {mode}: {report.rows} linhas")
    print(f"  marcadas como protocolo 1: {report.legacy_marked}")
    print(f"  Selected devolvidas a Pending: {report.selected_reset}")
    print(f"  já marcadas: {report.already_marked}")
    if report.backup_path is not None:
        print(f"  cópia de segurança: {report.backup_path}")
    return 0


def _run(args: argparse.Namespace) -> int:
    settings = RunSettings(
        csv_path=args.csv,
        store_root=args.store,
        work_dir=args.work_dir,
        workers=args.workers,
        limit=args.limit,
        crest_threads=args.crest_threads,
        crest_timeout_s=args.crest_timeout,
        mopac_timeout_s=args.mopac_timeout,
        keep_artifacts=args.keep_artifacts,
    )
    try:
        summary = run_training(settings)
    except NotMigratedError as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return 2
    except JobAlreadyRunning as exc:
        print(f"erro: {exc}", file=sys.stderr)
        return 3
    print(
        f"Concluído: {summary.attempted} tentadas, {summary.ok} ok, "
        f"{summary.failed} falhas{' (parado a pedido)' if summary.stopped else ''}"
    )
    return 0


def _stop(args: argparse.Namespace) -> int:
    if JobState(database_path(args.store)).request_stop():
        print("Parada solicitada: as moléculas em andamento terminam e o job encerra.")
        return 0
    print("Nenhum job em execução.", file=sys.stderr)
    return 1


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point; returns the process exit code."""
    args = _build_parser().parse_args(argv)
    handlers = {"status": _status, "migrate": _migrate, "run": _run, "stop": _stop}
    return handlers[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Rodar e ver passar**

Run: `poetry run pytest tests/training -v`
Expected: PASS. `poetry run python -m semi_imperium.training --help` lista os quatro subcomandos.

- [ ] **Step 5: Commit**

```bash
git add src/semi_imperium/training tests/training
git commit -m "feat(training): CLI status, migrate, run e stop" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 10: API do dashboard (status, start, stop) e remoção do catálogo

**Files:**
- Create: `src/grimperium/webapi/training.py`, `tests/webapi/test_training_api.py`
- Modify: `src/grimperium/webapi/app.py`, `models.py`, `config.py`
- Delete: `src/grimperium/webapi/catalog.py`; no `tests/webapi/test_app.py` apagar `test_routine_catalog_mirrors_semi_imperium_defaults` e atualizar `test_health_reports_read_only`

**Interfaces:**
- Produces (`models.py`): `Health(status: str, can_dispatch: bool)`; `TrainingProgress` (mesmos campos de `ledger.TrainingProgress`), `TrainingJob` (campos de `JobSnapshot` + `limit`), `TrainingFailure(mol_id, reason, at)`, `TrainingStatus(progress, job, geometries: int, failures: list[TrainingFailure])`, `TrainingStartRequest(workers: int = Field(3, ge=1, le=8), limit: int | None = Field(None, ge=1))`. Remover `RoutineCatalog` e correlatos.
- Produces (`training.py`): `TrainingLaunch(csv_path, store_root, workers, limit)`; `Launcher = Callable[[TrainingLaunch], None]`; `default_launcher(launch)`; `build_status(paths) -> TrainingStatus`; `start(paths, body, launcher) -> TrainingStatus`; `stop(paths) -> TrainingStatus`.
- Rotas: `GET /api/training/status`, `POST /api/training/start`, `POST /api/training/stop` (os POSTs exigem origem local). `create_app(paths=None, launcher=None)`. `DashboardPaths.default()` respeita `GRIMPERIUM_ROOT`.

- [ ] **Step 1: Escrever os testes**

`tests/webapi/test_training_api.py`:

```python
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from grimperium.webapi.app import create_app
from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.training import TrainingLaunch
from semi_imperium.training.db import database_path
from semi_imperium.training.job_state import JobState
from tests.training.conftest import row

SCHEMA_ROWS = [
    row("a", "OK", protocol="1"),
    row("b", "Pending", protocol=""),
    row("c", "Pending", protocol=""),
]


def _paths(tmp_path: Path, rows: list[dict[str, object]] | None = None) -> DashboardPaths:
    from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager

    schema = BatchCSVManager(csv_path=None).get_schema()
    pd.DataFrame(rows or SCHEMA_ROWS).reindex(columns=schema).to_csv(
        tmp_path / "pm7.csv", index=False
    )
    for name in ("exp.csv", "vs_cbs.csv"):
        (tmp_path / name).write_text("x\n")
    return DashboardPaths(
        pm7_csv=tmp_path / "pm7.csv",
        semi_imperium_store=tmp_path / "store",
        exp_csv=tmp_path / "exp.csv",
        exp_vs_cbs_csv=tmp_path / "vs_cbs.csv",
    )


class Recorder:
    def __init__(self) -> None:
        self.launches: list[TrainingLaunch] = []

    def __call__(self, launch: TrainingLaunch) -> None:
        self.launches.append(launch)


def _client(paths: DashboardPaths, recorder: Recorder) -> TestClient:
    return TestClient(create_app(paths, launcher=recorder))


def test_status_reports_progress_and_an_idle_job(tmp_path: Path) -> None:
    client = _client(_paths(tmp_path), Recorder())
    body = client.get("/api/training/status").json()
    assert body["progress"]["total"] == 3
    assert (body["progress"]["protocol1"], body["progress"]["pending"]) == (1, 2)
    assert body["progress"]["migrated"] is True
    assert body["job"]["state"] == "idle"
    assert body["geometries"] == 0 and body["failures"] == []
    assert not database_path(tmp_path / "store").exists()  # GET never creates it


def test_status_of_an_unmigrated_csv_is_not_an_error(tmp_path: Path) -> None:
    paths = _paths(tmp_path, [row("a", "OK"), row("b", "Pending")])
    body = _client(paths, Recorder()).get("/api/training/status").json()
    assert body["progress"]["migrated"] is False
    assert body["progress"]["ok_unmarked"] == 1


def test_status_with_a_missing_csv_is_503(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    paths.pm7_csv.unlink()
    response = _client(paths, Recorder()).get("/api/training/status")
    assert response.status_code == 503


def test_start_launches_with_the_requested_limits(tmp_path: Path) -> None:
    recorder = Recorder()
    response = _client(_paths(tmp_path), recorder).post(
        "/api/training/start", json={"workers": 2, "limit": 5}
    )
    assert response.status_code == 200
    (launch,) = recorder.launches
    assert (launch.workers, launch.limit) == (2, 5)
    assert launch.csv_path == tmp_path / "pm7.csv"


def test_start_is_refused_when_the_csv_is_not_migrated(tmp_path: Path) -> None:
    recorder = Recorder()
    paths = _paths(tmp_path, [row("a", "OK"), row("b", "Pending")])
    response = _client(paths, recorder).post("/api/training/start", json={})
    assert response.status_code == 409
    assert "migrate --apply" in response.json()["detail"]
    assert recorder.launches == []


def test_start_is_refused_while_a_job_is_alive(tmp_path: Path) -> None:
    recorder = Recorder()
    paths = _paths(tmp_path)
    JobState(database_path(paths.semi_imperium_store)).begin(
        pid=os.getpid(), workers=1, limit=None
    )
    response = _client(paths, recorder).post("/api/training/start", json={})
    assert response.status_code == 409
    assert recorder.launches == []


def test_a_stale_job_does_not_block_start(tmp_path: Path) -> None:
    recorder = Recorder()
    paths = _paths(tmp_path)
    JobState(database_path(paths.semi_imperium_store)).begin(
        pid=2_000_000_000, workers=1, limit=None
    )
    assert _client(paths, recorder).get("/api/training/status").json()["job"]["state"] == "stale"
    assert _client(paths, recorder).post("/api/training/start", json={}).status_code == 200
    assert len(recorder.launches) == 1


@pytest.mark.parametrize("body", [{"workers": 0}, {"workers": 9}, {"limit": 0}])
def test_start_validates_bounds(tmp_path: Path, body: dict[str, int]) -> None:
    recorder = Recorder()
    response = _client(_paths(tmp_path), recorder).post("/api/training/start", json=body)
    assert response.status_code == 422
    assert recorder.launches == []


def test_stop_flags_the_running_job(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    state = JobState(database_path(paths.semi_imperium_store))
    state.begin(pid=os.getpid(), workers=1, limit=None)
    response = _client(paths, Recorder()).post("/api/training/stop")
    assert response.status_code == 200
    assert state.stop_requested() is True


def test_stop_without_a_job_is_409(tmp_path: Path) -> None:
    response = _client(_paths(tmp_path), Recorder()).post("/api/training/stop")
    assert response.status_code == 409


def test_posts_from_foreign_origins_are_rejected(tmp_path: Path) -> None:
    recorder = Recorder()
    client = _client(_paths(tmp_path), recorder)
    foreign = client.post(
        "/api/training/start", json={}, headers={"Origin": "https://evil.example"}
    )
    assert foreign.status_code == 403
    local = client.post(
        "/api/training/start", json={}, headers={"Origin": "http://localhost:5173"}
    )
    assert local.status_code == 200
    assert len(recorder.launches) == 1


def test_health_says_the_api_can_dispatch(tmp_path: Path) -> None:
    body = _client(_paths(tmp_path), Recorder()).get("/api/health").json()
    assert body == {"status": "ok", "can_dispatch": True}
```

- [ ] **Step 2: Rodar e ver falhar**

Run: `poetry run pytest tests/webapi/test_training_api.py -v`
Expected: FAIL (`ImportError: cannot import name 'TrainingLaunch'`).

- [ ] **Step 3: Implementar modelos e configuração**

Em `models.py`: trocar `Health` por

```python
class Health(BaseModel):
    status: str
    can_dispatch: bool
```

acrescentar (importando `Field` de pydantic):

```python
class TrainingProgress(BaseModel):
    total: int
    protocol2: int
    protocol1: int
    ok_unmarked: int
    pending: int
    selected: int
    rerun: int
    running: int
    skipped: int
    migrated: bool


class TrainingJob(BaseModel):
    state: str
    pid: int | None
    started_at: str | None
    heartbeat_at: str | None
    finished_at: str | None
    workers: int
    limit: int | None
    done: int
    failed: int
    message: str
    stop_requested: bool


class TrainingFailure(BaseModel):
    mol_id: str
    reason: str
    at: str


class TrainingStatus(BaseModel):
    progress: TrainingProgress
    job: TrainingJob
    geometries: int
    failures: list[TrainingFailure]


class TrainingStartRequest(BaseModel):
    workers: int = Field(3, ge=1, le=8)
    limit: int | None = Field(None, ge=1)
```

e remover do arquivo `RoutineCatalog` e todas as classes que só ele usa (as que `catalog.py` importa).

Em `config.py`, `default()` passa a ler a raiz do ambiente:

```python
    @classmethod
    def default(cls, root: Path | None = None) -> DashboardPaths:
        """Paths of the active datasets under ``root/data``.

        ``GRIMPERIUM_ROOT`` points the API at another checkout-shaped folder
        (used to rehearse the routine on a copy of the data).
        """
        base = root or Path(os.environ.get("GRIMPERIUM_ROOT", _REPO_ROOT))
        data = base / "data"
        ...
```

(acrescentar `import os`; o resto do corpo continua igual, usando `data`.)

- [ ] **Step 4: Implementar `webapi/training.py`**

```python
"""Training endpoints: read status, launch and stop the routine subprocess.

The API never runs the science. ``start`` spawns the same command a person
would type (``python -m semi_imperium.training run``), detached, and returns;
the dashboard then polls ``status``.
"""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from fastapi import HTTPException

from grimperium.webapi import readers
from grimperium.webapi.config import DashboardPaths
from grimperium.webapi.models import (
    TrainingFailure,
    TrainingJob,
    TrainingProgress,
    TrainingStartRequest,
    TrainingStatus,
)

_MIGRATE_HINT = (
    "O CSV ainda não foi migrado para o protocolo 2. Rode no terminal: "
    "poetry run python -m semi_imperium.training migrate --apply"
)


@dataclass(frozen=True)
class TrainingLaunch:
    """What the launcher needs to start one run."""

    csv_path: Path
    store_root: Path
    workers: int
    limit: int | None


Launcher = Callable[[TrainingLaunch], None]


def default_launcher(launch: TrainingLaunch) -> None:
    """Start the routine detached; its output goes to ``training.log``."""
    command = [
        sys.executable,
        "-m",
        "semi_imperium.training",
        "--csv",
        str(launch.csv_path),
        "--store",
        str(launch.store_root),
        "run",
        "--workers",
        str(launch.workers),
    ]
    if launch.limit is not None:
        command += ["--limit", str(launch.limit)]
    launch.store_root.mkdir(parents=True, exist_ok=True)
    log = (launch.store_root / "training.log").open("ab")
    subprocess.Popen(  # noqa: S603 - fixed argv, no shell, local paths only
        command,
        stdin=subprocess.DEVNULL,
        stdout=log,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )


def build_status(paths: DashboardPaths) -> TrainingStatus:
    """Progress from the CSV, job and failures from SQLite."""
    from semi_imperium.training.db import database_path
    from semi_imperium.training.geometry_store import GeometryStore
    from semi_imperium.training.job_state import JobState
    from semi_imperium.training.ledger import read_progress

    try:
        progress = read_progress(paths.pm7_csv)
    except FileNotFoundError as exc:
        raise readers.DataUnavailable(paths.pm7_csv, "file not found") from exc
    db = database_path(paths.semi_imperium_store)
    state = JobState(db)
    snapshot = state.snapshot()
    return TrainingStatus(
        progress=TrainingProgress(**progress.__dict__),
        job=TrainingJob(
            state=snapshot.state,
            pid=snapshot.pid,
            started_at=snapshot.started_at,
            heartbeat_at=snapshot.heartbeat_at,
            finished_at=snapshot.finished_at,
            workers=snapshot.workers,
            limit=snapshot.limit,
            done=snapshot.done,
            failed=snapshot.failed,
            message=snapshot.message,
            stop_requested=snapshot.stop_requested,
        ),
        geometries=GeometryStore(db).count(),
        failures=[
            TrainingFailure(mol_id=mol, reason=reason, at=at)
            for mol, reason, at in state.recent_failures()
        ],
    )


def start(
    paths: DashboardPaths, body: TrainingStartRequest, launcher: Launcher
) -> TrainingStatus:
    """Launch the routine unless one is alive or the CSV is not migrated."""
    status = build_status(paths)
    if status.job.state in ("running", "stopping"):
        raise HTTPException(409, "Já existe um treinamento em execução.")
    if not status.progress.migrated:
        raise HTTPException(409, _MIGRATE_HINT)
    launcher(
        TrainingLaunch(
            csv_path=paths.pm7_csv,
            store_root=paths.semi_imperium_store,
            workers=body.workers,
            limit=body.limit,
        )
    )
    return build_status(paths)


def stop(paths: DashboardPaths) -> TrainingStatus:
    """Ask the running routine to stop dispatching new molecules."""
    from semi_imperium.training.db import database_path
    from semi_imperium.training.job_state import JobState

    if not JobState(database_path(paths.semi_imperium_store)).request_stop():
        raise HTTPException(409, "Nenhum treinamento em execução.")
    return build_status(paths)
```

- [ ] **Step 5: Ligar as rotas em `app.py` e remover o catálogo**

Em `app.py`: remover os imports de `routine_catalog` e `RoutineCatalog`, o endpoint `/api/routines/catalog` e a função `catalog`; atualizar a docstring do módulo para "Local dashboard API: reads the datasets and launches the training routine"; `Health(status="ok", can_dispatch=True)`; importar `Depends`, `HTTPException`, `training` e `TrainingStartRequest`, `TrainingStatus`; mudar a assinatura para `def create_app(paths: DashboardPaths | None = None, launcher: training.Launcher | None = None) -> FastAPI:` com `dispatch = launcher or training.default_launcher`; e acrescentar antes do `return app`:

```python
    allowed_origins = {
        f"http://{host}:{port}"
        for host in ("localhost", "127.0.0.1")
        for port in (5173, 4173, 8000)
    }

    def require_local_origin(request: Request) -> None:
        """Refuse browser POSTs from pages that are not the local dashboard."""
        origin = request.headers.get("origin")
        if origin is not None and origin not in allowed_origins:
            raise HTTPException(403, "Origem não permitida.")

    @app.get("/api/training/status")
    def training_status() -> TrainingStatus:
        return training.build_status(files)

    @app.post("/api/training/start", dependencies=[Depends(require_local_origin)])
    def training_start(body: TrainingStartRequest) -> TrainingStatus:
        return training.start(files, body, dispatch)

    @app.post("/api/training/stop", dependencies=[Depends(require_local_origin)])
    def training_stop() -> TrainingStatus:
        return training.stop(files)
```

Apagar `src/grimperium/webapi/catalog.py` (`git rm`). Em `tests/webapi/test_app.py`: apagar o teste do catálogo e trocar a asserção de saúde por `{"status": "ok", "can_dispatch": True}` (renomear o teste para `test_health_reports_dispatch_capability`).

- [ ] **Step 6: Rodar tudo da API**

Run: `timeout 300 poetry run pytest tests/webapi tests/training -v`
Expected: PASS. Conferir também `grep -rn "webapi.catalog\|routine_catalog\|RoutineCatalog" src tests` sem resultados.

- [ ] **Step 7: Portões e commit**

```bash
poetry run black src/grimperium/webapi src/semi_imperium/training tests/webapi tests/training
poetry run ruff check src/ tests/
poetry run mypy src/ --strict
git add -A src/grimperium/webapi tests/webapi
git commit -m "feat(webapi): rotas de treinamento (status, start, stop) e remoção do catálogo de rotinas" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 11: Dashboard — módulo Delta Entalpia e tela de Treinamento

**Files:**
- Modify: `apps/dashboard/src/{api.ts,useApi.ts,modules.ts,DeltaModule.tsx,ResultsView.tsx,SummaryView.tsx,App.tsx,App.test.tsx,style.css}`
- Rename: `GrimperiumModule.tsx` → `DeltaEnthalpyModule.tsx` (`git mv`)
- Create: `TrainingView.tsx`, `ComingSoonView.tsx`, `TrainingView.test.tsx`
- Delete: `RoutinesView.tsx`, `routines.ts`; `stats.test.ts` se testar só `routines.ts` (ver passo 1)

**Interfaces:**
- Consumes: rotas da Tarefa 10 (`GET /api/training/status`, `POST /api/training/start|stop`, `Health.can_dispatch`).
- Produces: módulo `delta-enthalpy` (`name: 'Delta Entalpia'`) com abas `summary`, `data` (`experimental`, `cbs`, `pm7`), `training`, `learning`, `calc`.

- [ ] **Step 1: Mapear o que o construtor de rotinas usa**

```bash
cd apps/dashboard
grep -rn "routines'\|from './routines'\|RoutinesView\|RoutineCatalog\|RoutineParameter\|RoutineProgram\|RoutineProperty" src | grep -v "^src/RoutinesView.tsx\|^src/routines.ts"
sed -n 1,12p src/stats.test.ts
```

Anotar cada ocorrência; todas serão removidas ou trocadas nos passos abaixo. Se `stats.test.ts` testar também `stats.ts` (não só `routines.ts`), remover apenas as importações e testes de `routines`.

- [ ] **Step 2: Escrever os testes (vitest) da nova tela**

`apps/dashboard/src/TrainingView.test.tsx`:

```tsx
import { cleanup, fireEvent, render, screen, waitFor } from '@testing-library/react'
import { afterEach, expect, test, vi } from 'vitest'
import { TrainingView } from './TrainingView'
import type { TrainingStatus } from './api'

const progress = {
  total: 29568, protocol2: 10, protocol1: 8232, ok_unmarked: 0, pending: 21000,
  selected: 0, rerun: 3, running: 1, skipped: 29, migrated: true,
}
const idle: TrainingStatus['job'] = { state: 'idle', pid: null, started_at: null, heartbeat_at: null, finished_at: null, workers: 0, limit: null, done: 0, failed: 0, message: '', stop_requested: false }

function status(overrides: Partial<TrainingStatus> = {}): TrainingStatus {
  return { progress, job: idle, geometries: 10, failures: [], ...overrides }
}

function json(data: unknown, code = 200): Response {
  return new Response(JSON.stringify(data), { status: code, headers: { 'Content-Type': 'application/json' } })
}

function mock(current: () => TrainingStatus) {
  const calls: { path: string; method: string; body: string | null }[] = []
  vi.stubGlobal('fetch', vi.fn((input: RequestInfo | URL, init?: RequestInit) => {
    const path = String(input)
    calls.push({ path, method: init?.method ?? 'GET', body: typeof init?.body === 'string' ? init.body : null })
    if (path === '/api/training/start' || path === '/api/training/stop') return Promise.resolve(json(current()))
    return Promise.resolve(json(current()))
  }))
  return calls
}

afterEach(() => { cleanup(); vi.unstubAllGlobals() })

test('mostra as contagens com denominador explícito', async () => {
  mock(() => status())
  render(<TrainingView />)
  const bar = await screen.findByRole('img', { name: /Andamento do treinamento/ })
  expect(bar).toHaveAccessibleName(/Protocolo 2: 10 de 29\.568/)
  expect(screen.getAllByText(/8\.232 preliminares/).length).toBeGreaterThan(0)
  expect(screen.getAllByText(/21\.003 na fila/).length).toBeGreaterThan(0)
})

test('inicia com workers e limite e mostra o estado do job', async () => {
  let current = status()
  const calls = mock(() => current)
  render(<TrainingView />)
  await screen.findByRole('img', { name: /Andamento do treinamento/ })
  fireEvent.change(screen.getByLabelText('Moléculas em paralelo'), { target: { value: '2' } })
  fireEvent.change(screen.getByLabelText('Limite de moléculas (opcional)'), { target: { value: '5' } })
  current = status({ job: { ...idle, state: 'running', workers: 2, limit: 5, pid: 123 } })
  fireEvent.click(screen.getByRole('button', { name: 'Iniciar treinamento' }))
  await waitFor(() => expect(calls.some((c) => c.path === '/api/training/start')).toBe(true))
  const start = calls.find((c) => c.path === '/api/training/start')!
  expect(start.method).toBe('POST')
  expect(JSON.parse(start.body ?? '{}')).toEqual({ workers: 2, limit: 5 })
  expect(await screen.findByText(/Em execução com 2 em paralelo/)).toBeInTheDocument()
  expect(screen.getByRole('button', { name: 'Parar' })).toBeEnabled()
  expect(screen.getByRole('button', { name: 'Iniciar treinamento' })).toBeDisabled()
})

test('bloqueia o início quando o CSV não foi migrado e mostra o comando', async () => {
  mock(() => status({ progress: { ...progress, migrated: false, ok_unmarked: 8232 } }))
  render(<TrainingView />)
  await screen.findByRole('img', { name: /Andamento do treinamento/ })
  expect(screen.getByRole('button', { name: 'Iniciar treinamento' })).toBeDisabled()
  expect(screen.getByText(/semi_imperium\.training migrate --apply/)).toBeInTheDocument()
})

test('lista as falhas recentes com o motivo', async () => {
  mock(() => status({ failures: [{ mol_id: 'mol_00042', reason: 'saddle: imaginary mode', at: '2026-10-12T10:00:00+00:00' }] }))
  render(<TrainingView />)
  expect(await screen.findByText('mol_00042')).toBeInTheDocument()
  expect(screen.getByText(/saddle: imaginary mode/)).toBeInTheDocument()
})

test('mostra o erro da API ao iniciar', async () => {
  vi.stubGlobal('fetch', vi.fn((input: RequestInfo | URL) =>
    Promise.resolve(String(input) === '/api/training/start' ? json({ detail: 'Já existe um treinamento em execução.' }, 409) : json(status()))))
  render(<TrainingView />)
  await screen.findByRole('img', { name: /Andamento do treinamento/ })
  fireEvent.click(screen.getByRole('button', { name: 'Iniciar treinamento' }))
  expect(await screen.findByRole('alert')).toHaveTextContent('Já existe um treinamento em execução.')
})
```

- [ ] **Step 3: Rodar e ver falhar**

Run: `cd apps/dashboard && npm test -- TrainingView`
Expected: FAIL (módulo `./TrainingView` inexistente).

- [ ] **Step 4: `api.ts` — tipos, `postJson` e `can_dispatch`**

Em `api.ts`: trocar `Health` por `{ status: string; can_dispatch: boolean }`; remover as interfaces `RoutineParameter`, `RoutineProgram`, `RoutineProperty`, `RoutineCatalog` (linhas ~154–190); acrescentar:

```ts
export interface TrainingProgress {
  total: number
  protocol2: number
  protocol1: number
  ok_unmarked: number
  pending: number
  selected: number
  rerun: number
  running: number
  skipped: number
  migrated: boolean
}

export interface TrainingJob {
  state: 'idle' | 'running' | 'stopping' | 'finished' | 'failed' | 'stale'
  pid: number | null
  started_at: string | null
  heartbeat_at: string | null
  finished_at: string | null
  workers: number
  limit: number | null
  done: number
  failed: number
  message: string
  stop_requested: boolean
}

export interface TrainingFailure { mol_id: string; reason: string; at: string }

export interface TrainingStatus {
  progress: TrainingProgress
  job: TrainingJob
  geometries: number
  failures: TrainingFailure[]
}
```

Refatorar `getJson` para compartilhar a leitura de erro (mantendo o texto das mensagens atuais) e criar `postJson`:

```ts
async function requestJson<T>(path: string, init: RequestInit): Promise<T> {
  let response: Response
  try {
    response = await fetch(path, init)
  } catch (error) {
    if (init.signal?.aborted) throw error
    throw new Error(`Não foi possível conectar à API (${path}). Verifique se o servidor está ativo.`, { cause: error })
  }
  if (!response.ok) {
    const body: unknown = await response.json().catch(() => null)
    const detail =
      typeof body === 'object' && body !== null && 'detail' in body
        ? String(body.detail)
        : response.statusText || 'erro desconhecido'
    const file =
      typeof body === 'object' && body !== null && 'path' in body ? ` (${String(body.path)})` : ''
    const hint = response.status === 500 || response.status === 502
      ? ' Verifique se a API está ativa na porta 8000.'
      : ''
    throw new Error(`${path}: ${response.status} — ${detail}${file}${hint}`)
  }
  return (await response.json()) as T
}

export function getJson<T>(path: string, signal: AbortSignal): Promise<T> {
  return requestJson<T>(path, { signal })
}

export function postJson<T>(path: string, body: unknown = {}): Promise<T> {
  return requestJson<T>(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  })
}
```

(Preservar exatamente o texto da `hint` e dos erros que `getJson` tem hoje; os testes existentes de erro dependem dele.)

- [ ] **Step 5: `useApi` com atualização periódica silenciosa**

Substituir `useApi.ts` por:

```ts
import { useCallback, useEffect, useState } from 'react'
import { getJson } from './api'

type ApiState<T> =
  | { status: 'loading' }
  | { status: 'error'; message: string }
  | { status: 'ready'; data: T }

interface Options {
  /** Atualiza a cada `pollMs` ms sem voltar ao estado "carregando". */
  pollMs?: number
}

export function useApi<T>(path: string, options: Options = {}): { state: ApiState<T>; reload: () => void } {
  const { pollMs } = options
  const [state, setState] = useState<ApiState<T>>({ status: 'loading' })
  const [revision, setRevision] = useState(0)

  useEffect(() => {
    const controller = new AbortController()
    const load = () => getJson<T>(path, controller.signal).then(
      (data) => setState({ status: 'ready', data }),
      (error: unknown) => {
        if (!controller.signal.aborted) {
          setState({ status: 'error', message: String(error instanceof Error ? error.message : error) })
        }
      },
    )
    void load()
    const timer = pollMs ? window.setInterval(() => { void load() }, pollMs) : undefined
    return () => {
      controller.abort()
      if (timer !== undefined) window.clearInterval(timer)
    }
  }, [path, revision, pollMs])

  const reload = useCallback(() => {
    setState({ status: 'loading' })
    setRevision((current) => current + 1)
  }, [])

  return { state, reload }
}
```

- [ ] **Step 6: Implementar `TrainingView.tsx`**

```tsx
import { useState } from 'react'
import { postJson } from './api'
import type { TrainingStatus } from './api'
import { PageHero, SectionHeading } from './StoryVisuals'
import { useApi } from './useApi'

const POLL_MS = 5000
const number = new Intl.NumberFormat('pt-BR')
const percent = new Intl.NumberFormat('pt-BR', { maximumFractionDigits: 1 })

const JOB_LABEL: Record<string, string> = {
  idle: 'Parado',
  running: 'Em execução',
  stopping: 'Parando (terminando as moléculas em andamento)',
  finished: 'Concluído',
  failed: 'Falhou',
  stale: 'Interrompido sem aviso',
}

export function TrainingView() {
  const { state, reload } = useApi<TrainingStatus>('/api/training/status', { pollMs: POLL_MS })
  const [workers, setWorkers] = useState('3')
  const [limit, setLimit] = useState('')
  const [busy, setBusy] = useState(false)
  const [actionError, setActionError] = useState<string | null>(null)

  const act = async (path: string, body?: unknown) => {
    setBusy(true)
    setActionError(null)
    try {
      await postJson<TrainingStatus>(path, body)
      reload()
    } catch (error) {
      setActionError(String(error instanceof Error ? error.message : error))
    } finally {
      setBusy(false)
    }
  }

  return (
    <div className="view">
      <PageHero
        variant="control" eyebrow="Treinamento · Semi-Imperium" title="Calcule as moléculas" accent="que ensinam o modelo."
        lead="Esta rotina leva cada molécula do SMILES até a entalpia de formação PM7 com o protocolo canônico e guarda a geometria otimizada. O painel só aciona e acompanha: o cálculo roda no back-end."
      />
      {state.status === 'loading' && <p role="status" className="loading-state">Carregando andamento…</p>}
      {state.status === 'error' && (
        <p role="alert" className="error-state">{state.message} <button type="button" onClick={reload}>Tentar de novo</button></p>
      )}
      {state.status === 'ready' && <Ready data={state.data} busy={busy} workers={workers} limit={limit}
        onWorkers={setWorkers} onLimit={setLimit} actionError={actionError}
        onStart={() => act('/api/training/start', { workers: Number(workers), limit: limit.trim() === '' ? null : Number(limit) })}
        onStop={() => act('/api/training/stop')} />}
    </div>
  )
}

interface ReadyProps {
  data: TrainingStatus
  busy: boolean
  workers: string
  limit: string
  actionError: string | null
  onWorkers: (value: string) => void
  onLimit: (value: string) => void
  onStart: () => void
  onStop: () => void
}

function Ready({ data, busy, workers, limit, actionError, onWorkers, onLimit, onStart, onStop }: ReadyProps) {
  const { progress, job } = data
  const queued = progress.pending + progress.selected + progress.rerun
  const segments = [
    { key: 'p2', label: 'Protocolo 2', count: progress.protocol2, className: 'training-bar__p2' },
    { key: 'p1', label: 'Preliminar (protocolo 1)', count: progress.protocol1, className: 'training-bar__p1' },
    { key: 'run', label: 'Em execução', count: progress.running, className: 'training-bar__run' },
    { key: 'queue', label: 'Na fila', count: queued, className: 'training-bar__queue' },
    { key: 'skip', label: 'Descartadas', count: progress.skipped, className: 'training-bar__skip' },
  ]
  const active = job.state === 'running' || job.state === 'stopping'
  const share = (count: number) => (progress.total === 0 ? 0 : (100 * count) / progress.total)
  const summary = `Andamento do treinamento: ${segments.map((s) => `${s.label}: ${number.format(s.count)} de ${number.format(progress.total)}`).join('; ')}`
  const canStart = !active && progress.migrated && !busy

  return (
    <>
      <section className="panel" aria-label="Andamento do treinamento">
        <SectionHeading eyebrow="Andamento" title={`${number.format(progress.protocol2)} de ${number.format(progress.total)} moléculas no protocolo 2`}
          badge={<span className="badge">{JOB_LABEL[job.state] ?? job.state}</span>}>
          {number.format(progress.protocol1)} preliminares do protocolo 1 ficam como verificação até serem recalculadas; {number.format(queued)} na fila.
        </SectionHeading>
        <div className="training-bar" role="img" aria-label={summary}>
          {segments.map((s) => s.count > 0 && <span key={s.key} className={s.className} style={{ width: `${share(s.count)}%` }} />)}
        </div>
        <ul className="training-legend">
          {segments.map((s) => (
            <li key={s.key}><i className={s.className} aria-hidden="true" />
              {s.key === 'p1' ? `${number.format(s.count)} preliminares` : s.key === 'queue' ? `${number.format(s.count)} na fila` : `${s.label}: ${number.format(s.count)}`}
              {' '}<small>({percent.format(share(s.count))}% de {number.format(progress.total)})</small>
            </li>
          ))}
        </ul>
      </section>

      <section className="panel" aria-label="Controle do treinamento">
        <SectionHeading eyebrow="Controle" title="Iniciar ou parar" />
        {!progress.migrated && (
          <p className="train-note">O CSV ainda não foi migrado para o protocolo 2. Rode no terminal: <code>poetry run python -m semi_imperium.training migrate --apply</code></p>
        )}
        <div className="training-controls">
          <label>Moléculas em paralelo
            <input type="number" min={1} max={8} value={workers} onChange={(e) => onWorkers(e.target.value)} disabled={active} />
          </label>
          <label>Limite de moléculas (opcional)
            <input type="number" min={1} value={limit} onChange={(e) => onLimit(e.target.value)} disabled={active} />
          </label>
          <button type="button" className="primary-button" disabled={!canStart} onClick={onStart}>Iniciar treinamento</button>
          <button type="button" className="ghost-button" disabled={!active || job.state === 'stopping' || busy} onClick={onStop}>Parar</button>
        </div>
        {active && <p role="status">Em execução com {job.workers} em paralelo: {number.format(job.done)} concluídas, {number.format(job.failed)} com falha nesta rodada.</p>}
        {job.state === 'stale' && <p role="status">O último treinamento foi interrompido sem aviso; as moléculas em andamento voltam para a fila no próximo início.</p>}
        {job.message && !active && <p role="status">{job.message}</p>}
        {actionError && <p role="alert" className="error-state">{actionError}</p>}
      </section>

      {data.failures.length > 0 && (
        <section className="panel" aria-label="Falhas recentes">
          <SectionHeading eyebrow="Falhas" title="Últimas moléculas que não passaram" />
          <ul className="training-failures">
            {data.failures.map((f) => (
              <li key={`${f.mol_id}-${f.at}`}><strong>{f.mol_id}</strong> <span>{f.reason}</span></li>
            ))}
          </ul>
        </section>
      )}
      <p className="train-note">{number.format(data.geometries)} geometrias finais guardadas no banco do treinamento.</p>
    </>
  )
}
```

Acrescentar ao `style.css` (usar as variáveis de cor que o arquivo já define; abrir o arquivo e reaproveitar os tokens de `.badge`/`.panel`; os nomes abaixo só dependem de `currentColor` e opacidade para não inventar tokens):

```css
.training-bar { display: flex; height: 14px; border-radius: 7px; overflow: hidden; background: color-mix(in srgb, currentColor 8%, transparent); margin: 12px 0; }
.training-bar span { display: block; height: 100%; }
.training-bar__p2, .training-legend .training-bar__p2 { background: #2f9e6b; }
.training-bar__p1, .training-legend .training-bar__p1 { background: #6b8fc9; }
.training-bar__run, .training-legend .training-bar__run { background: #d9a441; }
.training-bar__queue, .training-legend .training-bar__queue { background: #9aa3ad; }
.training-bar__skip, .training-legend .training-bar__skip { background: #c4574f; }
.training-legend { list-style: none; padding: 0; display: grid; gap: 4px; }
.training-legend i { display: inline-block; width: 10px; height: 10px; border-radius: 3px; margin-right: 8px; }
.training-controls { display: flex; flex-wrap: wrap; gap: 12px; align-items: end; }
.training-controls label { display: grid; gap: 4px; font-size: 0.9rem; }
.training-failures { list-style: none; padding: 0; display: grid; gap: 6px; }
.training-failures li { display: flex; gap: 10px; flex-wrap: wrap; }
```

(Conferir no navegador, passo 11, se as cores têm contraste suficiente nos dois temas; se o arquivo tiver tokens de cor próprios, trocar os hexadecimais por eles.)

- [ ] **Step 7: `ComingSoonView.tsx`**

```tsx
import { PageHero, SectionHeading } from './StoryVisuals'

/** Placeholder honesto de uma rotina cujo projeto ainda não foi construído. */
export function ComingSoonView({ eyebrow, title, accent, lead, depends }: {
  eyebrow: string; title: string; accent: string; lead: string; depends: string
}) {
  return (
    <div className="view">
      <PageHero variant="control" eyebrow={eyebrow} title={title} accent={accent} lead={lead} />
      <section className="panel" aria-label="Rotina em construção">
        <SectionHeading eyebrow="Status" title="Em breve" badge={<span className="badge">Em construção</span>}>{depends}</SectionHeading>
      </section>
    </div>
  )
}
```

- [ ] **Step 8: Reorganizar módulos e telas**

- `git mv apps/dashboard/src/GrimperiumModule.tsx apps/dashboard/src/DeltaEnthalpyModule.tsx` e reescrever:

```tsx
import { lazy } from 'react'
import { ComingSoonView } from './ComingSoonView'
import { SummaryView } from './SummaryView'
import type { ViewKey } from './SummaryView'
import type { ModuleViewProps } from './modules'

const ExperimentalView = lazy(() => import('./ExperimentalView').then((module) => ({ default: module.ExperimentalView })))
const CbsView = lazy(() => import('./CbsView').then((module) => ({ default: module.CbsView })))
const Pm7View = lazy(() => import('./Pm7View').then((module) => ({ default: module.Pm7View })))
const TrainView = lazy(() => import('./TrainView').then((module) => ({ default: module.TrainView })))
const TrainingView = lazy(() => import('./TrainingView').then((module) => ({ default: module.TrainingView })))

/** Telas do módulo Delta Entalpia: dados, treinamento, Delta Learning e cálculo. */
export function DeltaEnthalpyModule({ view, onNavigate }: ModuleViewProps) {
  switch (view as ViewKey) {
    case 'experimental': return <ExperimentalView />
    case 'cbs': return <CbsView />
    case 'pm7': return <Pm7View />
    case 'training': return <TrainingView />
    case 'learning': return <TrainView />
    case 'calc': return (
      <ComingSoonView eyebrow="Cálculo · Grimperium" title="Calcule uma molécula" accent="pelo nome."
        lead="Você dará o nome; a rotina reaproveita o valor já calculado ou roda o protocolo e aplica a correção do modelo treinado."
        depends="Depende do Delta Learning: o modelo corrigido ainda não existe." />
    )
    default: return <SummaryView onNavigate={onNavigate} />
  }
}
```

- `SummaryView.tsx`: `export type ViewKey = 'summary' | 'experimental' | 'cbs' | 'pm7' | 'training' | 'learning' | 'calc'`; o botão da linha ~170 passa a `onClick={() => onNavigate('training')}` com o texto `Acompanhar o treinamento →`.
- `TrainView.tsx`: trocar `eyebrow="Treinar · Modelo de correção Δ"` por `eyebrow="Delta Learning · Modelo de correção Δ"`; deixar o resto.
- `modules.ts`: trocar o import para `import { DeltaEnthalpyModule } from './DeltaEnthalpyModule'`; no `HUB`, `description` sem "rotinas de cálculo montadas aqui" (ex.: `'Hub de rotinas de cálculo de propriedades com métodos rápidos corrigidos por machine learning.'` pode ficar) e `tabs` só com `results` ('Resumos') e `about` ('Sobre'); substituir o módulo `grimperium` por:

```ts
  {
    key: 'delta-enthalpy',
    name: 'Delta Entalpia',
    property: 'Entalpia de formação (estado padrão)',
    description: 'PM7 corrigido por Δ-learning a partir de cálculos CBS: treinamento, Delta Learning e cálculo.',
    tabs: [
      { key: 'summary', label: 'Resumo', hint: 'Visão geral do módulo' },
      {
        key: 'data', label: 'Dados', hint: 'Tabelas usadas pelo módulo',
        children: [
          { key: 'experimental', label: 'Experimental', hint: 'Medidas de laboratório' },
          { key: 'cbs', label: 'CBS', hint: 'Referência de alta precisão' },
          { key: 'pm7', label: 'PM7', hint: 'Cálculo rápido (legado)' },
        ],
      },
      { key: 'training', label: 'Treinamento', hint: 'Semi-Imperium: calcula as moléculas PM7 do treino' },
      { key: 'learning', label: 'Delta Learning', hint: 'Treina o modelo de correção Δ' },
      { key: 'calc', label: 'Cálculo', hint: 'Grimperium: calcula uma molécula pelo nome' },
    ],
    Component: DeltaEnthalpyModule,
  },
```

- `DeltaModule.tsx`: remover o `lazy` de `RoutinesView`; `default: return <ResultsView onNavigate={onNavigate} />`.
- `ResultsView.tsx`: apagar o botão "Ir para Rotinas" e o `onNavigate` não usado (a prop sai da assinatura e de `DeltaModule`); trocar o parágrafo por: `O resumo reunirá, por propriedade, as moléculas calculadas. Acompanhe o treinamento em Delta Entalpia → Treinamento, na barra lateral.`
- `App.tsx`: o texto de saúde passa a `health.state.data.status === 'ok' && health.state.data.can_dispatch ? 'API ativa' : 'Resposta inesperada da API'`; o rodapé troca "o painel apenas dispara processos; os cálculos rodam no back-end" por "o painel aciona e acompanha as rotinas; os cálculos rodam no back-end".
- `git rm apps/dashboard/src/RoutinesView.tsx apps/dashboard/src/routines.ts` (e ajustar `stats.test.ts` conforme o passo 1).

- [ ] **Step 9: Atualizar `App.test.tsx`**

- No `mockApi`: remover `'/api/routines/catalog'` e a constante `catalog`; `'/api/health'` → `json({ status: 'ok', can_dispatch: true })` (e na linha ~139 também); acrescentar `'/api/training/status'` retornando um `TrainingStatus` como o do `TrainingView.test.tsx`.
- Trocar `renderGrimperium` por `renderDeltaEnthalpy` clicando em `{ name: /Delta Entalpia/ }`.
- Remover o teste `rotinas: monta, valida e salva uma rotina sem disparar`.
- Atualizar o texto `'API ativa · somente leitura'` para `'API ativa'`.
- Teste do hub: abas `['Resumos', 'Sobre']`, `Resumos` com `aria-current="page"`, sem o clique em "Ir para Rotinas" nem `Rotinas · Painel de controle`; terminar clicando em `Sobre`.
- Testes da barra lateral e do módulo: usar `Delta Entalpia` (breadcrumb `Delta Chemistry/Delta Entalpia`) e as abas `['Resumo', 'Dados', 'Treinamento', 'Delta Learning', 'Cálculo']`; no lugar de clicar em `Treinar`, clicar em `Delta Learning` e esperar `Cinco etapas, sem vazamento de dados`; acrescentar:

```tsx
  fireEvent.click(menu.getByRole('button', { name: 'Cálculo' }))
  expect(await screen.findByText('Em breve')).toBeInTheDocument()
  fireEvent.click(menu.getByRole('button', { name: 'Treinamento' }))
  expect(await screen.findByRole('img', { name: /Andamento do treinamento/ })).toBeInTheDocument()
```

- Em `AboutView.tsx`, se mencionar "Rotinas"/"Grimperium" como módulo (ver `grep -n "Rotina\|Grimperium" src/AboutView.tsx`), ajustar o texto para o novo nome mantendo o teste `Hoje só existe o módulo de entalpia de formação` verde (ou atualizar a frase nos dois lugares).

- [ ] **Step 10: Rodar os portões do front**

```bash
cd apps/dashboard
npm run lint && npm run typecheck && npm test && npm run build
```
Expected: tudo verde. Resolver avisos de variáveis não usadas deixadas pela remoção do construtor.

- [ ] **Step 11: Verificar no navegador (prévia)**

Subir as duas prévias (`grimperium-api`, `grimperium-dashboard` de `.claude/launch.json`) apontando a API para uma **cópia** dos dados: criar `/tmp/grimperium-rehearsal/data/` com cópia de `data/thermo_pm7.csv`, `thermo_exp_hf.csv`, `thermo_exp_hf_vs_cbs.csv` e `thermo_cbs_chon.csv`, e iniciar a API com `GRIMPERIUM_ROOT=/tmp/grimperium-rehearsal`. Abrir Delta Entalpia → Treinamento: conferir a barra, a legenda, o botão desabilitado com o comando de migração (CSV real ainda não migrado) e, no tema escuro e no claro, o contraste das cores. Capturar uma imagem.

- [ ] **Step 12: Commit**

```bash
git add -A apps/dashboard
git commit -m "feat(dashboard): módulo Delta Entalpia com tela de Treinamento; remove o construtor de rotinas" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 12: Documentação, portões completos e ensaio com ferramentas reais

**Files:**
- Modify: `CHANGELOG.md`, `apps/dashboard/README.md`, `docs/plans/2026-10-03-hub-restructure.md` (uma linha apontando para este plano)

- [ ] **Step 1: Documentar**

`CHANGELOG.md`, em `[Unreleased]` → `### Added` (criar a seção se não existir) e `### Changed`/`### Removed`, no formato de `.claude/CHANGELOG_TEMPLATE.md`:

- Rotina de treinamento Delta Entalpia (`semi_imperium.training`: protocolo 2 com CONFPASS + filtro, `migrate/run/status/stop`, geometrias em SQLite, coluna `protocol`, 62 colunas no CSV).
- API: `GET /api/training/status`, `POST /api/training/start|stop`; `Health.can_dispatch`; `GRIMPERIUM_ROOT`.
- Dashboard: módulo Delta Entalpia (Treinamento funcional; Delta Learning e Cálculo em breve); removido o construtor de rotinas e `/api/routines/catalog`.
- Registrar o resultado da linha de base da Tarefa 0 (o que passou, o que travou) e o que foi e não foi verificado.

`apps/dashboard/README.md`: atualizar a tabela de módulos (Delta Chemistry hub: Resumos e Sobre; Delta Entalpia: Resumo, Dados, Treinamento, Delta Learning, Cálculo), a frase "somente leitura" e acrescentar a seção:

```markdown
## Treinamento (Delta Entalpia)

O painel aciona e acompanha a rotina; o cálculo roda no back-end e também pelo terminal:

```bash
poetry run python -m semi_imperium.training status
poetry run python -m semi_imperium.training migrate            # ensaio
poetry run python -m semi_imperium.training migrate --apply    # grava (cria cópia .bak)
poetry run python -m semi_imperium.training run --workers 3 --limit 10
poetry run python -m semi_imperium.training stop
```

A migração é feita uma vez: marca as moléculas do protocolo 1 e devolve `Selected` à fila.
```

- [ ] **Step 2: Portões completos**

```bash
poetry run black src/ tests/
poetry run ruff check src/ tests/
poetry run python scripts/patch_rdkit_stubs.py
poetry run mypy src/ --strict
timeout 1500 poetry run pytest tests/ -v --cov=src/grimperium -x -q 2>&1 | tail -30
cd apps/dashboard && npm run lint && npm run typecheck && npm test && npm run build
```

Qualquer falha preexistente que não seja da branch: registrar com o teste e a mensagem no `CHANGELOG`/relatório final, sem "consertar" fora de escopo. Não declarar sucesso sem a saída de cada comando.

- [ ] **Step 3: Ensaio científico em cópia (sem tocar no `data/` real)**

```bash
REH=/tmp/grimperium-rehearsal
mkdir -p $REH/data && cp data/thermo_pm7.csv $REH/data/
poetry run python -m semi_imperium.training --csv $REH/data/thermo_pm7.csv --store $REH/data/semi_imperium migrate
poetry run python -m semi_imperium.training --csv $REH/data/thermo_pm7.csv --store $REH/data/semi_imperium migrate --apply
poetry run python -m semi_imperium.training --csv $REH/data/thermo_pm7.csv --store $REH/data/semi_imperium run --workers 1 --limit 1
poetry run python -m semi_imperium.training --csv $REH/data/thermo_pm7.csv --store $REH/data/semi_imperium status
```

Verificar e anotar:
1. A migração marcou 8.232 como protocolo 1 e devolveu 4.427 `Selected` a `Pending` (conferir com `status`).
2. A migração não alterou nenhum valor: comparar o CSV original com o migrado, coluna por coluna, ignorando só `status` e `protocol` e a molécula já processada no passo seguinte (rodar **antes** do `run`, ou excluir a linha processada):
```bash
poetry run python - <<'E'
import pandas as pd
a = pd.read_csv("data/thermo_pm7.csv", dtype=str, keep_default_na=False)
b = pd.read_csv("/tmp/grimperium-rehearsal/data/thermo_pm7.csv", dtype=str, keep_default_na=False)
assert list(a["mol_id"]) == list(b["mol_id"])
changed = [c for c in a.columns if c not in ("status",) and not a[c].equals(b[c])]
print("colunas alteradas além de status/protocol:", changed)
assert changed == [], changed
legacy = b["protocol"] == "1"
print("protocolo 1:", int(legacy.sum()), "| OK originais:", int((a["status"] == "OK").sum()))
E
```
Se `changed` listar colunas numéricas por reformatação de ponto flutuante (por exemplo `0.15234` virando `0.15234000000000001`), comparar numericamente (`pd.to_numeric(...).equals`) e registrar a reformatação no relatório; qualquer diferença de valor é um bloqueio.
3. A molécula processada: a linha ficou `OK` com `protocol = 2`, `H298_pm7` preenchido, `total_time`, `batch_id = protocol2-run-...`; o SQLite tem a geometria (`poetry run python -c "from semi_imperium.training.geometry_store import GeometryStore as G; g=G('/tmp/grimperium-rehearsal/data/semi_imperium/training.sqlite3'); print(g.count()); print(g.get(sorted(g.mol_ids())[0]).to_xyz())"`), a frequência mais baixa é positiva ou acima de −50 cm⁻¹ e `minimum_confirmed` é verdadeiro; o diretório `artifacts/<calculation_id>` foi apagado.
4. Se a molécula falhar (CREST, CONFPASS, MOPAC): ler o motivo em `status`/`training.log`, **não** mascarar; diagnosticar com `superpowers:systematic-debugging`. O porte do CONFPASS é experimental: se a falha vier dele, reportar ao usuário em vez de trocar a estratégia por conta própria.

- [ ] **Step 4: Ensaio pelo dashboard (ainda na cópia)**

Iniciar a API com `GRIMPERIUM_ROOT=/tmp/grimperium-rehearsal` e o dashboard; em Delta Entalpia → Treinamento: iniciar com limite 1 e 1 em paralelo, acompanhar o estado `Em execução` → `Concluído`, a barra e a contagem de geometrias mudarem, e testar `Parar` numa rodada com limite 3. Capturar imagens. Parar as prévias.

- [ ] **Step 5: Confirmação do usuário antes do `data/` real**

**Parar e perguntar ao usuário** (uma mensagem): "O ensaio na cópia passou/falhou assim: <resultados>. Posso rodar `migrate --apply` no `data/thermo_pm7.csv` real (cria cópia `.pre-protocol-*.bak`) e processar 1 molécula pelo dashboard?" Só prosseguir com a resposta afirmativa. Se confirmado: rodar `migrate` (ensaio), mostrar o relatório, depois `migrate --apply`, depois iniciar com limite 1 pelo dashboard e conferir os mesmos pontos do Passo 3 no arquivo real. Não rodar o lote completo.

- [ ] **Step 6: Commit final**

```bash
git add CHANGELOG.md apps/dashboard/README.md docs/plans/2026-10-03-hub-restructure.md
git commit -m "docs: changelog e README da rotina de treinamento Delta Entalpia" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
git status --short
git log --oneline main..HEAD
```

Entregar ao usuário: o resumo do que foi verificado e do que não foi, e a branch pronta para revisão (sem push).

---

## Acompanhamentos fora deste plano (registrar, não implementar)

- **Filtro do ML:** `ml/data_loader` treina com `status == 'OK' AND cbs_quality_flag == 'OK'`; com a coluna `protocol`, o projeto Delta Learning precisa decidir se as 8.232 do protocolo 1 entram (hoje entrariam).
- **Registros JSON do store:** cada molécula grava `data/semi_imperium/calculations/**.json`; só 2 estão versionados. Decidir se `data/semi_imperium/` entra no `.gitignore` (pendente da Etapa 0 do plano-mãe).
- **Recalcular as 8.232** com o protocolo 2 quando houver moléculas suficientes (decisão do usuário).
- **Envio único ao Google Planilhas** (projeto separado, pede confirmação antes de enviar qualquer arquivo).
- **Delta Learning** e **Grimperium (cálculo por nome)**.
- **Validar o porte do CONFPASS** contra a especificação antes de recomendar o lote completo.

## Self-Review (feito ao escrever)

- **Cobertura da spec:** protocolo 2 com CONFPASS + filtro (T1); 8.232 preservadas e marcadas, `Selected` → `Pending`, backup e ensaio (T5); só a geometria final em SQLite (T2, T3, T8, limpeza de artefatos em T8); CLI e dashboard acionam o mesmo comando (T9, T10, T11); progresso das 29 mil na aba (T11); renomeação para Delta Entalpia, abas reorganizadas, construtor removido (T11); retomada e `Running` órfão (T6, T7, T8); escritor único do CSV (T8); critério de sucesso (T12 passos 3–5); migração só com confirmação (T12 passo 5); lote completo fora de escopo.
- **Placeholders:** nenhum "TBD"; os pontos em que o código depende de uma verificação local (nome do arquivo de teste de `outcome_from_hamiltonian_result`, tokens de cor do CSS, texto do `AboutView`) trazem o comando exato para resolver.
- **Consistência de tipos:** `WorkerTask`/`WorkerResult` (T8) usados em `runner` e nos testes com os mesmos campos; `TrainingProgress` do ledger (T6) é espelhada por `TrainingProgress` do pydantic (T10) e do TypeScript (T11) com os mesmos 10 campos; `JobSnapshot` (T7) ↔ `TrainingJob` (T10/T11); `Health.can_dispatch` em T10 e T11.
