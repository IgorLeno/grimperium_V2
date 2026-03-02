---
name: grimperium-csv
description: Referência do schema CSV e regras de evolução para o Grimperium
user-invocable: false
---

# CSV Schema — Grimperium Reference

## Arquivo-chave

- `src/grimperium/crest_pm7/batch/csv_manager.py` — `BatchCSVManager`

## Schema Atual (Phase C — 54 colunas)

### Grupo 1: Identidade da Molécula
| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `name` | str | Identificador único da molécula |
| `smiles` | str | SMILES canônico |
| `inchi` | str | InChI string |

### Grupo 2: Status do Pipeline
| Coluna | Tipo | Valores |
|--------|------|---------|
| `status` | str | `PENDING`, `SELECTED`, `RUNNING`, `OK`, `RERUN`, `SKIP` |
| `crest_status` | str | `NOT_ATTEMPTED`, `PREOPT`, `SEARCH`, `OK`, `FAILED` |
| `mopac_status` | str | `NOT_ATTEMPTED`, `RUNNING`, `OK`, `FAILED` |

### Grupo 3: Resultados CREST
| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `crest_conformers_generated` | int | Total de confôrmeros gerados |
| `num_conformers_selected` | int | Confôrmeros selecionados para MOPAC |
| `crest_time` | float | Tempo CREST em segundos |
| `crest_error` | str | Mensagem de erro (se falhou) |

### Grupo 4: Resultados MOPAC
| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `H298_pm7` | float | Entalpia de formação mais estável (kcal/mol) |
| `abs_diff` | float | \|H298_pm7 - H298_ref\| |
| `rel_diff` | float | abs_diff / \|H298_ref\| |
| `H298_ref` | float | Referência experimental (kcal/mol) |
| `mopac_time` | float | Tempo total MOPAC agregado (s) |

### Grupo 5: Descritores PM7 (Phase C — 10 colunas)
| Coluna | Tipo | Unidade |
|--------|------|---------|
| `mopac_dipole_debye` | float | Debye |
| `mopac_ionization_potential_ev` | float | eV |
| `mopac_homo_ev` | float | eV |
| `mopac_lumo_ev` | float | eV |
| `mopac_gap_ev` | float | eV |
| `mopac_cosmo_area_a2` | float | Å² |
| `mopac_cosmo_volume_a3` | float | Å³ |
| `mopac_gradient_norm` | float | kcal/mol/Å |
| `mopac_num_scf_cycles` | int | — |
| `mopac_point_group` | str | e.g. "C1" |
| `mopac_time_s` | float | segundos |

## Operações Principais

### `write_row` — Nova molécula
```python
manager = BatchCSVManager(csv_path)
manager.write_row({
    "name": "aspirin",
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "status": "PENDING",
    # ... outros campos obrigatórios
})
```

### `update_row` — Atualizar após processamento
```python
# Atualiza apenas as colunas especificadas; mantém o resto
manager.update_row(
    molecule_name="aspirin",
    updates={
        "status": "OK",
        "H298_pm7": -145.3,
        "mopac_dipole_debye": 3.24,
        "mopac_homo_ev": -9.71,
        "mopac_lumo_ev": -0.83,
    }
)
```

## Regras de Evolução do Schema

### Adicionar Coluna Nova
1. Adicionar ao grupo correto em `RESULT_COLUMNS` (se for resultado) ou criar novo grupo
2. Adicionar valor padrão em `_empty_row()` (`np.nan` para float, `None` para str)
3. Atualizar lógica de `update_row` se necessário
4. **Não renomear colunas existentes** sem migration

### Remover Coluna
1. Remover de `RESULT_COLUMNS`
2. Adicionar ao histórico de breaking changes abaixo
3. Rodar migration script em dados existentes (se aplicável)

### Renomear Coluna
Tratar como: remover antiga + adicionar nova + migration de dados.

## Histórico de Breaking Changes

### Phase C (atual)
- **Removido:** `delta_1`, `delta_2`, `delta_3` (substituídas por `abs_diff`, `rel_diff`)
- **Renomeado:** `most_stable_hof` → `H298_pm7`
- **Adicionado:** 10 descritores MOPAC (`mopac_*`)
- **Adicionado:** `mopac_status`, `mopac_time`

### Phase B (histórico)
- Schema inicial com CREST + energia única

## Reset de Resultados

Quando um batch precisa ser reprocessado, `RESULT_COLUMNS` define quais colunas são zeradas:

```python
# Em BatchCSVManager — colunas zeradas no reset
RESULT_COLUMNS = [
    "crest_status", "crest_conformers_generated", "crest_time", "crest_error",
    "mopac_status", "num_conformers_selected", "mopac_time",
    "H298_pm7", "abs_diff", "rel_diff",
    "mopac_dipole_debye", "mopac_ionization_potential_ev",
    "mopac_homo_ev", "mopac_lumo_ev", "mopac_gap_ev",
    "mopac_cosmo_area_a2", "mopac_cosmo_volume_a3",
    "mopac_gradient_norm", "mopac_num_scf_cycles",
    "mopac_point_group", "mopac_time_s",
]
```
