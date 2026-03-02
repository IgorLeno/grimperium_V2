---
name: grimperium-mopac
description: Referência especializada para parsing de arquivos MOPAC no Grimperium
user-invocable: false
---

# MOPAC File Parsing — Grimperium Reference

## Arquivos-chave

- `src/grimperium/crest_pm7/mopac_descriptors.py` — Parser principal de `.aux`
- `src/grimperium/crest_pm7/mopac_optimizer.py` — Orquestrador de otimizações MOPAC

## Formatos de Arquivo MOPAC

| Extensão | Conteúdo |
|----------|----------|
| `.mop` | Input: geometria + keywords MOPAC (PM7, EPS=78.4, etc.) |
| `.out` | Output legível: energias, geometria final, erros |
| `.aux` | Output estruturado: todos os descritores, melhor para parsing |
| `.arc` | Archive: histórico de otimizações (geralmente ignorado) |

## Notação Fortran no `.aux`

O formato `.aux` usa notação Fortran `D` para expoentes — **não** notação científica padrão:

```
+0.577997D+02  →  57.7997   (D+02 = E+02 = × 10²)
-0.123456D-03  →  -0.000123456
```

**Conversão Python:**
```python
def parse_fortran_float(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))
```

## Descritores Disponíveis no `.aux`

| Descritor CSV | Fonte `.aux` | Unidade |
|---------------|--------------|---------|
| `mopac_dipole_debye` | `DIPOLE:DEBYE` | Debye |
| `mopac_ionization_potential_ev` | `IONIZATION_POTENTIAL` | eV |
| `mopac_homo_ev` | Calculado via `EIGENVALUES` + `NUM_ELECTRONS` | eV |
| `mopac_lumo_ev` | Calculado via `EIGENVALUES` + `NUM_ELECTRONS` | eV |
| `mopac_gap_ev` | `lumo - homo` | eV |
| `mopac_cosmo_area_a2` | `COSMO_AREA` | Å² |
| `mopac_cosmo_volume_a3` | `COSMO_VOLUME` | Å³ |
| `mopac_gradient_norm` | `GRADIENT_NORM` | kcal/mol/Å |
| `mopac_num_scf_cycles` | `NUM_SCF_CYCLES` | inteiro |
| `mopac_point_group` | `POINT_GROUP` | string (e.g., "C1") |
| `mopac_time_s` | `CPU_TIME` | segundos |

## Cálculo HOMO/LUMO

HOMO e LUMO não são extraídos diretamente — são calculados a partir de `EIGENVALUES` e `NUM_ELECTRONS`:

```python
# Padrão no .aux:
# EIGENVALUES[N]=  val1  val2  val3  ...  (um por orbital)
# NUM_ELECTRONS=   N_electrons

# HOMO = eigenvalue[N_electrons/2 - 1]  (0-indexed, alpha electrons)
# LUMO = eigenvalue[N_electrons/2]

def _calc_homo_lumo(eigenvalues: list[float], num_electrons: int) -> tuple[float, float]:
    homo_idx = num_electrons // 2 - 1
    lumo_idx = homo_idx + 1
    homo = eigenvalues[homo_idx]
    lumo = eigenvalues[lumo_idx] if lumo_idx < len(eigenvalues) else np.nan
    return homo, lumo
```

## Padrões de Regex para `.aux`

```python
import re

# Dipole (valor único em D notation)
dipole_match = re.search(r"DIPOLE:DEBYE\s*=\s*([+-]?\d+\.\d+[Dd][+-]\d+)", content)

# Eigenvalues (array multi-linha)
eigen_match = re.search(r"EIGENVALUES\[\d+\]\s*=\s*([\s\S]+?)(?=\n\w)", content)
# Split e converter cada valor com parse_fortran_float()

# Point group (string, não D notation)
pg_match = re.search(r"POINT_GROUP\s*=\s*(\S+)", content)
```

## Verificação de Convergência

Antes de extrair descritores, verifique se o cálculo convergiu:

```python
# No .out: linha "COMPUTATION SUCCESSFUL" indica sucesso
# No .aux: presença de HEAT_OF_FORMATION indica conclusão normal
# Falha: arquivo .out contém "FAILED TO CONVERGE" ou está vazio

def is_converged(aux_path: Path) -> bool:
    if not aux_path.exists() or aux_path.stat().st_size == 0:
        return False
    content = aux_path.read_text()
    return "HEAT_OF_FORMATION" in content
```

## Armadilhas Comuns

1. **D notation vs E notation** — nunca usar `float()` diretamente em valores `.aux`
2. **EIGENVALUES pode ser multi-linha** — usar regex com `[\s\S]+?` (não-guloso)
3. **NUM_ELECTRONS é o total** — dividir por 2 para obter índice de orbitais ocupados
4. **Arquivos `.arc` gerados automaticamente** — listados no `.claudeignore`, não processar
5. **COSMO_ valores** — só presentes se o cálculo usou solvente (keyword `EPS=`)
