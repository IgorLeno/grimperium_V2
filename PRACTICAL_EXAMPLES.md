# 💡 EXEMPLOS PRÁTICOS DE USO - PARA O CLAUDE CODE

## 📌 ANTES E DEPOIS (Snippets para copiar/colar)

---

## EXEMPLO 1: Configurar logging em um batch

### ❌ ANTES (sem logging)
```python
def execute_batch(self, batch_id: str, molecule_ids: List[str]):
    print(f"Starting batch {batch_id}")
    
    for mol_id in molecule_ids:
        print(f"Processing {mol_id}")
        # ... código aqui ...
        print("Done")
```

**Problema:** Sem rastreabilidade, sem timestamps, sem estrutura.

### ✅ DEPOIS (com FILE_2)
```python
from grimperium.crest_pm7.logging_enhancements import (
    setup_batch_logging,
    log_rdkit_start,
    log_rdkit_done,
)

def execute_batch(self, batch_id: str, molecule_ids: List[str]):
    # NOVO: Setup logging
    logger = setup_batch_logging(batch_id)
    logger.info(f"Starting batch: {batch_id}")
    
    for mol_id in molecule_ids:
        log_rdkit_start(logger, mol_id)
        # ... código RDKit ...
        log_rdkit_done(logger, mol_id, 
            nrotbonds=2, tpsa=45.5, aromatic_rings=1)
```

**Resultado:**
```
[13:04:26] [INFO] Starting batch: batch_0001
[13:04:26] [INFO] [mol_00001] 🧬 RDKit: Calculating descriptors...
[13:04:26] [INFO] [mol_00001]   ✓ nrotbonds=2.0, tpsa=45.5, aromatic_rings=1
```

---

## EXEMPLO 2: Usar paths centralizados

### ❌ ANTES (paths em /tmp)
```python
import tempfile
from pathlib import Path

def generate_conformers(mol_id: str):
    # Espalhado por vários lugares
    temp_dir = Path("/tmp/crest_pm7") / mol_id
    input_file = temp_dir / "input.xyz"
    output_file = temp_dir / "conformers.xyz"
    
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    # ... código aqui ...
    return str(output_file)
```

**Problemas:**
- Disperso em /tmp (não portável)
- Código duplicado em múltiplos arquivos
- Difícil limpar ou reorganizar
- Git fica poluído

### ✅ DEPOIS (centralizado com FILE_1)
```python
from grimperium.crest_pm7.paths import get_crest_temp_files

def generate_conformers(batch_id: str, mol_id: str):
    # NOVO: Centralizando paths
    crest_files = get_crest_temp_files(batch_id, mol_id)
    input_file = crest_files['input']
    output_file = crest_files['conformers']
    
    # Diretório já é criado automaticamente
    # ... código aqui ...
    return str(output_file)
```

**Benefícios:**
- Único lugar para definir estrutura
- Automaticamente cria diretórios
- Organizado: `./src/crest_pm7/tmp/batch_0001/mol_00001/`
- Fácil de limpar: `cleanup_batch("batch_0001")`

---

## EXEMPLO 3: Preencher CSV automaticamente

### ❌ ANTES (CSV vazio)
```python
def finish_molecule(mol_id: str, results: dict):
    # Cálculos feitos manualmente (ou não feitos!)
    
    # CSV fica assim:
    # mol_id, abs_diff, delta_1, delta_2, delta_3, ...
    # mol_001,        ,        ,        ,        , ...
    # ❌ Campos vazios!
    
    csv_df.loc[csv_df['mol_id'] == mol_id].update(results)
```

### ✅ DEPOIS (CSV completo com FILE_3)
```python
from grimperium.crest_pm7.csv_enhancements import (
    BatchSettingsCapture,
    CSVManagerExtensions,
)

def finish_molecule(batch_id: str, mol_id: str, 
                    h298_cbs: float, h298_pm7: float,
                    conformer_energies: list):
    
    # NOVO: Capturar settings uma vez
    batch_settings = BatchSettingsCapture.capture_batch_settings(self.pm7_config)
    
    # NOVO: Atualizar CSV automaticamente
    success = CSVManagerExtensions.update_molecule_with_mopac_results(
        csv_manager=self.csv_manager,
        mol_id=mol_id,
        h298_cbs=h298_cbs,
        h298_pm7=h298_pm7,
        mopac_hof_values=conformer_energies,
        batch_settings=batch_settings,
        batch_id=batch_id,
    )
    
    # CSV agora fica assim:
    # mol_id, abs_diff, abs_diff_%, delta_1, delta_2, delta_3, conformer_selected, v3, qm, ...
    # mol_001,     2.2,       12.57,   0.00,   0.45,   0.81,                0, T, F, ...
    # ✓ Completo!
```

---

## EXEMPLO 4: Suprimir warnings

### ❌ ANTES (com DtypeWarning)
```python
import pandas as pd

df = pd.read_csv("molecules.csv")

# Terminal cheio de warnings:
# DtypeWarning: Columns (abs_diff,delta_1,delta_2...) have mixed types.
# FutureWarning: Incompatible dtype for column...
```

### ✅ DEPOIS (limpo com FILE_2)
```python
import pandas as pd
from grimperium.crest_pm7.logging_enhancements import suppress_pandas_warnings

# NOVO: Uma linha no início da aplicação
suppress_pandas_warnings()

df = pd.read_csv("molecules.csv")

# Terminal limpo - sem warnings! ✓
```

---

## EXEMPLO 5: Calcular deltas (conceito importante!)

### 🔬 Entender o conceito

**Cenário:** Você tem uma molécula com 4 conformers otimizados pelo MOPAC.

```
Conformer 0: HOF = 0.42 kcal/mol  ← MELHOR (lowest energy)
Conformer 1: HOF = 0.87 kcal/mol
Conformer 2: HOF = 1.23 kcal/mol
Conformer 3: HOF = 1.89 kcal/mol
```

**Delta é a diferença em relação ao melhor:**

```
Δ1 = 0.42 - 0.42 = 0.00  ← O melhor sempre tem Δ = 0
Δ2 = 0.87 - 0.42 = 0.45  ← 2º melhor é 0.45 kcal/mol pior
Δ3 = 1.23 - 0.42 = 0.81  ← 3º melhor é 0.81 kcal/mol pior
```

### 📝 Usar em código

```python
from grimperium.crest_pm7.csv_enhancements import DeltaCalculations

# Seus dados de MOPAC
mopac_hof_values = [0.42, 0.87, 1.23, 1.89]

# Calcular deltas
delta_1, delta_2, delta_3, best_idx = \
    DeltaCalculations.calculate_deltas_and_select(mopac_hof_values)

print(f"Delta 1: {delta_1}")          # 0.0
print(f"Delta 2: {delta_2}")          # 0.45
print(f"Delta 3: {delta_3}")          # 0.81
print(f"Best conformer: {best_idx}")  # 0
```

---

## EXEMPLO 6: Pipeline completa (todas as 3 mudanças)

### 🚀 Caso de uso realista

```python
from grimperium.crest_pm7.paths import get_molecule_temp_dir, cleanup_batch
from grimperium.crest_pm7.logging_enhancements import (
    setup_batch_logging,
    log_rdkit_start, log_rdkit_done,
    log_crest_start, log_crest_done,
    log_mopac_start, log_mopac_done,
    suppress_pandas_warnings,
)
from grimperium.crest_pm7.csv_enhancements import (
    BatchSettingsCapture,
    CSVManagerExtensions,
)

# ==================== INÍCIO DO BATCH ====================
def run_batch_with_all_fixes(batch_id: str):
    
    # 1. Setup logging (FILE_2)
    logger = setup_batch_logging(batch_id)
    suppress_pandas_warnings()
    logger.info(f"🚀 Starting batch: {batch_id}")
    
    # 2. Capturar settings (FILE_3)
    batch_settings = BatchSettingsCapture.capture_batch_settings(self.pm7_config)
    
    for mol_id in ["mol_00001", "mol_00002", "mol_00003"]:
        
        # 3. Usar paths centralizados (FILE_1)
        mol_temp = get_molecule_temp_dir(batch_id, mol_id)
        logger.debug(f"Temp dir: {mol_temp}")
        
        # ==================== RDKit ====================
        log_rdkit_start(logger, mol_id)
        rdkit_result = self.run_rdkit(mol_id)
        log_rdkit_done(logger, mol_id, **rdkit_result)
        
        # ==================== CREST ====================
        log_crest_start(logger, mol_id)
        conformers = self.run_crest(batch_id, mol_id)
        log_crest_done(logger, mol_id, len(conformers), crest_time=4.2)
        
        # ==================== MOPAC ====================
        log_mopac_start(logger, mol_id, len(conformers))
        mopac_results = self.run_mopac(batch_id, mol_id, conformers)
        log_mopac_done(logger, mol_id, 
            best_conformer_idx=0,
            best_delta_energy=0.0,
            time_seconds=2.5
        )
        
        # ==================== UPDATE CSV ====================
        CSVManagerExtensions.update_molecule_with_mopac_results(
            csv_manager=self.csv_manager,
            mol_id=mol_id,
            h298_cbs=rdkit_result['h298_cbs'],
            h298_pm7=mopac_results['h298_pm7'],
            mopac_hof_values=mopac_results['conformer_energies'],
            batch_settings=batch_settings,
            batch_id=batch_id,
        )
    
    # ==================== LIMPEZA ====================
    cleanup_batch(batch_id)
    logger.info(f"✓ Batch completed successfully!")

# ==================== OUTPUT ====================
# [13:04:26] [INFO] 🚀 Starting batch: batch_0001
# [13:04:26] [DEBUG] Temp dir: ./src/crest_pm7/tmp/batch_0001/mol_00001
# [13:04:26] [INFO] [mol_00001] 🧬 RDKit: Calculating descriptors...
# [13:04:26] [INFO] [mol_00001]   ✓ nrotbonds=2.0, tpsa=45.5, aromatic_rings=1
# [13:04:26] [INFO] [mol_00001] 🔄 CREST: Starting conformer sampling...
# [13:04:30] [INFO] [mol_00001]   ✓ Generated 4 conformers in 4.2s
# [13:04:30] [INFO] [mol_00001] ⚛️  MOPAC: Optimizing 4 conformers...
# [13:04:32] [INFO] [mol_00001]   ✓ Selected conformer #0 with ΔE=0.0
# [13:04:32] [INFO] [mol_00001] Updated CSV with calculated deltas
# ... (repetir para mol_00002 e mol_00003)
# [13:04:45] [INFO] ✓ Batch completed successfully!
#
# CSV agora tem 11 campos preenchidos:
# abs_diff, abs_diff_%, delta_1, delta_2, delta_3,
# conformer_selected, v3, qm, nci, precise_scf, scf_threshold
```

---

## 🎯 PADRÕES PARA COPIAR/COLAR

### Padrão 1: Setup logging
```python
from grimperium.crest_pm7.logging_enhancements import setup_batch_logging

logger = setup_batch_logging(batch_id)
logger.info("Mensagem aqui")
logger.debug("Debug aqui")
logger.error("Erro aqui")
```

### Padrão 2: Registrar cada fase
```python
log_rdkit_start(logger, mol_id)
# ... código ...
log_rdkit_done(logger, mol_id, **props_dict)

log_crest_start(logger, mol_id)
# ... código ...
log_crest_done(logger, mol_id, num_conformers, time_s)

log_mopac_start(logger, mol_id, num_conformers)
# ... código ...
log_mopac_done(logger, mol_id, best_idx, delta, time_s)
```

### Padrão 3: Usar paths
```python
from grimperium.crest_pm7.paths import get_crest_temp_files, get_mopac_temp_files

files = get_crest_temp_files(batch_id, mol_id)
input_xyz = files['input']
output_xyz = files['conformers']

mopac_files = get_mopac_temp_files(batch_id, mol_id, conformer_idx)
mopac_input = mopac_files['input']
mopac_output = mopac_files['output']
```

### Padrão 4: Atualizar CSV
```python
from grimperium.crest_pm7.csv_enhancements import CSVManagerExtensions

success = CSVManagerExtensions.update_molecule_with_mopac_results(
    csv_manager=self.csv_manager,
    mol_id=mol_id,
    h298_cbs=-17.5,
    h298_pm7=-15.3,
    mopac_hof_values=[0.42, 0.87, 1.23],
    batch_settings=batch_settings,
    batch_id=batch_id,
)

if success:
    print("✓ CSV updated")
```

---

## ✅ CHECKLIST DURANTE IMPLEMENTAÇÃO

Ao integrar cada função, verifica se:

- [ ] Import funciona sem erros
- [ ] Função chamada com argumentos corretos
- [ ] Output aparece no log/console
- [ ] CSV é atualizado após chamada
- [ ] Sem warnings ou erros
- [ ] Teste com 1 molécula primeiro
- [ ] Depois teste com batch completo

---

## 📞 TROUBLESHOOTING RÁPIDO

| Problema | Solução |
|----------|---------|
| `ImportError: No module named 'grimperium.crest_pm7.paths'` | Verifique que `paths.py` foi copiado para `src/grimperium/crest_pm7/` |
| `DtypeWarning` ainda aparece | Chame `suppress_pandas_warnings()` no início da app |
| CSV não tem valores | Certifique que `update_molecule_with_mopac_results()` é chamado |
| Logs não aparecem | Chame `setup_batch_logging(batch_id)` antes de usar logger |
| `/tmp/crest_pm7` ainda é criado | Procure por `"/tmp"` no código e substitua por `get_*_temp_files()` |

---

## 🎉 PRONTO!

Você tem exemplos concretos e padrões prontos para copiar. Boa sorte na implementação!
