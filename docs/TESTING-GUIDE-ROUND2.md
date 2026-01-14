# 🧪 Testing Guide - Post Round 2 Fixes

**Data**: 2026-01-13  
**Status**: Ready to test  
**Objetivo**: Verificar que os 28 fixes (19+9) funcionam em produção

---

## 📋 PRÉ-REQUISITOS

✅ Round 1 (19 fixes) implementado  
✅ Round 2 (9 fixes) implementado  
✅ Linting passou (ruff)  
✅ Functional verification passou (9/9)  

---

## 🧪 FASE 1: SMOKE TEST (10 moléculas)

### Objetivo
Verificar que o sistema básico funciona com dataset mínimo.

### Comandos
```bash
# 1. Criar CSV de teste com 10 moléculas
python scripts/init_batch_csv.py \
  --input data/thermo_cbs_clean.csv \
  --output data/batch_test_10.csv \
  --limit 10

# 2. Verificar CSV criado

# Option 1: grep (reliable — counts non-empty lines)
grep -c . data/batch_test_10.csv
# Output: 11 (header + 10 molecules)

# Option 2: wc -l (counts newlines — may be off-by-1 without trailing newline)
wc -l data/batch_test_10.csv
# Output: 10 or 11 (depending on trailing newline in file)

# Option 3: Python (most reliable)
python -c "import pandas as pd; df = pd.read_csv('data/batch_test_10.csv'); print(f'Header + molecules: {len(df) + 1}')"
# Output: Header + molecules: 11

# Recommended: Use Option 1 (grep) or Option 3 (Python) for accuracy
head -n 2 data/batch_test_10.csv  # Ver header e primeira linha
```

### Verificações Esperadas
- ✅ CSV criado com 11 linhas (1 header + 10 moléculas)
- ✅ 33 colunas presentes
- ✅ status = "Pending" (todos)
- ✅ retry_count = 0 (todos)
- ✅ OK enum value = "OK" (uppercase) - FIX 2
- ✅ Nenhum erro de conversão - FIX 3
- ✅ Mensagens de erro determinísticas (se houver) - FIX 4

### Se Falhar
- Verificar se thermo_cbs_clean.csv existe
- Verificar se colunas "smiles" e "nheavy" existem
- Verificar logs de erro (sorted columns)

---

## 🧪 FASE 2: SCALE TEST (40 moléculas)

### Objetivo
Testar state transitions e JSON detail files.

### Comandos
```bash
# 1. Criar CSV de teste com 40 moléculas
python scripts/init_batch_csv.py \
  --input data/thermo_cbs_clean.csv \
  --output data/batch_test_40.csv \
  --limit 40

# 2. [Aqui você rodaria o batch processing]
# Exemplo (quando tiver CLI):
# python -m grimperium.cli batch-run \
#   --csv data/batch_test_40.csv \
#   --batch-size 10 \
#   --crest-timeout 30 \
#   --mopac-timeout 60
```

### Verificações Esperadas
- ✅ State transitions corretas:
  - PENDING → SELECTED
  - SELECTED → RUNNING
  - RUNNING → OK/RERUN/SKIP
- ✅ JSON detail files criados em data/conformer_details/
- ✅ File descriptors não vazam - FIX 1
- ✅ _safe_int funciona com NaN - FIX 3, 5
- ✅ CSV atualizado corretamente
- ✅ retry_count incrementado corretamente - FIX 5
- ✅ DEFAULT_MAX_OBSERVATIONS respeitado - FIX 6

### Monitoramento
```bash
# Durante execução, monitorar file descriptors (process-specific)

# Option 1: Monitor specific grimperium process
watch -n 1 'lsof -p $(pgrep -f "grimperium.cli" | head -1) 2>/dev/null | wc -l'

# Option 2: More specific pattern
watch -n 1 'lsof -p $(pgrep -f "batch-run" | head -1) 2>/dev/null | wc -l'

# Option 3: If multiple processes, sum all
watch -n 1 'for pid in $(pgrep -f "grimperium.cli"); do lsof -p $pid 2>/dev/null; done | wc -l'

# Option 4: Get PID first, then monitor
# In one terminal:
ps aux | grep grimperium.cli | grep -v grep  # Get PID
# In another terminal:
watch -n 1 'lsof -p <PID> 2>/dev/null | wc -l'  # Replace <PID>

# Verificar se JSON detail files existem
ls -l data/conformer_details/*.json | wc -l

# Verificar status counts
cut -d',' -f10 data/batch_test_40.csv | sort | uniq -c
```

### Se Falhar
- Verificar logs de erro
- Verificar se file descriptors estão vazando (lsof)
- Verificar se JSON files estão corrompidos
- Verificar se retry_count está sendo incrementado

---

## 🧪 FASE 3: PERFORMANCE ANALYSIS (100→1k→5k)

### Objetivo
Medir tempo por molécula, memória, file descriptor count.

### Comandos
```bash
# 1. Teste com 100 moléculas
python scripts/init_batch_csv.py \
  --input data/thermo_cbs_clean.csv \
  --output data/batch_test_100.csv \
  --limit 100

# 2. Monitorar recursos (platform-aware)

# Linux (GNU time with verbose output)
/usr/bin/time -v python -m grimperium.cli batch-run \
  --csv data/batch_test_100.csv \
  --batch-size 20 \
  > performance_100.log 2>&1

# macOS/BSD (Option 1: Use built-in time, less verbose)
time python -m grimperium.cli batch-run \
  --csv data/batch_test_100.csv \
  --batch-size 20 \
  > performance_100.log 2>&1

# macOS/BSD (Option 2: Install GNU time first)
# brew install gnu-time
gtime -v python -m grimperium.cli batch-run \
  --csv data/batch_test_100.csv \
  --batch-size 20 \
  > performance_100.log 2>&1

# Universal Python approach (works everywhere)
python << 'EOF'
import time, subprocess, sys

start = time.time()
proc = subprocess.run([
    sys.executable, '-m', 'grimperium.cli', 'batch-run',
    '--csv', 'data/batch_test_100.csv',
    '--batch-size', '20',
], capture_output=True, text=True)
elapsed = time.time() - start

# Save output
with open('performance_100.log', 'w') as f:
    f.write(proc.stdout)
    f.write(proc.stderr)

print(f'\n{"="*50}')
print(f'Execution time: {elapsed:.2f}s')
print(f'Exit code: {proc.returncode}')
print(f'{"="*50}')

sys.exit(proc.returncode)
EOF

# 3. Extrair métricas
grep "Maximum resident set size" performance_100.log
grep "File descriptors" performance_100.log
```

### Verificações Esperadas
- ✅ Tempo por molécula: < 5 minutos (média)
- ✅ Memória: < 2GB RSS
- ✅ File descriptors: < 100 (não vaza) - FIX 1
- ✅ Observations bounded: max 100 - FIX 6
- ✅ Nenhum crash em CSV corrompido - FIX 3

### Métricas para Coletar (com validação segura)

```python
import pandas as pd
import logging

LOG = logging.getLogger(__name__)

def load_and_analyze_batch(csv_path: str) -> dict:
    """Load batch CSV and compute statistics safely.
    
    Handles missing columns, invalid types, empty files gracefully.
    
    Args:
        csv_path: Path to batch CSV
        
    Returns:
        dict with keys: total, ok, rerun, skip, avg_time, success_rate, etc.
        Returns empty dict on error.
    """
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        LOG.error(f"Failed to load CSV {csv_path}: {e}")
        return {}
    
    # Validate not empty
    if df.empty:
        LOG.warning(f"CSV {csv_path} is empty")
        return {}
    
    # Validate required columns
    required_cols = ['status', 'total_execution_time']
    missing_cols = set(required_cols) - set(df.columns)
    if missing_cols:
        LOG.error(f"Missing required columns: {missing_cols}")
        return {}
    
    try:
        # Coerce total_execution_time to numeric (handles NaN, invalid)
        df['total_execution_time'] = pd.to_numeric(
            df['total_execution_time'],
            errors='coerce'
        )
        
        # Count statuses
        status_counts = df['status'].value_counts()
        
        stats = {
            'total': len(df),
            'ok': status_counts.get('OK', 0),
            'rerun': status_counts.get('Rerun', 0),
            'skip': status_counts.get('Skip', 0),
            'pending': status_counts.get('Pending', 0),
        }
        
        # Calculate timing statistics (safe for NaN)
        valid_times = df['total_execution_time'].dropna()
        if len(valid_times) > 0:
            stats['avg_time'] = float(valid_times.mean())
            stats['total_time'] = float(valid_times.sum())
            stats['min_time'] = float(valid_times.min())
            stats['max_time'] = float(valid_times.max())
        else:
            stats['avg_time'] = 0.0
            stats['total_time'] = 0.0
            stats['min_time'] = None
            stats['max_time'] = None
        
        # Calculate success rate (safe for zero total)
        if stats['total'] > 0:
            stats['success_rate'] = 100.0 * stats['ok'] / stats['total']
        else:
            stats['success_rate'] = 0.0
        
        return stats
        
    except Exception as e:
        LOG.error(f"Failed to analyze batch: {e}")
        return {}

# Usage
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    
    stats = load_and_analyze_batch('data/batch_test_100.csv')
    
    if stats:
        print(f"Total molecules: {stats['total']}")
        print(f"OK: {stats['ok']} ({stats['success_rate']:.1f}%)")
        print(f"Rerun: {stats['rerun']}")
        print(f"Skip: {stats['skip']}")
        print(f"Avg time: {stats['avg_time']:.2f}s")
        print(f"Total time: {stats['total_time']:.1f}s")
    else:
        print("Failed to load/analyze batch")
```

---

## 🧪 FASE 4: FULL PRODUCTION (30k moléculas)

### Objetivo
Rodar produção completa com monitoramento.

### Comandos
```bash
# 1. Criar CSV completo
python scripts/init_batch_csv.py \
  --input data/thermo_cbs_clean.csv \
  --output data/batch_30k.csv

#!/bin/bash
# Production batch orchestration for 30k molecules
# Dynamically calculates number of batches needed based on CSV size and batch_size

set -e  # Exit on error

CSV_PATH="data/batch_30k.csv"
BATCH_SIZE=500
CREST_TIMEOUT=30
MOPAC_TIMEOUT=60
LOG_FILE="production_30k.log"

# Validate CSV exists
if [ ! -f "$CSV_PATH" ]; then
    echo "ERROR: CSV file not found: $CSV_PATH"
    exit 1
fi

# Calculate total molecules
TOTAL_ROWS=$(wc -l < "$CSV_PATH")
TOTAL_MOLECULES=$((TOTAL_ROWS - 1))  # Subtract header row

# Calculate number of batches needed (ceiling division)
NUM_BATCHES=$(( (TOTAL_MOLECULES + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting batch orchestration"
echo "  CSV: $CSV_PATH"
echo "  Total molecules: $TOTAL_MOLECULES"
echo "  Batch size: $BATCH_SIZE"
echo "  Expected batches: $NUM_BATCHES"
echo "  Log file: $LOG_FILE"
echo ""

# Initialize log
cat > "$LOG_FILE" << EOF
=== Batch Orchestration Start ===
Date: $(date)
CSV: $CSV_PATH
Total molecules: $TOTAL_MOLECULES
Batch size: $BATCH_SIZE
Expected batches: $NUM_BATCHES

EOF

# Track statistics
TOTAL_PROCESSED=0
FAILED_BATCHES=0

# Loop through batches
for batch_num in $(seq 1 $NUM_BATCHES); do
    batch_id="batch_30k_$(printf '%03d' $batch_num)"
    
    # Display progress
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running batch $batch_num/$NUM_BATCHES (ID: $batch_id)"
    echo "  Expected molecules in this batch: $BATCH_SIZE (actual may be less for last batch)"
    
    # Log batch start
    echo "[Batch $batch_num/$NUM_BATCHES] $batch_id started at $(date)" >> "$LOG_FILE"
    
    # Run batch execution
    if python -m grimperium.cli batch-run \
        --csv "$CSV_PATH" \
        --batch-id "$batch_id" \
        --batch-size $BATCH_SIZE \
        --crest-timeout $CREST_TIMEOUT \
        --mopac-timeout $MOPAC_TIMEOUT \
        >> "$LOG_FILE" 2>&1; then
        
        echo "  ✅ Batch completed successfully"
        echo "[Batch $batch_num/$NUM_BATCHES] $batch_id completed successfully at $(date)" >> "$LOG_FILE"
        
        TOTAL_PROCESSED=$((TOTAL_PROCESSED + BATCH_SIZE))
    else
        echo "  ⚠️  Batch failed or encountered error"
        echo "[Batch $batch_num/$NUM_BATCHES] $batch_id FAILED at $(date)" >> "$LOG_FILE"
        
        FAILED_BATCHES=$((FAILED_BATCHES + 1))
    fi
    
    # Checkpoint: wait 5 seconds between batches
    if [ $batch_num -lt $NUM_BATCHES ]; then
        echo "  Waiting 5 seconds before next batch..."
        sleep 5
    fi
done

# Summary
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Batch orchestration completed"
echo "  Total batches: $NUM_BATCHES"
echo "  Failed batches: $FAILED_BATCHES"
echo "  Processed molecules: ~$TOTAL_PROCESSED"
echo "  Log: $LOG_FILE"

# Log summary
cat >> "$LOG_FILE" << EOF

=== Batch Orchestration Summary ===
Total batches: $NUM_BATCHES
Failed batches: $FAILED_BATCHES
Processed molecules: ~$TOTAL_PROCESSED
Completed at: $(date)

EOF

# Exit with error if any batch failed
if [ $FAILED_BATCHES -gt 0 ]; then
    echo "⚠️  WARNING: $FAILED_BATCHES batches failed"
    exit 1
fi

echo "✅ All batches completed successfully!"
exit 0
```

### Monitoramento Contínuo
```bash
# Terminal 1: Tail logs
tail -f production_30k.log

# Terminal 2: Monitorar recursos
watch -n 10 'ps aux | grep python | grep -v grep'

# Terminal 3: Status dashboard (no nested quotes)

# Option 1: Simple (no nested quotes)
watch -n 30 "cut -d, -f10 data/batch_30k.csv | sort | uniq -c"

# Option 2: Using awk (more robust)
watch -n 30 "awk -F, '{print \$10}' data/batch_30k.csv | sort | uniq -c"

# Option 3: Full path with explicit column name
watch -n 30 "cat data/batch_30k.csv | cut -d, -f10 | sort | uniq -c"
```

### Verificações Esperadas
- ✅ Nenhum crash por file descriptor leak - FIX 1
- ✅ CSVs antigos carregam corretamente - FIX 2
- ✅ CSV corrompido não causa crash - FIX 3
- ✅ Mensagens de erro determinísticas - FIX 4
- ✅ NaN handling correto - FIX 3, 5
- ✅ Memory bounded - FIX 6
- ✅ Type safety (sem exceções de tipo) - FIX 7, 8
- ✅ Documentação completa - FIX 9

### Análise Final
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/batch_30k.csv')

# Success rate
success_rate = (df['status'] == 'OK').sum() / len(df) * 100
print(f"Success rate: {success_rate:.1f}%")

# Status distribution
print("\nStatus distribution:")
print(df['status'].value_counts())

# Timing statistics
print("\nTiming statistics:")
print(df['total_execution_time'].describe())

# Quality grades
print("\nQuality grades:")
print(df['quality_grade'].value_counts())

# Plot: nheavy vs execution time
plt.scatter(df['nheavy'], df['total_execution_time'], alpha=0.3)
plt.xlabel('Number of heavy atoms')
plt.ylabel('Total execution time (s)')
plt.title('Execution time vs molecule size')
plt.savefig('execution_time_vs_nheavy.png')
```

---

## 🐛 DEBUGGING

### Se File Descriptor Leak (FIX 1)
```bash
# Monitorar FDs em tempo real
watch -n 1 'lsof -p $(pgrep -f python) | grep -E "json|tmp" | wc -l'

# Verificar se temp files estão sendo limpos
find data/conformer_details/ -name '.tmp_*' -ls
```

### Se CSV não Carrega (FIX 2)
```python
# Verificar enum values
from grimperium.crest_pm7.batch.enums import MoleculeStatus
print(f"OK value: '{MoleculeStatus.OK.value}'")
assert MoleculeStatus.OK.value == "OK", "Must be uppercase!"

# Verificar CSV
import pandas as pd
df = pd.read_csv('data/batch_test_10.csv')
print(df['status'].unique())
```

### Se Crash em CSV Corrompido (FIX 3)
```python
# Testar _safe_int com valores problemáticos
from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager
from pathlib import Path
import pandas as pd

mgr = BatchCSVManager(Path('/tmp/test.csv'))

# Test cases
test_values = [
    (pd.NA, "NaN"),
    (None, "None"),
    (float('nan'), "float NaN"),
    ("10", "string '10'"),
    ("3.5", "string '3.5'"),
    ("abc", "invalid string"),
    ("", "empty string"),
    ("  42  ", "whitespace string"),
]

for val, desc in test_values:
    result = mgr._safe_int(val, default=0)
    print(f"{desc:20s} → {result}")
```

---

## ✅ SUCCESS CRITERIA

### Smoke Test (Fase 1)
- [x] CSV criado com estrutura correta
- [x] 10 moléculas carregadas
- [x] Nenhum erro de parsing

### Scale Test (Fase 2)
- [ ] State transitions corretas
- [ ] JSON detail files criados
- [ ] File descriptors não vazam
- [ ] CSV atualizado corretamente

### Performance (Fase 3)
- [ ] Tempo por molécula < 5 min
- [ ] Memória < 2GB
- [ ] File descriptors < 100

### Production (Fase 4)
- [ ] 30k moléculas processadas
- [ ] Success rate > 80%
- [ ] Nenhum crash
- [ ] Documentação completa

---

## 📞 SUPORTE

Se encontrar problemas:

1. **Verificar logs**: `tail -f production_30k.log`
2. **Verificar status**: `cut -d',' -f10 data/batch_30k.csv | sort | uniq -c`
3. **Verificar FDs**: `lsof -p $(pgrep -f python) | wc -l`
4. **Consultar docs**: 
   - `docs/FIXES-2026-01-13.md` (todos os 28 fixes)
   - `docs/FIXES-ROUND2-SUMMARY.md` (Round 2 específico)

---

**Última atualização**: 2026-01-13  
**Status**: Ready to test 🚀
