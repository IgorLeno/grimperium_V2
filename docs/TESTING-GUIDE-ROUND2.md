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
wc -l data/batch_test_10.csv  # Deve mostrar 11 linhas (header + 10)
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
# Durante execução, monitorar file descriptors
watch -n 1 'lsof -p $(pgrep -f python) | wc -l'

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

# 2. Monitorar recursos
/usr/bin/time -v python -m grimperium.cli batch-run \
  --csv data/batch_test_100.csv \
  --batch-size 20 \
  > performance_100.log 2>&1

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

### Métricas para Coletar
```python
# Após batch processing
import pandas as pd
df = pd.read_csv('data/batch_test_100.csv')

# Success rate
success_rate = (df['status'] == 'OK').sum() / len(df) * 100
print(f"Success rate: {success_rate:.1f}%")

# Average time per molecule
avg_time = df['total_execution_time'].mean()
print(f"Average time: {avg_time:.1f}s")

# Status distribution
print(df['status'].value_counts())
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

# 2. Rodar em batches (exemplo: 500 moléculas por batch)
# Com checkpointing para poder retomar se falhar
for i in {1..60}; do
  echo "Batch $i/60"
  python -m grimperium.cli batch-run \
    --csv data/batch_30k.csv \
    --batch-size 500 \
    --crest-timeout 30 \
    --mopac-timeout 60 \
    >> production_30k.log 2>&1
  
  # Checkpoint: verificar status counts
  cut -d',' -f10 data/batch_30k.csv | sort | uniq -c
  
  sleep 5
done
```

### Monitoramento Contínuo
```bash
# Terminal 1: Tail logs
tail -f production_30k.log

# Terminal 2: Monitorar recursos
watch -n 10 'ps aux | grep python | grep -v grep'

# Terminal 3: Status dashboard
watch -n 30 'cut -d"," -f10 data/batch_30k.csv | sort | uniq -c'
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
