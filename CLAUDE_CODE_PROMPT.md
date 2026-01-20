# 🤖 PROMPT PARA CLAUDE CODE - IMPLEMENTAÇÃO DAS 3 CORREÇÕES

## 📋 CONTEXTO

Você está trabalhando no projeto **Grimperium V2**, uma framework de Delta Learning para predição de propriedades moleculares. O objetivo é implementar 3 correções críticas identificadas no Phase A test:

- **ISSUE #1:** CSV fields não são calculados (11 colunas vazias)
- **ISSUE #2:** Logging ausente + DtypeWarning
- **ISSUE #3:** Caminhos temporários em /tmp (não-portável)

## 🎯 ENTREGÁVEIS

Você receberá **3 arquivos Python prontos** que resolvem todos esses problemas:

1. **FILE_1_paths.py** - Gerenciamento centralizado de caminhos
2. **FILE_2_logging_enhancements.py** - Logging estruturado
3. **FILE_3_csv_enhancements.py** - Cálculos de deltas e preenchimento de CSV

## 🚀 INSTRUÇÕES DE IMPLEMENTAÇÃO

### FASE 1: PREPARAÇÃO (5 minutos)

```bash
# 1. Crie uma branch para as mudanças
git checkout -b feature/phase-a-fixes

# 2. Backup dos arquivos existentes
git stash

# 3. Copie os 3 arquivos Python para o projeto
cp FILE_1_paths.py src/grimperium/crest_pm7/paths.py
cp FILE_2_logging_enhancements.py src/grimperium/crest_pm7/logging_enhancements.py
cp FILE_3_csv_enhancements.py src/grimperium/crest_pm7/csv_enhancements.py

# 4. Crie a estrutura de diretórios
mkdir -p src/grimperium/crest_pm7/tmp
echo -e "*\n!.gitignore" > src/grimperium/crest_pm7/tmp/.gitignore

# 5. Verifique que os imports funcionam
python -c "from grimperium.crest_pm7.paths import get_molecule_temp_dir; print('✓ FILE_1 OK')"
python -c "from grimperium.crest_pm7.logging_enhancements import setup_batch_logging; print('✓ FILE_2 OK')"
python -c "from grimperium.crest_pm7.csv_enhancements import BatchSettingsCapture; print('✓ FILE_3 OK')"
```

---

### FASE 2: ATUALIZAÇÃO DE ARQUIVOS EXISTENTES

#### 2.1 - src/grimperium/crest_pm7/batch/execution_manager.py

**Mudanças necessárias:**

##### Step 1: Adicione os imports
```python
# No topo do arquivo, após os imports existentes
from grimperium.crest_pm7.paths import (
    get_molecule_temp_dir,
    get_crest_temp_files,
    get_mopac_temp_files,
)
from grimperium.crest_pm7.logging_enhancements import (
    setup_batch_logging,
    log_rdkit_start,
    log_rdkit_done,
    log_crest_start,
    log_crest_done,
    log_mopac_start,
    log_mopac_conformer_done,
    log_mopac_done,
    log_batch_summary,
    suppress_pandas_warnings,
)
from grimperium.crest_pm7.csv_enhancements import (
    BatchSettingsCapture,
    CSVManagerExtensions,
)
```

##### Step 2: Na função `execute_batch()`, adicione no início:
```python
def execute_batch(self, batch_id: str, molecule_ids: List[str]) -> Dict:
    """
    Execute PM7 calculations for a batch of molecules.
    """
    # NOVO: Configurar logging para o batch
    logger = setup_batch_logging(batch_id)
    logger.info(f"Starting batch: {batch_id} ({len(molecule_ids)} molecules)")
    
    # NOVO: Suprimir warnings de pandas
    suppress_pandas_warnings()
    
    # NOVO: Capturar configurações do batch
    batch_settings = BatchSettingsCapture.capture_batch_settings(self.pm7_config)
    logger.debug(f"Batch settings captured: {batch_settings}")
    
    # ... resto do código existente
```

##### Step 3: Para cada molécula, substitua chamadas de /tmp por paths.py:

**ANTES:**
```python
mol_dir = Path(f"/tmp/crest_pm7/{mol_id}")
crest_input = mol_dir / "crest_input.xyz"
crest_output = mol_dir / "crest_conformers.xyz"
```

**DEPOIS:**
```python
# NOVO: Use paths centralizados
crest_files = get_crest_temp_files(batch_id, mol_id)
crest_input = crest_files['input']
crest_output = crest_files['conformers']
```

##### Step 4: Adicione logging em cada etapa da pipeline:

**Ao iniciar RDKit:**
```python
log_rdkit_start(logger, mol_id)
# ... código RDKit ...
log_rdkit_done(logger, mol_id, 
    nrotbonds=rdkit_result['nrotbonds'],
    tpsa=rdkit_result['tpsa'],
    aromatic_rings=rdkit_result['aromatic_rings']
)
```

**Ao iniciar CREST:**
```python
log_crest_start(logger, mol_id)
# ... código CREST ...
log_crest_done(logger, mol_id, 
    num_conformers=len(conformers),
    time_seconds=crest_time
)
```

**Ao iniciar MOPAC:**
```python
log_mopac_start(logger, mol_id, len(conformers))
for conf_idx, conf in enumerate(conformers):
    # ... otimizar conformer ...
    log_mopac_conformer_done(logger, mol_id, 
        conf_idx=conf_idx,
        delta_energy=conf['delta_energy'],
        time_seconds=conf_time
    )
log_mopac_done(logger, mol_id,
    best_conformer_idx=best_idx,
    best_delta_energy=best_delta,
    time_seconds=total_mopac_time
)
```

##### Step 5: Ao final, atualize o CSV com os resultados:

**NOVO, após MOPAC:**
```python
# Atualizar CSV com resultados calculados
success = CSVManagerExtensions.update_molecule_with_mopac_results(
    csv_manager=self.csv_manager,
    mol_id=mol_id,
    h298_cbs=molecule_data['h298_cbs'],
    h298_pm7=mopac_result['h298_pm7'],
    mopac_hof_values=mopac_result['conformer_energies'],
    batch_settings=batch_settings,
    batch_id=batch_id,
)

if success:
    logger.info(f"[{mol_id}] CSV updated successfully")
else:
    logger.error(f"[{mol_id}] Failed to update CSV")
```

##### Step 6: Ao final do batch, adicione sumário:

```python
# Sumário do batch
log_batch_summary(logger, batch_id,
    total=len(molecule_ids),
    success=successful_count,
    failed=failed_count,
    skipped=skipped_count
)

return {
    'batch_id': batch_id,
    'total': len(molecule_ids),
    'successful': successful_count,
    'failed': failed_count,
    'logger': logger,
}
```

---

#### 2.2 - src/grimperium/cli/views/databases_view.py

**Mudança simples - 1 linha apenas:**

Localize a linha (~98):
```python
df = pd.read_csv(csv_path)
```

**MUDE PARA:**
```python
df = pd.read_csv(csv_path, low_memory=False)
```

**Por quê?** Isso suprime o DtypeWarning sobre colunas com tipos mistos.

---

#### 2.3 - src/grimperium/crest_pm7/conformer_generator.py

**Encontre todas as referencias a `/tmp/`:**

```bash
# Procure por padrões
grep -n "/tmp" src/grimperium/crest_pm7/conformer_generator.py
```

**ANTES:**
```python
temp_dir = Path("/tmp/crest_pm7") / mol_id
input_file = temp_dir / "input.xyz"
output_file = temp_dir / "conformers.xyz"
```

**DEPOIS:**
```python
# NOVO: Use paths centralizados
from grimperium.crest_pm7.paths import get_crest_temp_files

crest_files = get_crest_temp_files(batch_id, mol_id)
input_file = crest_files['input']
output_file = crest_files['conformers']
```

**Faça isso para TODAS as referências a /tmp**

---

#### 2.4 - .gitignore (root do projeto)

**Adicione no final:**
```
# Temporary directories for CREST/MOPAC (Phase A fixes)
src/grimperium/crest_pm7/tmp/*/
!src/grimperium/crest_pm7/tmp/.gitignore
```

**Por quê?** Ignora arquivos temporários mas permite o .gitignore existir no git.

---

### FASE 3: TESTES (30 minutos)

#### 3.1 - Testes de import
```bash
# Verifique que tudo importa sem erros
python -c "
from grimperium.crest_pm7.paths import get_molecule_temp_dir
from grimperium.crest_pm7.logging_enhancements import setup_batch_logging
from grimperium.crest_pm7.csv_enhancements import BatchSettingsCapture
print('✓ All imports successful')
"
```

#### 3.2 - Teste unitário dos 3 módulos
```bash
# Execute os testes embutidos
python src/grimperium/crest_pm7/paths.py
python src/grimperium/crest_pm7/logging_enhancements.py
python src/grimperium/crest_pm7/csv_enhancements.py
```

#### 3.3 - Teste com batch real (3 moléculas)
```bash
# Execute a CLI normalmente
python -m grimperium.cli.main

# Na CLI:
# 1. Selecione "Calculate PM7 Values"
# 2. Selecione 3 moléculas
# 3. Observe os logs (devem ter emojis 🧬🔄⚛️)
```

#### 3.4 - Verificações específicas
```bash
# 1. Verificar que CSV foi preenchido
grep "abs_diff" data/molecules_pm7/computed/batch_*.csv
# Resultado esperado: valores numéricos, não vazio

# 2. Verificar que /tmp/crest_pm7 NÃO foi criado
ls -la /tmp/crest_pm7 2>/dev/null
# Resultado esperado: não existe (OK!)

# 3. Verificar que arquivos estão em ./src/crest_pm7/tmp/
ls -la src/grimperium/crest_pm7/tmp/
# Resultado esperado: batch_0001/mol_XXXXX/ (organizado!)

# 4. Verificar que não há warnings
python -m grimperium.cli.main 2>&1 | grep -i warning
# Resultado esperado: nada (OK!)
```

---

### FASE 4: VALIDAÇÃO FINAL (15 minutos)

```bash
# 1. Rodar testes do projeto
pytest tests/ -v --cov

# 2. Verificar type hints
mypy src/grimperium/crest_pm7/ --strict

# 3. Verificar linting
ruff check src/grimperium/crest_pm7/

# 4. Verificar formatação
black --check src/grimperium/crest_pm7/
```

---

### FASE 5: COMMIT E MERGE

```bash
# 1. Revisar mudanças
git diff

# 2. Preparar commit
git add src/grimperium/crest_pm7/paths.py
git add src/grimperium/crest_pm7/logging_enhancements.py
git add src/grimperium/crest_pm7/csv_enhancements.py
git add src/grimperium/crest_pm7/batch/execution_manager.py
git add src/grimperium/cli/views/databases_view.py
git add src/grimperium/crest_pm7/conformer_generator.py
git add .gitignore
git add src/grimperium/crest_pm7/tmp/.gitignore

# 3. Commit com mensagem clara
git commit -m "fix: Phase A - implement 3 critical fixes

- ISSUE #1: Auto-calculate 11 CSV fields (FILE_3_csv_enhancements.py)
- ISSUE #2: Add structured logging + suppress DtypeWarning (FILE_2_logging_enhancements.py)
- ISSUE #3: Move temp paths from /tmp to project-local (FILE_1_paths.py)

Changes:
- New: 3 production-ready modules (993 LOC)
- Updated: execution_manager, databases_view, conformer_generator, .gitignore
- Tests: All passing, no warnings, CSV complete

Verified:
✓ CSV fields filled (all 11 columns)
✓ Detailed logging per tool (RDKit, CREST, MOPAC)
✓ Paths project-local (./src/crest_pm7/tmp/)
✓ No DtypeWarning or FutureWarning
✓ Type hints complete
✓ All imports working
"

# 4. Push para review
git push origin feature/phase-a-fixes

# 5. Criar Pull Request no GitHub
# - Title: "Phase A Fixes - 3 Critical Issues Resolved"
# - Description: Cole o commit message acima
```

---

## ✅ CHECKLIST DE SUCESSO

Após implementação, verifique:

- [ ] **CSV Completo**
  - [ ] Coluna `abs_diff` preenchida
  - [ ] Coluna `delta_1` preenchida
  - [ ] Coluna `delta_2` preenchida
  - [ ] Coluna `delta_3` preenchida
  - [ ] Coluna `conformer_selected` preenchida
  - [ ] Configurações (v3, qm, nci, etc.) preenchidas

- [ ] **Logging Detalhado**
  - [ ] RDKit logs aparecem com 🧬
  - [ ] CREST logs aparecem com 🔄
  - [ ] MOPAC logs aparecem com ⚛️
  - [ ] Sem DtypeWarning
  - [ ] Sem FutureWarning
  - [ ] Resumo de batch ao final

- [ ] **Caminhos Corretos**
  - [ ] Arquivos em `./src/crest_pm7/tmp/batch_XXX/mol_XXXXX/`
  - [ ] `/tmp/crest_pm7` NÃO criado
  - [ ] `.gitignore` previne commit de arquivos temporários
  - [ ] Fácil navegar pela estrutura

- [ ] **Qualidade de Código**
  - [ ] Todos os imports funcionam
  - [ ] Sem erros de type hint
  - [ ] Sem erros de linting
  - [ ] Formatação correta (black)
  - [ ] Testes passando

---

## 📚 RECURSOS ADICIONAIS

Se precisar de mais contexto durante a implementação:

1. **FILE_1_paths.py** - Leia as docstrings para exemplos de uso
2. **FILE_2_logging_enhancements.py** - Veja exemplos de output esperado
3. **FILE_3_csv_enhancements.py** - Entenda o conceito de "deltas de energia"
4. **INTEGRATION_GUIDE.md** - Guia detalhado passo-a-passo
5. **QUICK_REFERENCE.md** - Padrões e snippets para copiar/colar

---

## ⏱️ TEMPO ESTIMADO

- Preparação: 5 min
- Atualizações: 2-3 horas
- Testes: 30-45 min
- Commit: 15 min
- **TOTAL: ~6 horas**

---

## 🚀 BORA LÁ!

Você tem tudo que precisa. Os 3 arquivos estão prontos, as instruções são detalhadas e específicas.

**Próximo passo:** Comece com FASE 1 - Preparação (5 minutos)

Boa implementação! 🎉
