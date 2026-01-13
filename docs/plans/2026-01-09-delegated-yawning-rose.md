# CREST-PM7 Pipeline v0.3.2 — Plano Final de Implementação

**Status:** 100% especificado, pronto para implementação
**Versão:** 0.2.0 → 0.3.2
**Data:** 2026-01-09

---

## Quick Start (Próxima Sessão)

**Para continuar a implementação:**

```bash
# 1. Comando inicial na próxima sessão:
# "Continuar implementação do CREST-PM7 Pipeline v0.3.2 seguindo o plano em /home/igor/.claude/plans/delegated-yawning-rose.md"

# 2. Primeiro módulo a implementar:
# src/grimperium/crest_pm7/config.py (enums + PM7Config)

# 3. Ordem de implementação (18 itens total):
# Fase 1: config.py → logging_utils.py → validation.py
# Fase 2: energy_extractor.py → timeout_predictor.py → conformer_selector.py → threshold_monitor.py
# Fase 3: conformer_generator.py → mopac_optimizer.py
# Fase 4: molecule_processor.py → result_evaluator.py → pipeline.py → __init__.py
# Fase 5: diretórios → baselines → scripts → CHANGELOG.md
```

**Contexto do projeto:**
- Projeto: grimperium (química computacional)
- Local: `/home/igor/Projetos/grimperium`
- Novo módulo: `src/grimperium/crest_pm7/` (não existe ainda)
- Dependências: RDKit, sklearn, CREST (externo), MOPAC (externo), Open Babel (externo)

**Decisões já tomadas (não renegociar):**
- 5 regex patterns para HOF (HIGH→HIGH→MEDIUM→MEDIUM→LOW)
- HOF bounds: [-500, +500] kcal/mol
- Timeout: Huber regression, recalibrar a cada 50 moléculas
- Baseline: ±2.5 kcal/mol tolerância absoluta
- XYZ→SDF: Open Babel primeiro, RDKit fallback (SMILES obrigatório)
- Error handling: camadas não propagam exceções

---

## Objetivo

Pipeline robusto CREST + PM7 para processar 800 moléculas com:
- Geração de conformadores (CREST)
- Otimização PM7 (MOPAC)
- Extração robusta de HOF (5 regex)
- Monitoramento de qualidade (5 padrões)
- Baseline validation

---

## Arquitetura Final (12 Módulos)

```
src/grimperium/crest_pm7/
├── __init__.py                 # Exports (~35 linhas)
├── config.py                   # PM7Config + 6 enums (~120 linhas)
├── validation.py               # Environment validation (~70 linhas)
├── logging_utils.py            # Structured JSONL logging (~90 linhas)
├── energy_extractor.py         # 5 regex patterns + validate_hof (~150 linhas)
├── timeout_predictor.py        # Huber regression (~200 linhas)
├── conformer_selector.py       # get_num_conformers + ΔE (~100 linhas)
├── threshold_monitor.py        # 5 detection patterns (~250 linhas)
├── conformer_generator.py      # CREST + XYZ→SDF (~280 linhas)
├── mopac_optimizer.py          # MOPAC + robust parsing (~250 linhas)
├── molecule_processor.py       # PM7Result + ConformerData + MoleculeProcessor (~350 linhas)
├── result_evaluator.py         # ResultEvaluator (~120 linhas)
└── pipeline.py                 # CRESTPM7Pipeline orchestrator (~100 linhas)

data/molecules_pm7/
├── raw/
├── computed/
│   └── logs/
├── testing/
│   └── baselines/
│       ├── phase_a_molecules.csv
│       └── phase_a_expected.json
├── training_sets/
├── models/
└── analysis/

scripts/
├── phase_a_quick_test.py
└── utils/
    ├── baseline_validator.py
    └── baseline_generator.py
```

---

## Especificação 1: Energy Extractor

### 5 Padrões Regex (em ordem de tentativa)

```python
PATTERNS = [
    # Pattern 1: FINAL HEAT OF FORMATION (HIGH)
    # Linha: "FINAL HEAT OF FORMATION =      -17.93603 KCAL/MOL =      -75.04 KJ/MOL"
    {
        "name": "FINAL_HEAT_OF_FORMATION",
        "regex": r"FINAL\s+HEAT\s+OF\s+FORMATION\s*=\s*([-+]?\d+\.?\d*)\s*KCAL/MOL",
        "confidence": HOFConfidence.HIGH,
    },

    # Pattern 2: HEAT OF FORMATION sem FINAL (HIGH)
    # Linha: "HEAT OF FORMATION       =        -17.93603 KCAL/MOL"
    {
        "name": "HEAT_OF_FORMATION",
        "regex": r"HEAT\s+OF\s+FORMATION\s*=\s*([-+]?\d+\.?\d*)\s*KCAL/MOL",
        "confidence": HOFConfidence.HIGH,
    },

    # Pattern 3: Formato compacto (MEDIUM)
    # Linha: "HEAT OF FORMATION=-17.93603"
    {
        "name": "HEAT_OF_FORMATION_COMPACT",
        "regex": r"HEAT\s+OF\s+FORMATION\s*=\s*([-+]?\d+\.?\d+)",
        "confidence": HOFConfidence.MEDIUM,
    },

    # Pattern 4: Thermal energy 298K (MEDIUM)
    # Linha: "THERMAL ENERGY AT 298K = -12.345"
    {
        "name": "THERMAL_ENERGY_298K",
        "regex": r"THERMAL\s+ENERGY\s+AT\s+298K?\s*=\s*([-+]?\d+\.?\d*)",
        "confidence": HOFConfidence.MEDIUM,
    },

    # Pattern 5: Fallback - última ocorrência (LOW)
    # Qualquer "... -17.93603 KCAL/MOL ..."
    {
        "name": "FALLBACK_LAST_KCAL",
        "regex": r"([-+]?\d+\.?\d+)\s+KCAL/MOL",
        "confidence": HOFConfidence.LOW,
        "use_finditer_last": True,  # Pegar última match
    },
]
```

### Bounds e Validação

```python
HOF_MIN_BOUND = -500.0  # kcal/mol
HOF_MAX_BOUND = +500.0  # kcal/mol

def validate_hof(hof: float, nheavy: int) -> tuple[bool, str]:
    # Hard fail se fora de [-500, +500]
    if hof < HOF_MIN_BOUND or hof > HOF_MAX_BOUND:
        return False, f"HOF {hof:.2f} outside bounds [{HOF_MIN_BOUND}, {HOF_MAX_BOUND}]"

    # Soft warning se muito longe de esperado (~-10*nheavy)
    expected = -10.0 * nheavy
    if abs(hof - expected) > 300.0:
        LOG.warning(f"HOF {hof:.2f} deviates from expected ~{expected:.0f}")

    return True, "OK"
```

### Retornos

```python
# Sucesso: (hof=-17.93, method="FINAL_HEAT_OF_FORMATION", confidence=HIGH)
# Falha parsing: (None, None, None)
# Falha bounds: (None, "FINAL_HEAT_OF_FORMATION", HIGH)  # preserva diagnóstico
```

---

## Especificação 2: Timeout Predictor

### Modelo e Features

- **Modelo:** HuberRegressor (sklearn)
- **Feature única:** nheavy (número de átomos pesados)
- **Sem transformação:** tempo em segundos, linear

### Constantes

```python
DEFAULT_TIMEOUT = 300.0     # 5 minutos
MIN_SAMPLES_FOR_FIT = 20    # mínimo para treinar
TIMEOUT_MIN = 60.0          # mínimo absoluto
TIMEOUT_MAX = 3600.0        # máximo (1 hora)
```

### Lógica por Fase

```python
def predict(nheavy: int, num_conformers: int) -> tuple[float, TimeoutConfidence]:
    if n_samples < 20:
        # Phase A: heurística pura
        timeout = (120 + 5*nheavy) * (1 + 0.2*(num_conformers-1))
        return clamp(timeout), TimeoutConfidence.LOW

    elif n_samples < 50:
        # Phase B: modelo + 50% margem
        pred = model.predict([[nheavy]])[0]
        timeout = pred * num_conformers * 1.5
        return clamp(timeout), TimeoutConfidence.MEDIUM

    else:
        # Phase C+: modelo + 30% margem
        pred = model.predict([[nheavy]])[0]
        timeout = pred * num_conformers * 1.3
        return clamp(timeout), TimeoutConfidence.HIGH
```

### Recalibração

- **Quando:** a cada 50 moléculas (config.timeout_predictor_recalibrate_interval)
- **Como:** fit do zero com todos dados acumulados
- **Persistência:** `data/molecules_pm7/models/timeout_predictor.pkl`
- **Carregamento:** manual via load(), não automático

---

## Especificação 3: Baseline Validation

### Tolerância

```python
TOLERANCE_ABSOLUTE = 2.5  # kcal/mol

# Fórmula: hof_min <= hof_result <= hof_max
# Onde: hof_min = hof_value - 2.5, hof_max = hof_value + 2.5
```

### Critérios Molecule-Level Pass

Uma molécula passa se TODAS as condições são verdadeiras:
1. `result.success == True`
2. `result.most_stable_hof is not None`
3. `hof_min <= result.most_stable_hof <= hof_max`
4. `result.quality_grade in [A, B]`

### Critérios Phase A Pass (Global)

```python
PHASE_A_CRITERIA = {
    'min_success_rate': 1.0,           # 100% (3/3)
    'min_hof_extraction_rate': 1.0,    # 100% (3/3)
    'min_baseline_pass_rate': 1.0,     # 100% (3/3)
    'min_grade_ab_rate': 0.67,         # 67% (2/3)
    'require_zero_crashes': True,
}
```

### Exit Codes

- `exit 0`: Todas as condições satisfeitas
- `exit 1`: Qualquer condição falha

---

## Especificação 4: XYZ → SDF Fallback (RDKit)

### Estratégia

1. Open Babel primeiro: `obabel input.xyz -O output.sdf -h`
2. Se falhar, RDKit fallback com SMILES obrigatório

### RDKit Fallback (Estratégia 2A)

```python
def _convert_xyz_to_sdf_rdkit(xyz_file, sdf_file, smiles, mol_id, conf_index):
    # 1. Criar RDKit Mol do SMILES
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # 2. Parsear XYZ para coordenadas
    coords, n_atoms_xyz = parse_xyz(xyz_file)

    # 3. Verificar contagem de átomos
    if mol.GetNumAtoms() != n_atoms_xyz:
        return False  # FALHA

    # 4. Criar conformer com coordenadas do XYZ
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(x, y, z))
    mol.AddConformer(conf)

    # 5. Exportar SDF
    writer = Chem.SDWriter(str(sdf_file))
    writer.write(mol)
    writer.close()
    return True
```

### Decisão: SMILES Obrigatório

- Se SMILES não disponível: **FALHAR** (retornar False)
- Não tentar inferir conectividade apenas de coordenadas
- Razão: inferência é error-prone, resultados seriam suspeitos

---

## Especificação 5: Error Handling por Camada

### Camada Baixa (Utilities)

**Arquivos:** energy_extractor, timeout_predictor, conformer_selector, logging_utils

- **Exceções:** NÃO levanta
- **Retorno:** tuplas com status ou None
- **Logging:** DEBUG para detalhes, WARNING para problemas

### Camada Média (Wrappers)

**Arquivos:** conformer_generator, mopac_optimizer, validation

- **Exceções:** NÃO levanta (captura internamente)
- **Retorno:** objetos Result com status (SUCCESS/FAILED/etc.)
- **Logging:** INFO início/fim, WARNING problemas, ERROR falhas

### Camada Alta (Orquestração)

**Arquivos:** molecule_processor, result_evaluator, threshold_monitor, pipeline

- **Exceções:** NÃO levanta
- **Retorno:** sempre PM7Result/EvaluationResult válido
- **Comportamento:** converte falhas em result.success=False
- **Logging:** INFO progresso, WARNING problemas, ERROR falhas, CRITICAL pausar

---

## Especificação 6: Enums e Dataclasses

### Enums (config.py)

```python
class QualityGrade(str, Enum):
    A = "A"           # Confiável
    B = "B"           # Aceitável
    C = "C"           # Suspeito
    FAILED = "FAILED" # Erro crítico

class HOFConfidence(str, Enum):
    HIGH = "HIGH"     # Padrão 1-2
    MEDIUM = "MEDIUM" # Padrão 3-4
    LOW = "LOW"       # Padrão 5 (fallback)

class TimeoutConfidence(str, Enum):
    LOW = "LOW"       # < 20 samples
    MEDIUM = "MEDIUM" # 20-50 samples
    HIGH = "HIGH"     # 50+ samples

class CRESTStatus(str, Enum):
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"
    NOT_ATTEMPTED = "NOT_ATTEMPTED"

class MOPACStatus(str, Enum):
    SUCCESS = "SUCCESS"
    TIMEOUT = "TIMEOUT"
    SCF_FAILED = "SCF_FAILED"
    GEOMETRY_ERROR = "GEOMETRY_ERROR"
    NOT_ATTEMPTED = "NOT_ATTEMPTED"

class AlertLevel(str, Enum):
    INFO = "INFO"
    WARNING = "WARNING"
    CRITICAL = "CRITICAL"
```

### PM7Config (config.py)

```python
@dataclass
class PM7Config:
    # Phase control
    phase: str = "A"  # A | B | C | production

    # Executáveis
    crest_executable: str = "crest"
    mopac_executable: str = "mopac"

    # Caminhos
    temp_dir: Path = Path("/tmp/crest_pm7")
    output_dir: Path = Path("data/molecules_pm7/computed")

    # CREST
    max_conformers: int = 10
    energy_window: float = 6.0  # kcal/mol
    crest_timeout: float = 300.0  # segundos

    # MOPAC
    mopac_timeout_base: float = 120.0
    mopac_timeout_margin: float = 1.3

    # Conformer selection (calibrado em Phase C)
    nrotbonds_threshold_rigid_to_medium: int = 1
    nrotbonds_threshold_medium_to_flexible: int = 4

    # Timeout predictor
    timeout_predictor_recalibrate_interval: int = 50

    # Monitoring
    success_rate_warning: float = 0.80
    success_rate_critical: float = 0.70
    monitor_window_size: int = 50
```

### ConformerData (molecule_processor.py)

```python
@dataclass
class ConformerData:
    index: int
    mol_id: str

    # CREST
    crest_status: CRESTStatus = CRESTStatus.NOT_ATTEMPTED
    crest_geometry_file: Optional[Path] = None
    crest_error_message: Optional[str] = None

    # MOPAC
    mopac_status: MOPACStatus = MOPACStatus.NOT_ATTEMPTED
    mopac_output_file: Optional[Path] = None
    mopac_execution_time: float = 0.0
    mopac_timeout_used: Optional[float] = None
    mopac_error_message: Optional[str] = None

    # Energy
    energy_hof: Optional[float] = None
    hof_confidence: HOFConfidence = HOFConfidence.LOW
    hof_extraction_method: Optional[str] = None
    hof_extraction_successful: bool = False

    @property
    def is_successful(self) -> bool:
        return (self.mopac_status == MOPACStatus.SUCCESS and
                self.hof_extraction_successful and
                self.energy_hof is not None)
```

### PM7Result (molecule_processor.py)

```python
@dataclass
class PM7Result:
    # Identificação
    mol_id: str
    smiles: str
    timestamp: datetime = field(default_factory=datetime.now)
    phase: str = "A"

    # Metadados
    nheavy: Optional[int] = None
    nrotbonds: Optional[int] = None
    tpsa: Optional[float] = None
    aromatic_rings: Optional[int] = None
    has_heteroatoms: Optional[bool] = None

    # CREST
    crest_status: CRESTStatus = CRESTStatus.NOT_ATTEMPTED
    crest_conformers_generated: int = 0
    crest_time: Optional[float] = None
    crest_error: Optional[str] = None

    # MOPAC
    conformers: List[ConformerData] = field(default_factory=list)
    num_conformers_selected: Optional[int] = None
    total_execution_time: Optional[float] = None

    # Energy differences
    delta_e_12: Optional[float] = None
    delta_e_13: Optional[float] = None
    delta_e_15: Optional[float] = None

    # Timeout
    timeout_predicted: Optional[float] = None
    timeout_confidence: Optional[TimeoutConfidence] = None

    # Decisions (auditoria)
    decisions: List[str] = field(default_factory=list)

    # Quality
    quality_grade: QualityGrade = QualityGrade.FAILED
    issues: List[str] = field(default_factory=list)

    # Final
    success: bool = False
    error_message: Optional[str] = None

    @property
    def most_stable_hof(self) -> Optional[float]:
        successful = [c for c in self.conformers if c.is_successful]
        return min(c.energy_hof for c in successful) if successful else None

    def to_dict(self) -> dict:
        # Serialização JSON-safe (Path → str)
        ...
```

---

## Especificação 7: Baseline Data

### phase_a_molecules.csv

```csv
mol_id,smiles,nheavy,notes
methane,C,1,tiny molecule - no flexibility
ethane,CC,2,small molecule - 1 rotatable bond
benzene,c1ccccc1,6,aromatic - rigid structure
```

### phase_a_expected.json

```json
{
  "generation_date": "2026-01-09",
  "source": "MOPAC PM7 literature values",
  "tolerance_kcal_mol": 2.5,

  "molecules": {
    "methane": {
      "hof_value": -17.93,
      "hof_min": -20.43,
      "hof_max": -15.43,
      "quality_grade_expected": "A",
      "success_expected": true
    },
    "ethane": {
      "hof_value": -20.24,
      "hof_min": -22.74,
      "hof_max": -17.74,
      "quality_grade_expected": "A",
      "success_expected": true
    },
    "benzene": {
      "hof_value": -25.68,
      "hof_min": -28.18,
      "hof_max": -23.18,
      "quality_grade_expected": "A",
      "success_expected": true
    }
  },

  "phase_a_success_criteria": {
    "min_success_rate": 1.0,
    "min_hof_extraction_rate": 1.0,
    "min_baseline_pass_rate": 1.0,
    "min_grade_ab_rate": 0.67,
    "require_zero_crashes": true
  }
}
```

---

## Especificação 8: Logging

### Estrutura

```python
# JSONL (para análise):
# data/molecules_pm7/computed/logs/phase_A_20260109_120000.jsonl

# Texto (para debugging):
# data/molecules_pm7/computed/logs/phase_A_20260109_120000.log

# Console: stdout em tempo real
```

### Formato JSONL

```json
{"timestamp": "2026-01-09T12:00:00.123", "level": "INFO", "logger": "grimperium.crest_pm7.pipeline", "message": "mol_start", "mol_id": "methane", "smiles": "C"}
{"timestamp": "2026-01-09T12:00:45.456", "level": "INFO", "logger": "grimperium.crest_pm7.pipeline", "message": "mol_complete", "mol_id": "methane", "grade": "A", "hof": -17.93}
```

### Política de Retenção (v0.4.0+)

- Max file size: 50 MB
- Backup count: 5
- Compression: gzip
- Retention: 30 dias

---

## Ordem de Implementação

### Fase 1: Core (sem dependências)
1. `config.py` - enums + PM7Config
2. `logging_utils.py` - StructuredLogHandler
3. `validation.py` - validate_environment

### Fase 2: Processing
4. `energy_extractor.py` - 5 regex + validate_hof
5. `timeout_predictor.py` - Huber regression
6. `conformer_selector.py` - get_num_conformers + ΔE
7. `threshold_monitor.py` - 5 patterns + Alert

### Fase 3: Wrappers
8. `conformer_generator.py` - CREST + XYZ→SDF
9. `mopac_optimizer.py` - MOPAC + robust parsing

### Fase 4: Orquestração
10. `molecule_processor.py` - PM7Result + ConformerData + processor
11. `result_evaluator.py` - ResultEvaluator
12. `pipeline.py` - CRESTPM7Pipeline
13. `__init__.py` - exports

### Fase 5: Teste
14. Criar diretórios
15. Baseline data (CSV + JSON)
16. `baseline_validator.py`
17. `phase_a_quick_test.py`
18. `CHANGELOG.md`

---

## Verificação Final

```bash
# 1. Estrutura
ls -la src/grimperium/crest_pm7/
ls -la data/molecules_pm7/testing/baselines/

# 2. Imports
python -c "from grimperium.crest_pm7 import CRESTPM7Pipeline, PM7Config; print('OK')"

# 3. Environment
python -c "from grimperium.crest_pm7 import validate_environment, PM7Config; print(validate_environment(PM7Config()))"

# 4. Phase A
python scripts/phase_a_quick_test.py

# 5. Resultados
cat data/molecules_pm7/computed/phase_a_results.json

# 6. Logs
head -20 data/molecules_pm7/computed/logs/phase_A_*.jsonl

# 7. Testes existentes
pytest tests/ -x
```

---

## Resumo de Decisões Críticas

| Área | Decisão |
|------|---------|
| HOF Regex | 5 padrões: HIGH→HIGH→MEDIUM→MEDIUM→LOW |
| HOF Bounds | [-500, +500] kcal/mol |
| Timeout inicial | 300s heurística (Phase A) |
| Timeout modelo | Huber, fit com 20+ samples, recalibrar cada 50 |
| Baseline tolerância | ±2.5 kcal/mol absoluto |
| Phase A pass | 100% success, 100% HOF, 100% baseline, 67% A/B |
| XYZ→SDF | Open Babel + RDKit (SMILES obrigatório) |
| Error handling | Camadas não propagam exceções |
| Timeline | 20-24 horas |

---

**PRÓXIMO PASSO:** Implementar módulos na ordem especificada, começando por `config.py`.
