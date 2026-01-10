# Grimperium Data Module

O módulo `data` fornece funcionalidades para carregamento, processamento e fusão de dados químicos e moleculares.

## Componentes

### Loader

Sistema de carregamento de dados de múltiplas fontes.

```python
from grimperium.data.loader import DataLoader

# Carregar dados
loader = DataLoader()
data = loader.load_from_file("dataset.csv")

# Carregar dados com configuração customizada
data = loader.load_from_file(
    "dataset.csv",
    format="csv",
    columns=["smiles", "property"]
)
```

**Formatos suportados:**
- CSV
- JSON
- HDF5
- Pickle
- Custom formats via plugins

### Semiempirical

Interface para dados de métodos semiempíricos (PM6, PM7, etc.).

```python
from grimperium.data.semiempirical import SemiempiricalLoader

# Carregar dados semiempíricos
se_loader = SemiempiricalLoader()
se_data = se_loader.load_pm6_data("pm6_results.json")

# Processar propriedades
properties = se_loader.extract_properties(se_data)
```

**Propriedades suportadas:**
- Energias de formação
- Geometrias otimizadas
- Propriedades eletrônicas
- Cargas atômicas

### Fusion

Sistema de fusão de dados de múltiplas fontes.

```python
from grimperium.data.fusion import DataFusion

# Criar fusão de dados
fusion = DataFusion()

# Adicionar fontes de dados
fusion.add_source("dft", dft_data, weight=0.6)
fusion.add_source("semiempirical", se_data, weight=0.4)

# Fundir dados
merged_data = fusion.merge()
```

**Estratégias de fusão:**
- Weighted averaging
- Consensus-based
- Hierarchical fusion
- Custom fusion strategies

## Arquivos

- `loader.py` - Sistema de carregamento de dados
- `semiempirical.py` - Interface para dados semiempíricos
- `fusion.py` - Sistema de fusão de dados

## Workflow Típico

```python
from grimperium.data import DataLoader, SemiempiricalLoader, DataFusion

# 1. Carregar dados DFT
loader = DataLoader()
dft_data = loader.load_from_file("dft_dataset.csv")

# 2. Carregar dados semiempíricos
se_loader = SemiempiricalLoader()
se_data = se_loader.load_pm6_data("pm6_results.json")

# 3. Fundir dados
fusion = DataFusion()
fusion.add_source("dft", dft_data)
fusion.add_source("pm6", se_data)
final_data = fusion.merge()

# 4. Usar para treinamento
from grimperium.core import DeltaLearner
learner = DeltaLearner()
learner.fit(final_data)
```

## Status

✅ **Production Ready**
- Múltiplos formatos suportados
- Validação de dados integrada
- Testes extensivos

## Dependências

- Pandas
- NumPy
- h5py (para HDF5)
- RDKit (para processamento molecular)

## Ver também

- [Core Module](README_core.md) - Algoritmos de aprendizado
- [Models Module](README_models.md) - Modelos de ML
- [Documentation](../docs/build/html/grimperium.data.html) - API completa
