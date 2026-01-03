🗂️ DATASET CONTEXT & DELTA STRATEGY - Grimperium v0.2.0

📊 Dataset Overview

**Tamanho:** 52.837 moléculas únicas  
**Colunas:** 59  
**Valores ausentes:** 0 (zero!)  
**Foco:** Termodinâmica computacional de alta precisão (molecular, estado gasoso)

---

### Blocos de Dados

**1. Identificação & Indexação**
- `Unnamed: 0` → índice sequencial
- `smiles` → representação SMILES (crucial para ML)
- `xyz` → coordenadas cartesianas otimizadas (permite descritores geométricos)

**2. Eletrônica Fundamental**
- `charge` → carga total (maioria zero, alguns íons)
- `multiplicity` → spin (maioria singletos, alguns radicais)
- `nheavy` → número de átomos pesados (correlação forte com custo computacional)

**3. Propriedades Termodinâmicas a 298 K** (NÚCLEO CIENTÍFICO)
- `H298_cbs` → Entalpia CBS (referência de alta precisão)
- `H298_b3` → Entalpia B3LYP/DFT (método comum, menos preciso que CBS)
- `S298` → Entropia molar padrão

**4. Parâmetros Auxiliares**
- `A`, `B` → coeficientes para funções termodinâmicas contínuas

**5. Capacidade Calorífica Temperatura-Dependente** (BLOCO RARO & VALIOSO)
- `cp_1` até `cp_45` → 45 pontos de Cp em progresso ordenado de temperatura
- Crescimento suave e monotônico (fisicamente consistente)

---

### Química Representada

✅ **Presente:**
- Hidrocarbonetos alifáticos e aromáticos
- Compostos aromáticos policíclicos
- Álcoois, éteres, aldeídos, cetonas, ácidos carboxílicos
- Ésteres, aminas, amidas, nitrilas, heterociclos
- Compostos sulfurados simples
- Organoclorados/bromados/fluorados

❌ **Ausente:**
- Metais de transição
- Complexos de coordenação
- Sólidos, polímeros, redes cristalinas
- Sais iônicos extensos

---

### Qualidade & Confiabilidade

- ✅ Zero NaNs
- ✅ Comportamento termodinâmico suave e consistente
- ✅ Correlações claras entre tamanho, Cp, entropia
- ✅ Compatibilidade eletrônica (multiplicidade/carga/tipo)
- **Conclusão:** Dataset de altíssima confiabilidade científica

---

## 🎯 Delta-Learning Strategy

### Contexto do Usuário

**Objetivo:** Criar delta (Δ) para o **melhor modelo semiempírico disponível atualmente** para cálculo de entalpia de formação.

### O que é Delta-Learning aqui?

```
Delta = H298_CBS - H298_SEMIEMPIRICAL_BEST

Ou seja:
├─ H298_CBS: referência de altíssima precisão (limite quasi-exato)
├─ H298_SEMIEMPIRICAL_BEST: modelo rápido mas com erro sistemático
└─ Delta: aprender a corrigir a diferença via ML
```

**Vantagem:** Treinar um modelo ML para corrigir erros semiempíricos é muito mais rápido que treinar do zero em CBS.

---

### Fluxo de Dados Revisado

```
CHEMPERIUM + PM7 DATA
    ↓
FEATURE ENGINEERING
├─ Tabular: nheavy, charge, multiplicity
├─ Morgan Fingerprints (256 bits)
└─ RDKit Descriptors: MW, TPSA, LogP
    ↓
DELTA COMPUTATION
└─ delta = H298_CBS - H298_PM7
    ↓
MODEL TRAINING
├─ KernelRidgeRegressor (RBF kernel)
└─ XGBoostRegressor (gradient boosting)
    ↓
ENSEMBLE COMBINATION
├─ Weighted average
└─ Final prediction: H298_CBS ≈ H298_PM7 + delta_ensemble(X)
    ↓
VALIDATION
├─ K-fold CV
├─ Hold-out test set
└─ Metrics: RMSE, MAE, R²
```

---

### 3 Estratégias de Delta-Learning Mapeadas

**OPÇÃO A (Simples - ESCOLHIDA)**
```
Treinar um modelo ML direto no delta: y = H298_CBS - H298_PM7
Predição: H298_CBS ≈ H298_PM7_dado + delta_ML(X)
Vantagem: Simples, direto, interpretável
```

**OPÇÃO B (Ensemble Delta)**
```
Base learner: KRR(X, H298_PM7) → aprende "offset"
Delta learner: XGB(X, delta) → aprende correção
Predição: KRR(X) + XGB(X)
Vantagem: Combina aprendizado de offset + correção
```

**OPÇÃO C (Multi-Delta Comparative)**
```
delta_b3 = H298_CBS - H298_B3
delta_semiemp = H298_CBS - H298_PM7
Treinar ambos, comparar ganho percentual
Vantagem: Explorar dados B3LYP já disponíveis
```

---

## 📌 Próximos Passos Imediatos

**Antes de implementar, confirme:**

1. ✅ PM7 é o semiempírico escolhido
2. ✅ Dados PM7 serão calculados via CREST + MOPAC
3. ✅ Features híbridas (tabular + Morgan + RDKit)
4. ✅ Estratégia A (delta simples)
5. ✅ Validação vs CBS com RMSE/MAE/R²

Todos confirmados! → Vamos para o scaffolding!
