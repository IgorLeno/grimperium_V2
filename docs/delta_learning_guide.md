# Delta-Learning Guide

This guide explains the delta-learning approach used in Grimperium for improving semiempirical thermodynamic predictions.

## What is Delta-Learning?

Delta-learning is a machine learning strategy where we train models to predict the **difference (delta)** between two levels of theory, rather than predicting absolute values directly.

```
δ = y_high_accuracy - y_low_accuracy
```

In Grimperium:
```
δ = H298_CBS - H298_PM7
```

Where:
- **H298_CBS**: High-accuracy CBS (Complete Basis Set) enthalpy
- **H298_PM7**: Fast semiempirical PM7 enthalpy
- **δ**: The correction that ML learns to predict

## Why Delta-Learning?

### 1. Physical Motivation

Semiempirical methods like PM7 capture most of the physics:

- Correct general trends
- Reasonable absolute values
- ~70-80% of the variance explained

The remaining ~20-30% error is systematic and learnable.

### 2. Statistical Advantages

Training on deltas instead of absolute values:

| Aspect | Direct Learning | Delta Learning |
|--------|-----------------|----------------|
| Target range | -200 to +50 kcal/mol | -10 to +10 kcal/mol |
| Target variance | High | Low |
| Outlier sensitivity | High | Low |
| Extrapolation | Risky | Safer |

### 3. Interpretability

Delta predictions have clear physical meaning:
- δ > 0: PM7 underestimates binding
- δ < 0: PM7 overestimates binding
- |δ| large: PM7 struggles with this chemistry

## The Grimperium Approach

### Strategy A: Simple Delta (Implemented)

```
1. Compute delta: δ = H298_CBS - H298_PM7
2. Train model: f(X) → δ
3. Predict: H298_CBS ≈ H298_PM7 + f(X)
```

**Advantages:**
- Simple and interpretable
- Single model to train
- Easy to debug

### Strategy B: Ensemble Delta (Future)

```
1. Base learner: KRR(X) → offset
2. Delta learner: XGB(X) → correction
3. Combine: prediction = KRR(X) + XGB(X)
```

### Strategy C: Multi-Delta Comparative (Future)

```
1. δ_PM7 = H298_CBS - H298_PM7
2. δ_B3 = H298_CBS - H298_B3
3. Compare which delta is easier to learn
```

## Mathematical Framework

### Delta Computation

Given:
- CBS reference: y_cbs ∈ ℝⁿ
- PM7 values: y_pm7 ∈ ℝⁿ
- Features: X ∈ ℝⁿˣᵈ

Delta:
```
δ = y_cbs - y_pm7
```

### Model Training

Minimize:
```
L(θ) = Σᵢ (f_θ(xᵢ) - δᵢ)²
```

Where f_θ is the ML model (KRR or XGBoost).

### Prediction

For new molecule with features x_new and PM7 value y_pm7_new:
```
ŷ_cbs = y_pm7_new + f_θ(x_new)
```

## Typical Delta Distribution

Based on PM7 vs CBS for organic molecules:

```
         │
         │      ┌───┐
         │      │   │
         │     ┌┘   └┐
         │    ┌┘     └┐
         │   ┌┘       └┐
         │  ┌┘         └┐
         │ ┌┘           └┐
    ─────┴─┴─────────────┴─────
        -10  -5   0   +5  +10
              δ (kcal/mol)
```

- **Mean**: ~-3 kcal/mol (PM7 less negative than CBS)
- **Std**: ~3 kcal/mol
- **Range**: typically -10 to +5 kcal/mol

## When Delta-Learning Works Best

### Good Cases

1. **Systematic errors**: PM7 consistently over/underestimates
2. **Size-dependent errors**: Error scales with molecule size
3. **Functional group errors**: Specific groups have consistent bias

### Challenging Cases

1. **Random errors**: No pattern to learn
2. **Outliers**: Very unusual chemistry
3. **Edge cases**: Far from training distribution

## Validation Strategy

### Cross-Validation

```python
from sklearn.model_selection import KFold

kf = KFold(n_splits=5, shuffle=True, random_state=42)
for train_idx, val_idx in kf.split(X):
    # Train on train_idx, evaluate on val_idx
    pass
```

### Hold-Out Test Set

```python
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, delta, test_size=0.2, random_state=42
)
```

### Metrics

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| RMSE | √(Σ(y-ŷ)²/n) | Standard error magnitude |
| MAE | Σ|y-ŷ|/n | Average absolute error |
| R² | 1 - SS_res/SS_tot | Variance explained |

### Success Criteria

```
RMSE(PM7+δ_ML vs CBS) << RMSE(PM7 vs CBS)
```

Typical improvement:
- PM7 raw: 4-6 kcal/mol RMSE
- PM7 + δ_ML: 0.5-2 kcal/mol RMSE

## Code Example

```python
import numpy as np
from grimperium.core import DeltaLearner
from grimperium.core.metrics import rmse, mae

# Load data
learner = DeltaLearner()
learner.load_data("chemperium.csv", "pm7.csv")

# Compute features
learner.compute_features()

# Train
learner.train()

# Evaluate
metrics = learner.evaluate()

# Compare with baseline
h298_pm7 = test_data["H298_pm7"]
h298_cbs = test_data["H298_cbs"]

# Baseline error (PM7 raw)
baseline_rmse = rmse(h298_cbs, h298_pm7)

# Delta-corrected error
delta_pred = learner.predict(X_test, return_delta=True)
corrected = h298_pm7 + delta_pred
corrected_rmse = rmse(h298_cbs, corrected)

print(f"PM7 raw RMSE: {baseline_rmse:.2f} kcal/mol")
print(f"PM7+δ RMSE: {corrected_rmse:.2f} kcal/mol")
print(f"Improvement: {baseline_rmse - corrected_rmse:.2f} kcal/mol")
```

## Best Practices

### 1. Data Quality

- Ensure CBS values are high-quality
- Validate PM7 calculations converged
- Check for outliers before training

### 2. Feature Selection

- Include molecular size (nheavy)
- Use fingerprints for substructure patterns
- Add physicochemical descriptors

### 3. Model Selection

- Start with simple KRR (smooth corrections)
- Add XGBoost for complex patterns
- Use ensemble for robustness

### 4. Hyperparameter Tuning

- Tune KRR alpha via cross-validation
- Use early stopping for XGBoost
- Optimize ensemble weights on validation set

## Common Pitfalls

### 1. Data Leakage

Don't use CBS values as features!

### 2. Overfitting

- Use regularization (alpha in KRR, early stopping in XGB)
- Monitor validation metrics
- Keep model complexity in check

### 3. Extrapolation

- Test on chemically diverse molecules
- Be cautious with very different chemistry
- Consider uncertainty estimates

## References

1. Ramakrishnan et al. (2015). "Big Data Meets Quantum Chemistry"
2. von Lilienfeld et al. (2020). "Quantum Machine Learning"
3. Stewart (2013). "PM7 Optimization for Organometallic Compounds"
