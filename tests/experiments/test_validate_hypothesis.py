"""
Hypothesis Validation Experiment.

Core Hypothesis:
    Delta-learning (training on CBS-PM7 correction) produces better predictions
    than direct learning (training directly on CBS values).

Experimental Design:
    - Model Delta: DeltaLearner (learns y_delta = y_cbs - y_pm7)
    - Model Direct: DeltaLearningEnsemble (learns y_cbs directly)
    - Same architecture: w_krr=0.5, w_xgb=0.5
    - Same features: 10D enriched features
    - Same split: 80/20, random_state=42

Decision Gate:
    1. rmse_delta < rmse_direct (hypothesis validation)
    2. rmse_delta < 20.0 kcal/mol (quality threshold)

Expected Results:
    - rmse_delta: 8-12 kcal/mol
    - rmse_direct: 15-30 kcal/mol
    - Test PASSES if both criteria met, FAILS otherwise

Usage:
    pytest tests/experiments/test_validate_hypothesis.py -v -s
"""

import numpy as np
import pytest
from sklearn.model_selection import train_test_split

from grimperium.core.delta_learning import DeltaLearner
from grimperium.models.delta_ensemble import DeltaLearningEnsemble
from grimperium.core.metrics import compute_all_metrics


def test_decision_gate_delta_vs_direct(real_data_1k, synthetic_data_1k):
    """
    Validate hypothesis: delta-learning outperforms direct learning.

    This test implements the decision gate experiment to validate
    that learning corrections (delta) is more effective than learning
    absolute values (direct) for CBS energy prediction.
    """

    # ============================================================
    # STEP 1: Data Loading with Fallback
    # ============================================================
    try:
        X, y_cbs, y_pm7, y_delta = real_data_1k
        data_source = "real_data_1k"
    except Exception as e:
        print(f"⚠️  Real data load failed: {e}")
        print("📊 Falling back to synthetic_data_1k")
        X, y_cbs, y_pm7, y_delta = synthetic_data_1k
        data_source = "synthetic_data_1k"

    # Log experiment header and data info
    print(f"\n{'='*60}")
    print(f"🔬 HYPOTHESIS VALIDATION EXPERIMENT")
    print(f"{'='*60}")
    print(f"Data source: {data_source}")
    print(f"Samples: {X.shape[0]}")
    print(f"Features: {X.shape[1]}")
    print(f"y_cbs range: [{y_cbs.min():.2f}, {y_cbs.max():.2f}] kcal/mol")
    print(f"y_pm7 range: [{y_pm7.min():.2f}, {y_pm7.max():.2f}] kcal/mol")
    print(f"delta range: [{y_delta.min():.2f}, {y_delta.max():.2f}] kcal/mol")

    # ============================================================
    # STEP 2: Train/Test Split (Same Split for Both Models)
    # ============================================================
    X_train, X_test, y_cbs_train, y_cbs_test, y_pm7_train, y_pm7_test = train_test_split(
        X, y_cbs, y_pm7,
        test_size=0.2,
        random_state=42
    )

    print(f"\n{'─'*60}")
    print(f"📊 DATA SPLIT")
    print(f"{'─'*60}")
    print(f"Training samples: {len(X_train)}")
    print(f"Testing samples: {len(X_test)}")

    # ============================================================
    # STEP 3: Model Delta (Reference Approach)
    # ============================================================
    print(f"\n{'─'*60}")
    print(f"🔷 MODEL 1: DELTA LEARNING (Reference)")
    print(f"{'─'*60}")
    print("Training DeltaLearner on y_delta = y_cbs - y_pm7...")

    # Initialize and train
    model_delta = DeltaLearner()  # Uses default w_krr=0.5, w_xgb=0.5
    model_delta.fit(X_train, y_cbs_train, y_pm7_train)

    # Evaluate on test set
    metrics_delta = model_delta.evaluate(X_test, y_cbs_test, y_pm7_test)

    # Extract key metrics
    rmse_delta = metrics_delta['rmse']
    mae_delta = metrics_delta['mae']
    r2_delta = metrics_delta['r2']

    print(f"✅ Training complete")
    print(f"RMSE: {rmse_delta:.4f} kcal/mol")
    print(f"MAE:  {mae_delta:.4f} kcal/mol")
    print(f"R²:   {r2_delta:.4f}")

    # ============================================================
    # STEP 4: Model Direct (Baseline)
    # ============================================================
    print(f"\n{'─'*60}")
    print(f"🔶 MODEL 2: DIRECT LEARNING (Baseline)")
    print(f"{'─'*60}")
    print("Training ensemble DIRECTLY on y_cbs (no delta)...")

    # Initialize and train directly on CBS
    model_direct = DeltaLearningEnsemble(w_krr=0.5, w_xgb=0.5)
    model_direct.fit(X_train, y_cbs_train)  # Train on CBS directly

    # Predict and evaluate
    y_pred_direct = model_direct.predict(X_test)
    metrics_direct = compute_all_metrics(y_cbs_test, y_pred_direct)

    # Extract key metrics
    rmse_direct = metrics_direct['rmse']
    mae_direct = metrics_direct['mae']
    r2_direct = metrics_direct['r2']

    print(f"✅ Training complete")
    print(f"RMSE: {rmse_direct:.4f} kcal/mol")
    print(f"MAE:  {mae_direct:.4f} kcal/mol")
    print(f"R²:   {r2_direct:.4f}")

    # ============================================================
    # STEP 5: Comparison & Decision Gate
    # ============================================================
    print(f"\n{'='*60}")
    print(f"⚖️  COMPARISON: DELTA vs DIRECT")
    print(f"{'='*60}")

    # Calculate improvements
    rmse_improvement = rmse_direct - rmse_delta
    mae_improvement = mae_direct - mae_delta
    r2_improvement = r2_delta - r2_direct

    print(f"{'Metric':<15} {'Delta':<15} {'Direct':<15} {'Improvement':<15}")
    print(f"{'-'*60}")
    print(f"{'RMSE':<15} {rmse_delta:<15.4f} {rmse_direct:<15.4f} {rmse_improvement:<15.4f}")
    print(f"{'MAE':<15} {mae_delta:<15.4f} {mae_direct:<15.4f} {mae_improvement:<15.4f}")
    print(f"{'R²':<15} {r2_delta:<15.4f} {r2_direct:<15.4f} {r2_improvement:<15.4f}")

    # Decision gate criteria
    criterion_1 = rmse_delta < rmse_direct
    criterion_2 = rmse_delta < 20.0

    print(f"\n{'─'*60}")
    print(f"🚦 DECISION GATE")
    print(f"{'─'*60}")
    print(f"Criterion 1: RMSE_delta < RMSE_direct")
    print(f"  → {rmse_delta:.4f} < {rmse_direct:.4f} = {criterion_1}")
    print(f"Criterion 2: RMSE_delta < 20.0 kcal/mol")
    print(f"  → {rmse_delta:.4f} < 20.0 = {criterion_2}")

    gate_pass = criterion_1 and criterion_2

    if gate_pass:
        print(f"\n✅ DECISION GATE: PASS")
        print(f"   Hypothesis validated: Delta-learning outperforms direct learning")
        print(f"   RMSE improvement: {rmse_improvement:.4f} kcal/mol")
    else:
        print(f"\n❌ DECISION GATE: FAIL")
        if not criterion_1:
            print(f"   ⚠️  Delta learning did NOT outperform direct learning")
            print(f"   RMSE_delta ({rmse_delta:.4f}) >= RMSE_direct ({rmse_direct:.4f})")
        if not criterion_2:
            print(f"   ⚠️  RMSE quality threshold not met")
            print(f"   RMSE_delta ({rmse_delta:.4f}) >= 20.0 kcal/mol")

    print(f"{'='*60}\n")

    # ============================================================
    # STEP 6: Assert gate pass (test fails if hypothesis not validated)
    # ============================================================
    assert gate_pass, (
        f"Decision gate failed: "
        f"criterion_1={criterion_1}, "
        f"criterion_2={criterion_2}"
    )
