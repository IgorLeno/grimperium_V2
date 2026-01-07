"""
BATCH 3 - DECISION GATE TEST (Hypothesis Validation)

HYPOTHESIS:
-----------
"Delta-learning is superior to direct learning because it learns a correction
of lower variance in a realistic regime."

METHOD:
-------
- Data: 1000 molecules with realistic distribution (-1000 to +1000 kcal/mol)
- Split: 80/20, same train/test (no distribution shift)
- Models: DeltaLearner vs DeltaLearningEnsemble directly on y_cbs
- Metrics: RMSE, MAE, R2 (interpretable scale)

SUCCESS CRITERIA (Decision Gate):
---------------------------------
1. rmse_delta < rmse_direct (delta wins in absolute terms)
2. rmse_delta < 20 kcal/mol (acceptable quality)

EXPECTED RESULTS:
-----------------
- rmse_delta: 10-15 kcal/mol
- rmse_direct: 20-50 kcal/mol
- Improvement: 2-5x (NOT 73x - that's an outlier artifact)

OUTCOME:
--------
If both criteria pass: Hypothesis validated in realistic regime
If any fails: Requires adjustments in features or architecture

WHY THIS DESIGN:
----------------
The original implementation used unfiltered data with H298_cbs outliers
up to -325407 kcal/mol. Random sampling created severe distribution shift
(train mean=-594, test mean=+21), causing Direct model to fail with
RMSE=1008 and R2=-1.09. While this demonstrates delta-learning robustness,
it's not a fair hypothesis test.

This version uses filtered data (-1000 to +1000 kcal/mol) to:
1. Remove the 0.9% extreme outliers
2. Ensure train/test have similar distributions
3. Provide meaningful comparison (2-5x improvement, not 73x)
4. Generate interpretable metrics for Week 2/3 decisions

For robustness testing with extreme data, see test_stress_distribution_shift.py
"""

import numpy as np
from sklearn.model_selection import train_test_split

from grimperium.core.delta_learning import DeltaLearner
from grimperium.models.delta_ensemble import DeltaLearningEnsemble
from grimperium.core.metrics import compute_all_metrics


def test_decision_gate_delta_vs_direct(real_data_1k_filtered):
    """
    Validate core hypothesis: delta-learning outperforms direct learning
    (realistic regime).

    Main BATCH 3 decision gate test. Uses filtered data to ensure a fair
    comparison without distribution shift artifacts.

    Expected outcome:
    - rmse_delta: 10-15 kcal/mol
    - rmse_direct: 20-50 kcal/mol
    - Improvement: 2-5x

    If this test passes, the hypothesis is validated for realistic use cases.
    """

    print("\n" + "=" * 70)
    print("BATCH 3 - DECISION GATE TEST (Hypothesis Validation)")
    print("=" * 70)

    # ================================================================
    # STEP 1: Load data with realistic distribution
    # ================================================================
    print("\n[STEP 1] Loading data (realistic, no extreme outliers)...")
    X, y_cbs, y_pm7 = real_data_1k_filtered

    print(f"  Data shape: X={X.shape}, y_cbs={y_cbs.shape}, y_pm7={y_pm7.shape}")
    print(f"  y_cbs statistics: mean={np.mean(y_cbs):.1f}, " f"std={np.std(y_cbs):.1f}")
    print(f"  y_cbs range: [{np.min(y_cbs):.1f}, {np.max(y_cbs):.1f}]")

    # ================================================================
    # STEP 2: Train/test split (80/20, no distribution shift)
    # ================================================================
    print("\n[STEP 2] Train/test split (80/20, random_state=42)...")
    (
        X_train,
        X_test,
        y_cbs_train,
        y_cbs_test,
        y_pm7_train,
        y_pm7_test,
    ) = train_test_split(X, y_cbs, y_pm7, test_size=0.2, random_state=42)

    print(f"  Train: {len(X_train)} samples")
    print(f"  Test:  {len(X_test)} samples")

    # Validate: no distribution shift
    train_mean = np.mean(y_cbs_train)
    test_mean = np.mean(y_cbs_test)
    mean_diff = abs(train_mean - test_mean)

    print(f"  Train y_cbs mean: {train_mean:.1f}")
    print(f"  Test y_cbs mean:  {test_mean:.1f}")
    print(
        f"  Mean difference: {mean_diff:.1f} " "(should be < 100 for fair comparison)"
    )

    # With filtered data, distribution shift should be minimal
    assert mean_diff < 100, (
        f"Distribution shift detected: {mean_diff:.1f} > 100. " "Check data filtering."
    )

    # ================================================================
    # STEP 3: Model Delta (reference approach)
    # ================================================================
    print("\n[STEP 3] Training DELTA model (learns y_delta = y_cbs - y_pm7)...")
    model_delta = DeltaLearner()
    model_delta.fit(X_train, y_cbs_train, y_pm7_train)

    metrics_delta = model_delta.evaluate(X_test, y_cbs_test, y_pm7_test)
    rmse_delta = metrics_delta["rmse"]
    mae_delta = metrics_delta["mae"]
    r2_delta = metrics_delta["r2"]

    print(f"  RMSE: {rmse_delta:.2f} kcal/mol (expected: 10-15)")
    print(f"  MAE:  {mae_delta:.2f} kcal/mol")
    print(f"  R2:   {r2_delta:.4f}")

    # ================================================================
    # STEP 4: Model Direct (baseline)
    # ================================================================
    print("\n[STEP 4] Training DIRECT model (learns y_cbs directly)...")
    model_direct = DeltaLearningEnsemble(w_krr=0.5, w_xgb=0.5)
    model_direct.fit(X_train, y_cbs_train)

    y_pred_direct = model_direct.predict(X_test)
    metrics_direct = compute_all_metrics(y_cbs_test, y_pred_direct)
    rmse_direct = metrics_direct["rmse"]
    mae_direct = metrics_direct["mae"]
    r2_direct = metrics_direct["r2"]

    print(f"  RMSE: {rmse_direct:.2f} kcal/mol (expected: 20-50)")
    print(f"  MAE:  {mae_direct:.2f} kcal/mol")
    print(f"  R2:   {r2_direct:.4f}")

    # ================================================================
    # STEP 5: Comparison
    # ================================================================
    print("\n[STEP 5] Comparison: Delta vs Direct...")
    print("-" * 70)
    print(f"{'Metric':<15} {'Delta':<15} {'Direct':<15} {'Winner':<10}")
    print("-" * 70)

    rmse_improvement = rmse_direct - rmse_delta
    rmse_improvement_pct = (
        (rmse_improvement / rmse_direct) * 100 if rmse_direct > 0 else 0
    )
    rmse_ratio = rmse_direct / rmse_delta if rmse_delta > 0 else float("inf")

    rmse_winner = "DELTA" if rmse_delta < rmse_direct else "DIRECT"
    mae_winner = "DELTA" if mae_delta < mae_direct else "DIRECT"
    r2_winner = "DELTA" if r2_delta > r2_direct else "DIRECT"

    print(
        f"{'RMSE':<15} {rmse_delta:<15.2f} {rmse_direct:<15.2f} " f"{rmse_winner:<10}"
    )
    print(f"{'MAE':<15} {mae_delta:<15.2f} {mae_direct:<15.2f} {mae_winner:<10}")
    print(f"{'R2':<15} {r2_delta:<15.4f} {r2_direct:<15.4f} {r2_winner:<10}")
    print("-" * 70)
    print(
        f"RMSE improvement: {rmse_improvement:.2f} kcal/mol "
        f"({rmse_improvement_pct:.1f}%)"
    )
    print(f"Delta is {rmse_ratio:.1f}x better than Direct")

    # ================================================================
    # STEP 6: Decision Gate
    # ================================================================
    print("\n[STEP 6] DECISION GATE...")
    print("-" * 70)

    criterion_1 = rmse_delta < rmse_direct
    criterion_2 = rmse_delta < 20.0
    gate_pass = criterion_1 and criterion_2

    print("  Criterion 1: rmse_delta < rmse_direct")
    print(f"    {rmse_delta:.2f} < {rmse_direct:.2f}? -> {criterion_1}")

    print("\n  Criterion 2: rmse_delta < 20.0 kcal/mol")
    print(f"    {rmse_delta:.2f} < 20.0? -> {criterion_2}")

    print("-" * 70)

    if gate_pass:
        print("  DECISION GATE: PASS")
        print("  Hypothesis validated in realistic regime!")
        print(
            f"  Delta-learning is {rmse_improvement_pct:.1f}% better than "
            "direct learning"
        )
    else:
        print("  DECISION GATE: FAIL")
        if not criterion_1:
            print(
                "  Criterion 1 failed: "
                f"Delta ({rmse_delta:.2f}) >= Direct ({rmse_direct:.2f})"
            )
        if not criterion_2:
            print(f"  Criterion 2 failed: Delta RMSE ({rmse_delta:.2f}) >= 20.0")

    print("=" * 70 + "\n")

    # ================================================================
    # STEP 7: Assert
    # ================================================================
    assert gate_pass, (
        f"Decision gate failed: c1={criterion_1} "
        f"(rmse_delta={rmse_delta:.2f}, rmse_direct={rmse_direct:.2f}), "
        f"c2={criterion_2}"
    )

    return {
        "rmse_delta": rmse_delta,
        "rmse_direct": rmse_direct,
        "rmse_improvement": rmse_improvement,
        "rmse_improvement_pct": rmse_improvement_pct,
        "rmse_ratio": rmse_ratio,
        "metrics_delta": metrics_delta,
        "metrics_direct": metrics_direct,
        "gate_pass": gate_pass,
        "regime": "realistic (filtered, no distribution shift)",
    }


def test_synthetic_fallback(synthetic_data_1k):
    """
    Fallback test using synthetic data.

    Use when real data is unavailable or for fast CI runs.
    Synthetic data has controlled distribution, so results may differ
    from real data.
    """
    print("\n" + "=" * 70)
    print("FALLBACK TEST (Synthetic Data)")
    print("=" * 70)

    X, y_cbs, y_pm7 = synthetic_data_1k

    # Quick validation with synthetic data
    (
        X_train,
        X_test,
        y_cbs_train,
        y_cbs_test,
        y_pm7_train,
        y_pm7_test,
    ) = train_test_split(X, y_cbs, y_pm7, test_size=0.2, random_state=42)

    # Delta model
    model_delta = DeltaLearner()
    model_delta.fit(X_train, y_cbs_train, y_pm7_train)
    metrics_delta = model_delta.evaluate(X_test, y_cbs_test, y_pm7_test)

    # Direct model
    model_direct = DeltaLearningEnsemble()
    model_direct.fit(X_train, y_cbs_train)
    y_pred_direct = model_direct.predict(X_test)
    metrics_direct = compute_all_metrics(y_cbs_test, y_pred_direct)

    print(f"\n  Delta RMSE: {metrics_delta['rmse']:.2f}")
    print(f"  Direct RMSE: {metrics_direct['rmse']:.2f}")
    print(f"  Delta wins: {metrics_delta['rmse'] < metrics_direct['rmse']}")
    print("=" * 70 + "\n")

    # Synthetic data should also show delta advantage
    assert (
        metrics_delta["rmse"] < metrics_direct["rmse"]
    ), "Delta should outperform direct even on synthetic data"
