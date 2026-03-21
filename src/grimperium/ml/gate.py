# src/grimperium/ml/gate.py
"""Quality gate for trained DeltaLearner models.

Gate criteria reflect the minimum acceptable performance given 1 537 molecules
with delta mean ≈ 5 kcal/mol, std ≈ 6.45 kcal/mol:

    MAE  ≤ 3.5 kcal/mol  — naive mean-delta predictor scores ~5.15 kcal/mol;
                            3.5 marks a clear regression floor.
    R²   ≥ 0.97           — requires real variance explanation beyond the mean.
    RMSE ≤ 5.0 kcal/mol  — complements MAE; detects outlier amplification.

max_error and MAPE are intentionally excluded:
    - max_error is driven by organic radicals ([C], [CH]) — legitimate dataset members.
    - MAPE is numerically undefined when H298 is near zero (heats of formation cross zero).
"""

from __future__ import annotations

GATE_MAE_MAX: float = 3.5  # kcal/mol
GATE_R2_MIN: float = 0.97
GATE_RMSE_MAX: float = 5.0  # kcal/mol


def evaluate_gate(metrics: dict[str, float]) -> bool:
    """Return True if *metrics* satisfy all quality-gate criteria.

    Args:
        metrics: Dict with keys ``mae``, ``r2``, ``rmse`` (and optionally
                 others that are ignored). Missing keys default to worst-case
                 values (inf / 0.0) so a sparse dict fails the gate.

    Returns:
        True only when MAE ≤ GATE_MAE_MAX, R² ≥ GATE_R2_MIN,
        and RMSE ≤ GATE_RMSE_MAX simultaneously.
    """
    mae = float(metrics.get("mae", float("inf")))
    r2 = float(metrics.get("r2", 0.0))
    rmse = float(metrics.get("rmse", float("inf")))

    return mae <= GATE_MAE_MAX and r2 >= GATE_R2_MIN and rmse <= GATE_RMSE_MAX
