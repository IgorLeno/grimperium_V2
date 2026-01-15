"""
Reproduz os benchmarks do relatório técnico: BATCH 3 (Hypothesis Validation).

Este script roda exatamente o mesmo protocolo codificado em:
- tests/experiments/test_validate_hypothesis.py (regime realista/filtrado)
- tests/experiments/test_stress_distribution_shift.py
  (regime extremo/não filtrado)
- tests/experiments/conftest.py (filtro, amostragem, PM7 mock e seeds)

Além do "point estimate" (seed=42), o script calcula incerteza via bootstrap
no conjunto de teste (IC 95% por percentis) para:
- RMSE (Delta Learning)
- RMSE (Direct Model)
- R² (Delta Learning)
- R² (Direct Model)
- Razão (Direct/Delta) no regime extremo (robustness ratio)

Uso:
  python scripts/benchmarks/reproduce_batch3.py --data-path thermo_cbs_opt.csv
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from grimperium.core.delta_learning import DeltaLearner
from grimperium.core.metrics import compute_all_metrics
from grimperium.data.loader import ChemperiumLoader
from grimperium.models.delta_ensemble import DeltaLearningEnsemble


class Metrics(TypedDict):
    """Métricas retornadas por `compute_all_metrics`."""

    rmse: float
    mae: float
    r2: float
    mape: float
    max_error: float


class FitPredictOutput(TypedDict):
    """
    Saída tipada de `_fit_predict`.

    Mantém shape estável e elimina `type: ignore` em acessos por chave.
    """

    y_true: np.ndarray
    y_pred_delta: np.ndarray
    y_pred_direct: np.ndarray
    metrics_delta: Metrics
    metrics_direct: Metrics
    rmse_ratio: float


class SummaryStats(TypedDict):
    mean: float
    std: float
    ci95_lo: float
    ci95_hi: float


class BootstrapCI95(TypedDict):
    rmse_delta: SummaryStats
    rmse_direct: SummaryStats
    r2_delta: SummaryStats
    r2_direct: SummaryStats
    rmse_ratio: SummaryStats


class Env(TypedDict):
    git_commit: str
    python: str
    platform: str
    cpu: str
    mem_gb: float | None
    numpy: str
    pandas: str
    sklearn: str
    xgboost: str


class RegimeConfigJSON(TypedDict):
    filtered: bool
    filter_range_kcal_mol: list[int] | None
    n_samples: int
    test_size: float
    sample_seed: int
    split_seed: int
    pm7_mock_seed: int
    pm7_mock_magnitude_bias_std: float
    bootstrap_n: int
    bootstrap_seed: int


class PointEstimate(TypedDict):
    rmse_delta: float
    rmse_direct: float
    r2_delta: float
    r2_direct: float
    rmse_ratio: float


class RegimeResults(TypedDict):
    config: RegimeConfigJSON
    point_estimate: PointEstimate
    bootstrap_ci95: BootstrapCI95


class Results(TypedDict):
    env: Env
    protocol_seed: int
    regimes: dict[str, RegimeResults]


def _get_git_commit() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
        return out
    except Exception:
        return "unknown"


def _get_cpu_model() -> str:
    try:
        txt = Path("/proc/cpuinfo").read_text(
            encoding="utf-8", errors="ignore"
        )
        m = re.search(r"^model name\s*:\s*(.+)$", txt, re.M)
        return m.group(1).strip() if m else "unknown"
    except Exception:
        return "unknown"


def _get_mem_total_gb() -> float | None:
    try:
        meminfo = Path("/proc/meminfo").read_text(encoding="utf-8")
        for line in meminfo.splitlines():
            if line.startswith("MemTotal:"):
                kb = int(line.split()[1])
                return kb / 1024 / 1024
    except Exception:
        return None
    return None


def _create_realistic_mock_pm7(
    y_cbs: np.ndarray,
    X_basic: np.ndarray,
    *,
    seed: int | None = 42,
    magnitude_bias_std: float = 0.5,
) -> np.ndarray:
    """
    Copiado 1:1 de `tests/experiments/conftest.py::create_realistic_mock_pm7`.

    Nota: usa `np.random.seed(seed)` (global) por compatibilidade com o
    protocolo atual dos experimentos.
    """
    if seed is not None:
        np.random.seed(seed)
    nheavy = X_basic[:, 0]

    # Component 1: Base bias (PM7 systematically overestimates)
    base_bias = -5.0

    # Component 2: Size-dependent error
    size_error = (1 + np.sqrt(np.maximum(nheavy, 1))) * np.random.normal(
        0, 1.5, len(y_cbs)
    )

    # Component 3: Magnitude-dependent bias (2% de |y_cbs| por padrão).
    magnitude_bias = (np.abs(y_cbs) / 50) * np.random.normal(
        0, magnitude_bias_std, len(y_cbs)
    )

    # Component 4: Random Gaussian noise
    gaussian_noise = np.random.normal(0, 7, len(y_cbs))

    total_error = base_bias + size_error + magnitude_bias + gaussian_noise
    return y_cbs - total_error


def _create_enriched_features(X_basic: np.ndarray) -> np.ndarray:
    """
    Copiado 1:1 de `tests/experiments/conftest.py::create_enriched_features`.
    """
    nheavy = X_basic[:, 0:1]
    charge = X_basic[:, 1:2]
    mult = X_basic[:, 2:3]

    return np.hstack(
        [
            nheavy,
            charge,
            mult,
            nheavy**2,
            np.sqrt(np.abs(nheavy) + 1),
            nheavy * charge,
            nheavy * mult,
            charge**2,
            mult**2,
            np.ones_like(nheavy),
        ]
    )


@dataclass(frozen=True)
class RegimeConfig:
    name: str
    filtered: bool
    magnitude_bias_std: float


def _load_xy(
    *,
    data_path: Path,
    n_samples: int,
    sample_seed: int,
    filtered: bool,
    magnitude_bias_std: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Match the fixture's behavior: use ChemperiumLoader + full dataset file.
    df = ChemperiumLoader.load_thermo_cbs_opt(path=data_path, validate=True)

    if filtered:
        df = df[(df["H298_cbs"] >= -1000) & (df["H298_cbs"] <= 1000)]

    df_sample = df.sample(n=n_samples, random_state=sample_seed)

    X_basic = (
        df_sample[["nheavy", "charge", "multiplicity"]].values.astype(float)
    )
    y_cbs = df_sample["H298_cbs"].values.astype(float)
    X = _create_enriched_features(X_basic)
    y_pm7 = _create_realistic_mock_pm7(
        y_cbs, X_basic, seed=sample_seed, magnitude_bias_std=magnitude_bias_std
    )
    return X, y_cbs, y_pm7


def _fit_predict(
    *,
    X: np.ndarray,
    y_cbs: np.ndarray,
    y_pm7: np.ndarray,
    split_seed: int,
    test_size: float,
) -> FitPredictOutput:
    (
        X_train,
        X_test,
        y_cbs_train,
        y_cbs_test,
        y_pm7_train,
        y_pm7_test,
    ) = train_test_split(
        X, y_cbs, y_pm7, test_size=test_size, random_state=split_seed
    )

    model_delta = DeltaLearner()
    model_delta.fit(X_train, y_cbs_train, y_pm7_train)
    y_pred_delta = model_delta.predict(X_test, y_pm7_test)
    metrics_delta = compute_all_metrics(y_cbs_test, y_pred_delta)

    model_direct = DeltaLearningEnsemble(w_krr=0.5, w_xgb=0.5)
    model_direct.fit(X_train, y_cbs_train)
    y_pred_direct = model_direct.predict(X_test)
    metrics_direct = compute_all_metrics(y_cbs_test, y_pred_direct)

    # Guard contra divisão por zero (ou valores extremamente pequenos).
    denom = float(metrics_delta["rmse"])
    if abs(denom) < 1e-12:
        rmse_ratio = float("inf")
    else:
        rmse_ratio = float(metrics_direct["rmse"]) / denom

    return {
        "y_true": y_cbs_test,
        "y_pred_delta": y_pred_delta,
        "y_pred_direct": y_pred_direct,
        "metrics_delta": metrics_delta,
        "metrics_direct": metrics_direct,
        "rmse_ratio": rmse_ratio,
    }


def _bootstrap_ci(
    *,
    y_true: np.ndarray,
    y_pred_delta: np.ndarray,
    y_pred_direct: np.ndarray,
    n_bootstrap: int,
    seed: int,
) -> BootstrapCI95:
    rng = np.random.default_rng(seed)
    n = len(y_true)
    idx = np.arange(n)

    rmse_deltas: list[float] = []
    rmse_directs: list[float] = []
    r2_deltas: list[float] = []
    r2_directs: list[float] = []
    ratios: list[float] = []

    for _ in range(n_bootstrap):
        b = rng.choice(idx, size=n, replace=True)
        m_delta = compute_all_metrics(y_true[b], y_pred_delta[b])
        m_direct = compute_all_metrics(y_true[b], y_pred_direct[b])
        rmse_deltas.append(float(m_delta["rmse"]))
        rmse_directs.append(float(m_direct["rmse"]))
        r2_deltas.append(float(m_delta["r2"]))
        r2_directs.append(float(m_direct["r2"]))
        denom = float(m_delta["rmse"])
        if abs(denom) < 1e-12:
            ratios.append(float("inf"))
        else:
            ratios.append(float(m_direct["rmse"]) / denom)

    def summarize(xs: list[float]) -> SummaryStats:
        arr = np.asarray(xs, dtype=float)
        lo, hi = np.percentile(arr, [2.5, 97.5])
        return {
            "mean": float(arr.mean()),
            "std": float(arr.std(ddof=1)),
            "ci95_lo": float(lo),
            "ci95_hi": float(hi),
        }

    return {
        "rmse_delta": summarize(rmse_deltas),
        "rmse_direct": summarize(rmse_directs),
        "r2_delta": summarize(r2_deltas),
        "r2_direct": summarize(r2_directs),
        "rmse_ratio": summarize(ratios),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-path",
        type=Path,
        default=Path("thermo_cbs_opt.csv"),
        help=(
            "Caminho para o CSV do Chemperium (padrão: ./thermo_cbs_opt.csv)."
        ),
    )
    parser.add_argument("--n-samples", type=int, default=1000)
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument(
        "--seed", type=int, default=42, help="Seed do protocolo."
    )
    parser.add_argument("--bootstrap-n", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=42)
    args = parser.parse_args()

    env: Env = {
        "git_commit": _get_git_commit(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "cpu": _get_cpu_model(),
        "mem_gb": _get_mem_total_gb(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "sklearn": "unknown",
        "xgboost": "unknown",
    }
    try:
        import sklearn  # type: ignore

        env["sklearn"] = sklearn.__version__
    except Exception:
        pass
    try:
        import xgboost  # type: ignore

        env["xgboost"] = xgboost.__version__
    except Exception:
        pass

    regimes = [
        RegimeConfig(
            name="realistic_filtered",
            filtered=True,
            magnitude_bias_std=0.5,
        ),
        RegimeConfig(
            name="extreme_unfiltered",
            filtered=False,
            magnitude_bias_std=3.0,
        ),
    ]

    results: Results = {"env": env, "protocol_seed": args.seed, "regimes": {}}

    for reg in regimes:
        X, y_cbs, y_pm7 = _load_xy(
            data_path=args.data_path,
            n_samples=args.n_samples,
            sample_seed=args.seed,
            filtered=reg.filtered,
            magnitude_bias_std=reg.magnitude_bias_std,
        )
        fit = _fit_predict(
            X=X,
            y_cbs=y_cbs,
            y_pm7=y_pm7,
            split_seed=args.seed,
            test_size=args.test_size,
        )
        ci = _bootstrap_ci(
            y_true=fit["y_true"],
            y_pred_delta=fit["y_pred_delta"],
            y_pred_direct=fit["y_pred_direct"],
            n_bootstrap=args.bootstrap_n,
            seed=args.bootstrap_seed,
        )

        filter_range = [-1000, 1000] if reg.filtered else None
        results["regimes"][reg.name] = {
            "config": {
                "filtered": reg.filtered,
                "filter_range_kcal_mol": filter_range,
                "n_samples": args.n_samples,
                "test_size": args.test_size,
                "sample_seed": args.seed,
                "split_seed": args.seed,
                "pm7_mock_seed": args.seed,
                "pm7_mock_magnitude_bias_std": reg.magnitude_bias_std,
                "bootstrap_n": args.bootstrap_n,
                "bootstrap_seed": args.bootstrap_seed,
            },
            "point_estimate": {
                "rmse_delta": float(fit["metrics_delta"]["rmse"]),
                "rmse_direct": float(fit["metrics_direct"]["rmse"]),
                "r2_delta": float(fit["metrics_delta"]["r2"]),
                "r2_direct": float(fit["metrics_direct"]["r2"]),
                "rmse_ratio": float(fit["rmse_ratio"]),
            },
            "bootstrap_ci95": ci,
        }

    print(json.dumps(results, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
