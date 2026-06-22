"""Visualization charts for Phase G.

Generates publication-quality dark-themed PNG charts from ML predictions:
- Parity plot (predicted vs CBS reference)
- Delta histogram (prediction error distribution)
- Residuals plot (residuals vs predicted values)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import numpy.typing as npt  # noqa: E402
import pandas as pd  # noqa: E402

logger = logging.getLogger(__name__)

_BG_COLOR = "#1a1a1a"
FloatArray = npt.NDArray[np.float64]


@dataclass
class ChartGenerationResult:
    """Result of chart generation with file paths and metrics."""

    parity_plot: Path
    delta_histogram: Path
    residuals_plot: Path
    n_points: int
    rmse: float
    r2: float
    abs_error_histogram: Path | None = None
    abs_error_boxplot: Path | None = None
    abs_error_vs_molwt: Path | None = None


def generate_charts(
    csv_path: Path | str,
    output_dir: Path | str,
    extra_charts: bool = False,
) -> ChartGenerationResult:
    """Generate all 3 analysis charts from a CSV with predictions.

    Args:
        csv_path: Path to CSV containing ``H298_cbs`` and ``H298_predicted``.
        output_dir: Directory where PNG files will be saved.

    Returns:
        ChartGenerationResult with file paths and computed metrics.

    Raises:
        FileNotFoundError: If ``csv_path`` does not exist.
        ValueError: If CSV lacks ``H298_predicted`` column or has no valid rows.
    """
    csv_path = Path(csv_path)
    output_dir = Path(output_dir)

    if not csv_path.exists():
        msg = f"CSV file not found: {csv_path}"
        raise FileNotFoundError(msg)

    df = pd.read_csv(csv_path)

    if "H298_predicted" not in df.columns:
        msg = "CSV does not contain 'H298_predicted' column. Run predictions first."
        raise ValueError(msg)

    # Keep only rows with both CBS and predicted values
    valid = df.dropna(subset=["H298_cbs", "H298_predicted"])
    if len(valid) == 0:
        msg = "No rows with both H298_cbs and H298_predicted values."
        raise ValueError(msg)

    y_cbs = valid["H298_cbs"].to_numpy(dtype=np.float64)
    y_pred = valid["H298_predicted"].to_numpy(dtype=np.float64)

    # Compute metrics
    residuals = y_pred - y_cbs
    rmse = float(np.sqrt(np.mean(residuals**2)))
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y_cbs - np.mean(y_cbs)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate charts
    parity_path = _plot_parity(y_cbs, y_pred, rmse, r2, output_dir)
    histogram_path = _plot_delta_histogram(residuals, output_dir)
    residuals_path = _plot_residuals(y_pred, residuals, output_dir)

    n_points = len(valid)
    n_extra = 0
    abs_error_histogram: Path | None = None
    abs_error_boxplot: Path | None = None
    abs_error_vs_molwt: Path | None = None

    if extra_charts:
        abs_errors = np.abs(residuals)
        abs_error_histogram = _plot_abs_error_histogram(abs_errors, output_dir)
        abs_error_boxplot = _plot_abs_error_boxplot(abs_errors, output_dir)
        if "rdkit_mol_weight" in valid.columns:
            mol_weights = valid["rdkit_mol_weight"].to_numpy(dtype=np.float64)
            abs_error_vs_molwt = _plot_abs_error_vs_molwt(
                abs_errors, mol_weights, output_dir
            )
        n_extra = 2 + (1 if abs_error_vs_molwt is not None else 0)

    total = 3 + n_extra
    logger.info("Generated %d charts from %d data points", total, n_points)

    return ChartGenerationResult(
        parity_plot=parity_path,
        delta_histogram=histogram_path,
        residuals_plot=residuals_path,
        n_points=n_points,
        rmse=rmse,
        r2=r2,
        abs_error_histogram=abs_error_histogram,
        abs_error_boxplot=abs_error_boxplot,
        abs_error_vs_molwt=abs_error_vs_molwt,
    )


def _plot_parity(
    y_cbs: FloatArray,
    y_pred: FloatArray,
    rmse: float,
    r2: float,
    output_dir: Path,
) -> Path:
    """Create parity plot (predicted vs CBS reference).

    Args:
        y_cbs: CBS reference values.
        y_pred: Predicted values.
        rmse: Root mean squared error.
        r2: Coefficient of determination.
        output_dir: Output directory for the PNG.

    Returns:
        Path to the saved PNG file.
    """
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(8, 7), dpi=150)
    fig.set_facecolor(_BG_COLOR)
    ax.set_facecolor(_BG_COLOR)

    ax.scatter(y_cbs, y_pred, c="#00D9FF", alpha=0.7, s=40, edgecolors="none")

    # y=x reference line
    all_vals = np.concatenate([y_cbs, y_pred])
    lo = float(np.min(all_vals)) - 5
    hi = float(np.max(all_vals)) + 5
    ax.plot([lo, hi], [lo, hi], "--", color="#666666", linewidth=1, label="y = x")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)

    # Annotation
    ax.text(
        0.05,
        0.95,
        f"R\u00b2 = {r2:.4f}\nRMSE = {rmse:.3f} kcal/mol",
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox={"boxstyle": "round", "facecolor": "#333333", "alpha": 0.8},
        color="white",
    )

    ax.set_xlabel("CBS Reference (kcal/mol)", fontsize=12)
    ax.set_ylabel("Predicted (kcal/mol)", fontsize=12)
    ax.set_title("Parity Plot: Predicted vs CBS", fontsize=14, fontweight="bold")
    ax.legend(loc="lower right")

    path = output_dir / "parity_plot.png"
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def _plot_delta_histogram(
    delta: FloatArray,
    output_dir: Path,
) -> Path:
    """Create histogram of prediction errors (delta = predicted - CBS).

    Args:
        delta: Array of prediction errors.
        output_dir: Output directory for the PNG.

    Returns:
        Path to the saved PNG file.
    """
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fig.set_facecolor(_BG_COLOR)
    ax.set_facecolor(_BG_COLOR)

    ax.hist(
        delta,
        bins=min(50, max(5, len(delta))),
        color="#FF00FF",
        alpha=0.7,
        edgecolor="#CC00CC",
    )

    # Reference lines
    ax.axvline(x=0, color="#666666", linestyle="--", linewidth=1, label="Zero")
    mean_delta = float(np.mean(delta))
    ax.axvline(
        x=mean_delta,
        color="#FFFF00",
        linestyle="-.",
        linewidth=1,
        label=f"Mean = {mean_delta:.3f}",
    )

    ax.set_xlabel("\u0394 (Predicted \u2212 CBS) [kcal/mol]", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Delta Distribution", fontsize=14, fontweight="bold")
    ax.legend()

    path = output_dir / "delta_histogram.png"
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def _plot_residuals(
    y_pred: FloatArray,
    residuals: FloatArray,
    output_dir: Path,
) -> Path:
    """Create residuals vs predicted values scatter plot.

    Args:
        y_pred: Predicted values.
        residuals: Residual errors (predicted - CBS).
        output_dir: Output directory for the PNG.

    Returns:
        Path to the saved PNG file.
    """
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fig.set_facecolor(_BG_COLOR)
    ax.set_facecolor(_BG_COLOR)

    ax.scatter(y_pred, residuals, c="#00FF80", alpha=0.5, s=40, edgecolors="none")

    # Reference lines
    xlims = ax.get_xlim()
    ax.axhline(y=0, color="#666666", linestyle="--", linewidth=1)
    ax.axhline(y=5, color="#FF6666", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.axhline(y=-5, color="#FF6666", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.set_xlim(xlims)

    ax.set_xlabel("Predicted H298 (kcal/mol)", fontsize=12)
    ax.set_ylabel("Residual (kcal/mol)", fontsize=12)
    ax.set_title("Residuals vs Predicted", fontsize=14, fontweight="bold")

    path = output_dir / "residuals_plot.png"
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def _plot_abs_error_histogram(
    abs_errors: FloatArray,
    output_dir: Path,
) -> Path:
    """Histogram of absolute errors.

    Args:
        abs_errors: Array of absolute prediction errors.
        output_dir: Output directory for the PNG.

    Returns:
        Path to the saved PNG file.
    """
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fig.set_facecolor(_BG_COLOR)
    ax.set_facecolor(_BG_COLOR)

    ax.hist(
        abs_errors,
        bins=min(50, max(5, len(abs_errors))),
        color="#00D9FF",
        alpha=0.7,
        edgecolor="#009AB8",
    )
    mean_err = float(np.mean(abs_errors))
    ax.axvline(
        x=mean_err,
        color="#FFFF00",
        linestyle="-.",
        linewidth=1,
        label=f"Mean = {mean_err:.3f}",
    )
    ax.legend()

    ax.set_xlabel("|Error| (kcal/mol)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Absolute Error Distribution", fontsize=14, fontweight="bold")

    path = output_dir / "abs_error_histogram.png"
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def _plot_abs_error_boxplot(
    abs_errors: FloatArray,
    output_dir: Path,
) -> Path:
    """Boxplot of absolute errors.

    Args:
        abs_errors: Array of absolute prediction errors.
        output_dir: Output directory for the PNG.

    Returns:
        Path to the saved PNG file.
    """
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(6, 7), dpi=150)
    fig.set_facecolor(_BG_COLOR)
    ax.set_facecolor(_BG_COLOR)

    bp = ax.boxplot(
        abs_errors,
        patch_artist=True,
        notch=False,
        vert=True,
        widths=0.5,
    )
    bp["boxes"][0].set_facecolor("#00D9FF")
    bp["boxes"][0].set_alpha(0.7)
    for element in ("whiskers", "caps", "medians"):
        for line in bp[element]:
            line.set_color("#FFFFFF")
    for flier in bp["fliers"]:
        flier.set(marker="o", color="#FF6666", alpha=0.5, markersize=4)

    ax.set_ylabel("|Error| (kcal/mol)", fontsize=12)
    ax.set_xticks([])
    ax.set_title("Absolute Error Boxplot", fontsize=14, fontweight="bold")

    path = output_dir / "abs_error_boxplot.png"
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def _plot_abs_error_vs_molwt(
    abs_errors: FloatArray,
    mol_weights: FloatArray,
    output_dir: Path,
) -> Path:
    """Scatter of absolute error vs molecular weight.

    Args:
        abs_errors: Array of absolute prediction errors.
        mol_weights: Array of molecular weights (RDKit).
        output_dir: Output directory for the PNG.

    Returns:
        Path to the saved PNG file.
    """
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(9, 6), dpi=150)
    fig.set_facecolor(_BG_COLOR)
    ax.set_facecolor(_BG_COLOR)

    ax.scatter(mol_weights, abs_errors, c="#FF00FF", alpha=0.5, s=30, edgecolors="none")

    ax.set_xlabel("Molecular Weight (g/mol)", fontsize=12)
    ax.set_ylabel("|Error| (kcal/mol)", fontsize=12)
    ax.set_title("|Error| vs Molecular Weight", fontsize=14, fontweight="bold")

    path = output_dir / "abs_error_vs_molwt.png"
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path
