"""
Evaluation metrics for Grimperium.

This module provides metrics for evaluating delta-learning models:
    - RMSE: Root Mean Squared Error (primary metric)
    - MAE: Mean Absolute Error
    - R²: Coefficient of Determination
    - MAPE: Mean Absolute Percentage Error

All metrics are designed to work with both delta values
and full H298 predictions.

Example:
    >>> from grimperium.core.metrics import rmse, mae, r2_score
    >>> y_true = [1.0, 2.0, 3.0]
    >>> y_pred = [1.1, 2.2, 2.9]
    >>> print(f"RMSE: {rmse(y_true, y_pred):.4f}")
    >>> print(f"MAE: {mae(y_true, y_pred):.4f}")
    >>> print(f"R²: {r2_score(y_true, y_pred):.4f}")

"""

from typing import Optional, Union

import numpy as np


def rmse(
    y_true: Union[np.ndarray, list],
    y_pred: Union[np.ndarray, list],
    sample_weight: Optional[np.ndarray] = None,
) -> float:
    """
    Compute Root Mean Squared Error.

    RMSE = sqrt(mean((y_true - y_pred)²))

    Args:
        y_true: True values
        y_pred: Predicted values
        sample_weight: Optional sample weights

    Returns:
        RMSE value (lower is better)

    Example:
        >>> rmse([1, 2, 3], [1.1, 2.0, 2.9])
        0.0816...

    """
    y_true_arr: np.ndarray = np.asarray(y_true)
    y_pred_arr: np.ndarray = np.asarray(y_pred)

    squared_errors = (y_true_arr - y_pred_arr) ** 2

    if sample_weight is not None:
        sample_weight_arr: np.ndarray = np.asarray(sample_weight)
        return float(np.sqrt(np.average(squared_errors, weights=sample_weight_arr)))

    return float(np.sqrt(np.mean(squared_errors)))


def mae(
    y_true: Union[np.ndarray, list],
    y_pred: Union[np.ndarray, list],
    sample_weight: Optional[np.ndarray] = None,
) -> float:
    """
    Compute Mean Absolute Error.

    MAE = mean(|y_true - y_pred|)

    Args:
        y_true: True values
        y_pred: Predicted values
        sample_weight: Optional sample weights

    Returns:
        MAE value (lower is better)

    Example:
        >>> mae([1, 2, 3], [1.1, 2.0, 2.9])
        0.0666...

    """
    y_true_arr: np.ndarray = np.asarray(y_true)
    y_pred_arr: np.ndarray = np.asarray(y_pred)

    absolute_errors = np.abs(y_true_arr - y_pred_arr)

    if sample_weight is not None:
        sample_weight_arr: np.ndarray = np.asarray(sample_weight)
        return float(np.average(absolute_errors, weights=sample_weight_arr))

    return float(np.mean(absolute_errors))


def r2_score(
    y_true: Union[np.ndarray, list],
    y_pred: Union[np.ndarray, list],
    sample_weight: Optional[np.ndarray] = None,
) -> float:
    """
    Compute R² (Coefficient of Determination).

    R² = 1 - SS_res / SS_tot
    where:
        SS_res = sum((y_true - y_pred)²)
        SS_tot = sum((y_true - mean(y_true))²)

    Args:
        y_true: True values
        y_pred: Predicted values
        sample_weight: Optional sample weights

    Returns:
        R² value (higher is better, max 1.0)

    Example:
        >>> r2_score([1, 2, 3], [1.1, 2.0, 2.9])
        0.985...

    """
    y_true_arr: np.ndarray = np.asarray(y_true)
    y_pred_arr: np.ndarray = np.asarray(y_pred)

    if sample_weight is not None:
        sample_weight_arr: np.ndarray = np.asarray(sample_weight)
        y_mean = np.average(y_true_arr, weights=sample_weight_arr)
        ss_res: float = float(np.sum(sample_weight_arr * (y_true_arr - y_pred_arr) ** 2))
        ss_tot: float = float(np.sum(sample_weight_arr * (y_true_arr - y_mean) ** 2))
    else:
        y_mean = np.mean(y_true_arr)
        ss_res = float(np.sum((y_true_arr - y_pred_arr) ** 2))
        ss_tot = float(np.sum((y_true_arr - y_mean) ** 2))

    if ss_tot == 0:
        return 1.0 if ss_res == 0 else 0.0

    return float(1 - (ss_res / ss_tot))


def mean_absolute_percentage_error(
    y_true: Union[np.ndarray, list],
    y_pred: Union[np.ndarray, list],
    epsilon: float = 1e-10,
) -> float:
    """
    Compute Mean Absolute Percentage Error.

    MAPE = mean(|y_true - y_pred| / |y_true|) * 100

    Args:
        y_true: True values
        y_pred: Predicted values
        epsilon: Small value to avoid division by zero

    Returns:
        MAPE value in percentage (lower is better)

    Note:
        MAPE can be misleading when y_true contains values near zero.

    """
    y_true_arr: np.ndarray = np.asarray(y_true)
    y_pred_arr: np.ndarray = np.asarray(y_pred)

    # Avoid division by zero
    denominator = np.maximum(np.abs(y_true_arr), epsilon)
    percentage_errors = np.abs((y_true_arr - y_pred_arr) / denominator)

    return float(np.mean(percentage_errors) * 100)


def max_error(
    y_true: Union[np.ndarray, list],
    y_pred: Union[np.ndarray, list],
) -> float:
    """
    Compute Maximum Absolute Error.

    Args:
        y_true: True values
        y_pred: Predicted values

    Returns:
        Maximum absolute error

    """
    y_true_arr: np.ndarray = np.asarray(y_true)
    y_pred_arr: np.ndarray = np.asarray(y_pred)

    absolute_errors = np.abs(y_true_arr - y_pred_arr)

    return float(np.max(absolute_errors))


def compute_all_metrics(
    y_true: Union[np.ndarray, list],
    y_pred: Union[np.ndarray, list],
    sample_weight: Optional[np.ndarray] = None,
) -> dict[str, float]:
    """
    Compute all metrics at once.

    Args:
        y_true: True values
        y_pred: Predicted values
        sample_weight: Optional sample weights

    Returns:
        Dictionary with all metrics:
            - rmse
            - mae
            - r2
            - mape
            - max_error

    Example:
        >>> metrics = compute_all_metrics([1, 2, 3], [1.1, 2.0, 2.9])
        >>> print(metrics)
        {'rmse': 0.0816, 'mae': 0.0666, 'r2': 0.985, ...}

    """
    return {
        'rmse': rmse(y_true, y_pred, sample_weight),
        'mae': mae(y_true, y_pred, sample_weight),
        'r2': r2_score(y_true, y_pred, sample_weight),
        'mape': mean_absolute_percentage_error(y_true, y_pred),
        'max_error': max_error(y_true, y_pred),
    }


def compare_methods(
    y_true: np.ndarray,
    predictions: dict[str, np.ndarray],
) -> dict[str, dict[str, float]]:
    """
    Compare multiple prediction methods.

    Args:
        y_true: True values
        predictions: Dict mapping method name to predictions

    Returns:
        Nested dict: {method: {metric: value}}

    Example:
        >>> y_true = np.array([1, 2, 3])
        >>> preds = {
        ...     "PM7": np.array([0.9, 1.8, 2.7]),
        ...     "PM7+delta": np.array([1.05, 2.1, 2.95])
        ... }
        >>> results = compare_methods(y_true, preds)
        >>> print(results["PM7+delta"]["rmse"])

    """
    results = {}

    for method_name, y_pred in predictions.items():
        results[method_name] = compute_all_metrics(y_true, y_pred)

    return results
