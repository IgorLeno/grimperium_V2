"""
Database analyzer for GRIMPERIUM.

Provides qualitative CSV analysis for molecular databases,
including status counts, quality grade distribution, outlier
detection, and pipeline health metrics.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class AnalysisReport:
    """Complete analysis report for a molecular database CSV."""

    csv_path: Path
    total_rows: int
    status_counts: dict[str, int]
    orphan_running: int
    ok_count: int
    cbs_suspect_count: int
    cbs_ok_count: int

    # Null fields in OK rows
    null_h298_pm7: int
    null_target_delta: int

    # Quality grade distribution (from conformer_details JSON)
    grade_distribution: dict[str, int]
    grade_ab_rate: float

    # H298_pm7 stats (OK rows)
    h298_pm7_mean: float | None
    h298_pm7_std: float | None
    h298_pm7_min: float | None
    h298_pm7_max: float | None
    h298_pm7_positive_count: int
    h298_pm7_extreme_count: int  # H298_pm7/nheavy outside [-80, +20]

    # target_delta stats (OK rows)
    target_delta_mean: float | None
    target_delta_std: float | None
    target_delta_min: float | None
    target_delta_max: float | None
    target_delta_outliers_200: int
    target_delta_outliers_500: int

    # Electronic descriptors (OK rows)
    homo_mean: float | None
    homo_std: float | None
    homo_suspicious: int  # |HOMO| > 100 eV
    lumo_mean: float | None
    lumo_suspicious: int
    gap_mean: float | None
    gap_negative_count: int  # gap < 0
    gap_suspicious: int  # gap > 100

    # Conformers
    conformers_zero_count: int
    conformers_one_count: int
    conformers_mean: float | None

    # Timing (dynamic thresholds)
    crest_time_mean: float | None
    crest_time_max: float | None
    crest_timeout_count: int
    mopac_timeout_count: int

    # Pipeline health
    reruns_distribution: dict[int, int]
    max_reruns_count: int  # rows with reruns == 3
    crest_failed_count: int
    mopac_failed_count: int

    # Alerts + outliers
    alerts: list[str]
    top_delta_outliers: list[dict[str, Any]]
    top_homo_suspicious: list[dict[str, Any]]


class DatabaseAnalyzer:
    """Analyzes a molecular database CSV for quality and pipeline health."""

    def __init__(self) -> None:
        pass

    def analyze(
        self,
        csv_path: Path,
        detail_dir: Path | None = None,
    ) -> AnalysisReport:
        """Run full analysis on a CSV file.

        Args:
            csv_path: Path to the database CSV.
            detail_dir: Optional path to conformer_details directory
                        containing per-molecule JSON files with quality_grade.

        Returns:
            AnalysisReport with all computed metrics.
        """
        df = pd.read_csv(csv_path, low_memory=False)
        alerts: list[str] = []

        total_rows = len(df)
        status_counts = self._status_counts(df)
        ok_mask = (
            df["status"] == "OK"
            if "status" in df.columns
            else pd.Series([False] * total_rows, index=df.index)
        )
        ok_df = df[ok_mask]
        ok_count = len(ok_df)

        # CBS quality flag
        cbs_suspect = 0
        cbs_ok = 0
        if "cbs_quality_flag" in ok_df.columns:
            cbs_suspect = int((ok_df["cbs_quality_flag"] == "SUSPECT").sum())
            cbs_ok = int((ok_df["cbs_quality_flag"] == "OK").sum())
        if cbs_suspect > 0:
            alerts.append(
                f"⚠️ {cbs_suspect} molecules with SUSPECT CBS reference "
                "(excluded from training)"
            )

        orphan_running = self._orphan_running(df)
        if orphan_running > 0:
            alerts.append(
                f"⚠️ {orphan_running} RUNNING molecules without timestamp (orphaned)"
            )

        # Null fields in OK rows
        null_h298 = (
            int(ok_df["H298_pm7"].isna().sum()) if "H298_pm7" in ok_df.columns else 0
        )
        null_delta = (
            int(ok_df["target_delta_kcalmol"].isna().sum())
            if "target_delta_kcalmol" in ok_df.columns
            else 0
        )
        if null_h298 > 0:
            alerts.append(f"⚠️ {null_h298} OK rows with null H298_pm7")
        if null_delta > 0:
            alerts.append(f"⚠️ {null_delta} OK rows with null target_delta")

        # Grade distribution
        grade_dist, grade_ab = self._grade_distribution(ok_df, detail_dir)
        if grade_dist and grade_ab < 0.67:
            alerts.append(f"⚠️ A+B grade rate {grade_ab:.1%} is below 67% threshold")

        # H298_pm7 stats
        h298 = self._h298_stats(ok_df)
        if h298["positive_count"] > 0:
            alerts.append(f"⚠️ {h298['positive_count']} OK rows with positive H298_pm7")
        if h298["extreme_count"] > 0:
            alerts.append(
                f"⚠️ {h298['extreme_count']} molecules with extreme "
                "H298_pm7/nheavy (outside [-80, +20] kcal/mol)"
            )

        # Target delta stats
        delta = self._target_delta_stats(ok_df)
        if delta["outliers_500"] > 0:
            alerts.append(
                f"🚨 {delta['outliers_500']} molecules with "
                "|target_delta| > 500 kcal/mol"
            )
        elif delta["outliers_200"] > 0:
            alerts.append(
                f"⚠️ {delta['outliers_200']} molecules with "
                "|target_delta| > 200 kcal/mol"
            )

        # Electronic descriptors
        elec = self._electronic_stats(ok_df)
        if elec["gap_negative"] > 0:
            alerts.append(
                f"🚨 {elec['gap_negative']} molecules with " "negative HOMO-LUMO gap"
            )
        if elec["homo_suspicious"] > 0:
            alerts.append(f"⚠️ {elec['homo_suspicious']} molecules with |HOMO| > 100 eV")

        # Conformers
        conf = self._conformer_stats(ok_df)

        # Timing
        timing = self._timing_stats(ok_df)
        if timing["crest_timeout"] > 0:
            alerts.append(f"⚠️ {timing['crest_timeout']} molecules hit CREST timeout")
        if timing["mopac_timeout"] > 0:
            alerts.append(f"⚠️ {timing['mopac_timeout']} molecules hit MOPAC timeout")

        # Pipeline health (scoped to OK rows to avoid counting PENDING/RUNNING)
        pipeline = self._pipeline_health(ok_df)
        if pipeline["max_reruns"] > 0:
            alerts.append(
                f"⚠️ {pipeline['max_reruns']} molecules with reruns=3 "
                "(permanent SKIP candidates)"
            )
        if pipeline["mopac_failed"] > 0:
            alerts.append(
                f"⚠️ {pipeline['mopac_failed']} rows with " "abnormal mopac_status"
            )

        # Top outliers
        top_delta = self._top_delta_outliers(ok_df, n=10)
        top_homo = self._top_homo_suspicious(ok_df, n=10)

        return AnalysisReport(
            csv_path=csv_path,
            total_rows=total_rows,
            status_counts=status_counts,
            orphan_running=orphan_running,
            ok_count=ok_count,
            cbs_suspect_count=cbs_suspect,
            cbs_ok_count=cbs_ok,
            null_h298_pm7=null_h298,
            null_target_delta=null_delta,
            grade_distribution=grade_dist,
            grade_ab_rate=grade_ab,
            h298_pm7_mean=h298["mean"],
            h298_pm7_std=h298["std"],
            h298_pm7_min=h298["min"],
            h298_pm7_max=h298["max"],
            h298_pm7_positive_count=h298["positive_count"],
            h298_pm7_extreme_count=h298["extreme_count"],
            target_delta_mean=delta["mean"],
            target_delta_std=delta["std"],
            target_delta_min=delta["min"],
            target_delta_max=delta["max"],
            target_delta_outliers_200=delta["outliers_200"],
            target_delta_outliers_500=delta["outliers_500"],
            homo_mean=elec["homo_mean"],
            homo_std=elec["homo_std"],
            homo_suspicious=elec["homo_suspicious"],
            lumo_mean=elec["lumo_mean"],
            lumo_suspicious=elec["lumo_suspicious"],
            gap_mean=elec["gap_mean"],
            gap_negative_count=elec["gap_negative"],
            gap_suspicious=elec["gap_suspicious"],
            conformers_zero_count=conf["zero"],
            conformers_one_count=conf["one"],
            conformers_mean=conf["mean"],
            crest_time_mean=timing["crest_mean"],
            crest_time_max=timing["crest_max"],
            crest_timeout_count=timing["crest_timeout"],
            mopac_timeout_count=timing["mopac_timeout"],
            reruns_distribution=pipeline["reruns_dist"],
            max_reruns_count=pipeline["max_reruns"],
            crest_failed_count=pipeline["crest_failed"],
            mopac_failed_count=pipeline["mopac_failed"],
            alerts=alerts,
            top_delta_outliers=top_delta,
            top_homo_suspicious=top_homo,
        )

    # ------------------------------------------------------------------
    # Private analysis methods
    # ------------------------------------------------------------------

    @staticmethod
    def _status_counts(df: pd.DataFrame) -> dict[str, int]:
        if "status" not in df.columns:
            return {}
        return dict(df["status"].value_counts())

    @staticmethod
    def _orphan_running(df: pd.DataFrame) -> int:
        """Count RUNNING rows without a timestamp (likely orphaned)."""
        if "status" not in df.columns:
            return 0
        running = df[df["status"] == "RUNNING"]
        if running.empty:
            return 0
        if "timestamp" not in df.columns:
            return len(running)
        return int(running["timestamp"].isna().sum())

    @staticmethod
    def _grade_distribution(
        ok_df: pd.DataFrame,
        detail_dir: Path | None,
    ) -> tuple[dict[str, int], float]:
        """Load quality grades from per-molecule JSON files."""
        if detail_dir is None or not detail_dir.exists():
            return {}, 0.0

        mol_ids: set[str] = set()
        if "mol_id" in ok_df.columns:
            mol_ids = set(ok_df["mol_id"].dropna().astype(str))

        grades: dict[str, int] = {}
        total_graded = 0

        for json_file in detail_dir.glob("*.json"):
            try:
                with open(json_file, encoding="utf-8") as f:
                    data = json.load(f)
                mol_id = data.get("mol_id", json_file.stem)
                if mol_ids and str(mol_id) not in mol_ids:
                    continue
                grade = data.get("quality_grade")
                if grade:
                    grade_str = str(grade).upper()
                    grades[grade_str] = grades.get(grade_str, 0) + 1
                    total_graded += 1
            except (json.JSONDecodeError, OSError):
                continue

        if total_graded == 0:
            return grades, 0.0

        ab_count = grades.get("A", 0) + grades.get("B", 0)
        return grades, ab_count / total_graded

    @staticmethod
    def _h298_stats(ok_df: pd.DataFrame) -> dict[str, Any]:
        col = "H298_pm7"
        result: dict[str, Any] = {
            "mean": None,
            "std": None,
            "min": None,
            "max": None,
            "positive_count": 0,
            "extreme_count": 0,
        }
        if col not in ok_df.columns or ok_df[col].dropna().empty:
            return result

        vals = ok_df[col].dropna()
        result["mean"] = float(vals.mean())
        result["std"] = float(vals.std())
        result["min"] = float(vals.min())
        result["max"] = float(vals.max())
        result["positive_count"] = int((vals > 0).sum())

        # Extreme: normalize by nheavy
        if "nheavy" in ok_df.columns:
            valid = ok_df[[col, "nheavy"]].dropna()
            valid = valid[valid["nheavy"] > 0]
            if not valid.empty:
                per_heavy = valid[col] / valid["nheavy"]
                result["extreme_count"] = int(
                    ((per_heavy < -80) | (per_heavy > 20)).sum()
                )

        return result

    @staticmethod
    def _target_delta_stats(ok_df: pd.DataFrame) -> dict[str, Any]:
        col = "target_delta_kcalmol"
        result: dict[str, Any] = {
            "mean": None,
            "std": None,
            "min": None,
            "max": None,
            "outliers_200": 0,
            "outliers_500": 0,
        }
        if col not in ok_df.columns or ok_df[col].dropna().empty:
            return result

        # Exclude SUSPECT CBS rows — their deltas are nonsensical (e.g. −145k)
        df = ok_df
        if "cbs_quality_flag" in df.columns:
            df = df[df["cbs_quality_flag"] == "OK"]
        vals = df[col].dropna()
        result["mean"] = float(vals.mean())
        result["std"] = float(vals.std())
        result["min"] = float(vals.min())
        result["max"] = float(vals.max())
        result["outliers_200"] = int((vals.abs() > 200).sum())
        result["outliers_500"] = int((vals.abs() > 500).sum())
        return result

    @staticmethod
    def _electronic_stats(ok_df: pd.DataFrame) -> dict[str, Any]:
        result: dict[str, Any] = {
            "homo_mean": None,
            "homo_std": None,
            "homo_suspicious": 0,
            "lumo_mean": None,
            "lumo_suspicious": 0,
            "gap_mean": None,
            "gap_negative": 0,
            "gap_suspicious": 0,
        }
        if "mopac_homo_ev" in ok_df.columns:
            vals = ok_df["mopac_homo_ev"].dropna()
            if not vals.empty:
                result["homo_mean"] = float(vals.mean())
                result["homo_std"] = float(vals.std())
                result["homo_suspicious"] = int((vals.abs() > 100).sum())

        if "mopac_lumo_ev" in ok_df.columns:
            vals = ok_df["mopac_lumo_ev"].dropna()
            if not vals.empty:
                result["lumo_mean"] = float(vals.mean())
                result["lumo_suspicious"] = int((vals.abs() > 100).sum())

        if "mopac_gap_ev" in ok_df.columns:
            vals = ok_df["mopac_gap_ev"].dropna()
            if not vals.empty:
                result["gap_mean"] = float(vals.mean())
                result["gap_negative"] = int((vals < 0).sum())
                result["gap_suspicious"] = int((vals > 100).sum())

        return result

    @staticmethod
    def _conformer_stats(ok_df: pd.DataFrame) -> dict[str, Any]:
        col = "crest_conformers_generated"
        result: dict[str, Any] = {"zero": 0, "one": 0, "mean": None}
        if col not in ok_df.columns:
            return result
        vals = ok_df[col].dropna()
        if vals.empty:
            return result
        result["zero"] = int((vals == 0).sum())
        result["one"] = int((vals == 1).sum())
        result["mean"] = float(vals.mean())
        return result

    @staticmethod
    def _timing_stats(ok_df: pd.DataFrame) -> dict[str, Any]:
        result: dict[str, Any] = {
            "crest_mean": None,
            "crest_max": None,
            "crest_timeout": 0,
            "mopac_timeout": 0,
        }
        if "crest_time_s" in ok_df.columns:
            vals = ok_df["crest_time_s"].dropna()
            if not vals.empty:
                result["crest_mean"] = float(vals.mean())
                result["crest_max"] = float(vals.max())
                # Dynamic threshold from assigned_crest_timeout (min → sec)
                if "assigned_crest_timeout" in ok_df.columns:
                    both = ok_df[["crest_time_s", "assigned_crest_timeout"]].dropna()
                    if not both.empty:
                        result["crest_timeout"] = int(
                            (
                                both["crest_time_s"]
                                >= both["assigned_crest_timeout"] * 60
                            ).sum()
                        )

        if (
            "mopac_time_s" in ok_df.columns
            and "assigned_mopac_timeout" in ok_df.columns
        ):
            both = ok_df[["mopac_time_s", "assigned_mopac_timeout"]].dropna()
            if not both.empty:
                result["mopac_timeout"] = int(
                    (both["mopac_time_s"] >= both["assigned_mopac_timeout"] * 60).sum()
                )

        return result

    @staticmethod
    def _pipeline_health(df: pd.DataFrame) -> dict[str, Any]:
        result: dict[str, Any] = {
            "reruns_dist": {},
            "max_reruns": 0,
            "crest_failed": 0,
            "mopac_failed": 0,
        }
        if "reruns" in df.columns:
            vals = df["reruns"].dropna().astype(int)
            result["reruns_dist"] = dict(vals.value_counts().sort_index())
            result["max_reruns"] = int((vals == 3).sum())

        if "crest_status" in df.columns:
            has_status = df["crest_status"].notna() & (df["crest_status"] != "")
            healthy = df["crest_status"].str.lower().isin({"ok", "success"})
            result["crest_failed"] = int((~healthy & has_status).sum())

        if "mopac_status" in df.columns:
            has_status = df["mopac_status"].notna() & (df["mopac_status"] != "")
            healthy = df["mopac_status"].str.lower() == "ok"
            result["mopac_failed"] = int((~healthy & has_status).sum())

        return result

    @staticmethod
    def _top_delta_outliers(
        ok_df: pd.DataFrame,
        n: int = 10,
    ) -> list[dict[str, Any]]:
        col = "target_delta_kcalmol"
        if col not in ok_df.columns:
            return []
        # Exclude SUSPECT CBS rows from outlier ranking
        filtered = ok_df
        if "cbs_quality_flag" in filtered.columns:
            filtered = filtered[filtered["cbs_quality_flag"] == "OK"]
        valid = filtered[filtered[col].notna()].copy()
        if valid.empty:
            return []
        valid["abs_delta"] = valid[col].abs()
        top = valid.nlargest(n, "abs_delta")
        results: list[dict[str, Any]] = []
        for _, row in top.iterrows():
            entry: dict[str, Any] = {"delta": float(row[col])}
            if "mol_id" in ok_df.columns:
                entry["mol_id"] = str(row["mol_id"])
            if "smiles" in ok_df.columns:
                entry["smiles"] = str(row["smiles"])
            if "nheavy" in ok_df.columns and pd.notna(row["nheavy"]):
                entry["nheavy"] = int(row["nheavy"])
            results.append(entry)
        return results

    @staticmethod
    def _top_homo_suspicious(
        ok_df: pd.DataFrame,
        n: int = 10,
    ) -> list[dict[str, Any]]:
        col = "mopac_homo_ev"
        if col not in ok_df.columns:
            return []
        valid = ok_df[ok_df[col].notna()].copy()
        suspicious = valid[valid[col].abs() > 100]
        if suspicious.empty:
            return []
        top = suspicious.head(n)
        results: list[dict[str, Any]] = []
        for _, row in top.iterrows():
            entry: dict[str, Any] = {"homo_ev": float(row[col])}
            if "mol_id" in ok_df.columns:
                entry["mol_id"] = str(row["mol_id"])
            if "smiles" in ok_df.columns:
                entry["smiles"] = str(row["smiles"])
            results.append(entry)
        return results
