"""Per-arm metrics and the decision rule fixed before the experiment ran.

Only the latest line of each (mol_id, arm) counts. The delta
``H298_cbs - dHf_PM7`` uses verified minima only; saddle, unverified,
failed and error outcomes are counted separately and never enter it.

Comparison sets (paired, so no arm is judged on easier molecules):

* delta: molecules whose every arm produced a verified delta;
* recovery: molecules whose every arm produced an RMSD.

Decision rule (spec section 5): among arms whose recovery fraction
(``min_rmsd < 0.5 A``) is at most 5 percentage points below arm A's,
take the lowest delta standard deviation; arms within 0.1 kcal/mol of
that lowest are tied and the cheapest wins. Cost = mean selection plus
MOPAC wall time per molecule (CREST is shared by all arms); a remaining
tie goes to the earlier arm in A, B, C order.
"""

from __future__ import annotations

import csv
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from conformer_validation.arms import ARM_NAMES
from conformer_validation.metrics import (
    RECOVERY_THRESHOLD_ANGSTROM,
    median,
    median_absolute_deviation,
    sample_std,
)

RECOVERY_MARGIN = 0.05
TIE_KCAL_MOL = 0.1
OUTCOME_CLASSES = ("verified", "saddle", "unverified", "failed", "error")


@dataclass(frozen=True)
class ArmSummary:
    """One line of ``summary.csv``."""

    arm: str
    molecules: int
    verified: int
    saddle: int
    unverified: int
    failed: int
    error: int
    delta_n_all: int
    delta_std_all: float | None
    delta_mad_all: float | None
    delta_n_common: int
    delta_std_common: float | None
    delta_mad_common: float | None
    recovery_n_common: int
    recovery_fraction_common: float | None
    final_rmsd_median: float | None
    mean_overlap_with_a: float | None
    identical_to_a: int
    optimizations_total: int
    force_total: int
    selection_seconds_total: float
    mopac_seconds_total: float
    cost_seconds_per_molecule: float | None


@dataclass(frozen=True)
class Decision:
    """Outcome of the decision rule, with the reasoning spelled out."""

    winner: str | None
    eligible: tuple[str, ...]
    tied: tuple[str, ...]
    explanation: str


def latest_rows(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep the last line written for every (mol_id, arm)."""
    latest: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        latest[(row["mol_id"], row["arm"])] = row
    return list(latest.values())


def outcome_class(row: dict[str, Any]) -> str:
    """verified / saddle / unverified / failed, or error."""
    if row.get("error"):
        return "error"
    return str(row["final_class"])


def _common(
    by_arm: dict[str, dict[str, dict[str, Any]]], arms: Sequence[str], field: str
) -> set[str]:
    sets = [
        {mol for mol, row in by_arm[arm].items() if row.get(field) is not None}
        for arm in arms
    ]
    return set.intersection(*sets) if sets else set()


def summarize(rows: Iterable[dict[str, Any]]) -> tuple[list[ArmSummary], Decision]:
    """Per-arm summaries and the decision over the present arms."""
    latest = latest_rows(rows)
    arms = [arm for arm in ARM_NAMES if any(row["arm"] == arm for row in latest)]
    by_arm: dict[str, dict[str, dict[str, Any]]] = {arm: {} for arm in arms}
    for row in latest:
        by_arm[row["arm"]][row["mol_id"]] = row

    common_delta = _common(by_arm, arms, "delta_kcal_mol")
    common_rmsd = _common(by_arm, arms, "min_rmsd")
    summaries: list[ArmSummary] = []
    for arm in arms:
        arm_rows = by_arm[arm]
        ok = [row for row in arm_rows.values() if not row.get("error")]
        classes = Counter(outcome_class(row) for row in arm_rows.values())
        deltas = [
            float(row["delta_kcal_mol"])
            for row in ok
            if row.get("delta_kcal_mol") is not None
        ]
        deltas_common = [
            float(arm_rows[mol]["delta_kcal_mol"]) for mol in sorted(common_delta)
        ]
        recovered = [
            float(arm_rows[mol]["min_rmsd"]) < RECOVERY_THRESHOLD_ANGSTROM
            for mol in sorted(common_rmsd)
        ]
        overlaps = [
            (int(row["overlap_with_A"]), len(row["selected_indices"]))
            for row in ok
            if row.get("overlap_with_A") is not None
        ]
        selection = sum(float(row["selection_seconds"]) for row in ok)
        mopac = sum(float(row["mopac_seconds"]) for row in ok)
        summaries.append(
            ArmSummary(
                arm=arm,
                molecules=len(arm_rows),
                verified=classes["verified"],
                saddle=classes["saddle"],
                unverified=classes["unverified"],
                failed=classes["failed"],
                error=classes["error"],
                delta_n_all=len(deltas),
                delta_std_all=sample_std(deltas),
                delta_mad_all=median_absolute_deviation(deltas),
                delta_n_common=len(deltas_common),
                delta_std_common=sample_std(deltas_common),
                delta_mad_common=median_absolute_deviation(deltas_common),
                recovery_n_common=len(recovered),
                recovery_fraction_common=(
                    sum(recovered) / len(recovered) if recovered else None
                ),
                final_rmsd_median=median(
                    [
                        float(row["final_rmsd"])
                        for row in ok
                        if row.get("final_rmsd") is not None
                    ]
                ),
                mean_overlap_with_a=(
                    sum(o for o, _ in overlaps) / len(overlaps) if overlaps else None
                ),
                identical_to_a=sum(o == size for o, size in overlaps),
                optimizations_total=sum(int(row["n_optimizations"]) for row in ok),
                force_total=sum(int(row["n_force"]) for row in ok),
                selection_seconds_total=selection,
                mopac_seconds_total=mopac,
                cost_seconds_per_molecule=(selection + mopac) / len(ok) if ok else None,
            )
        )
    return summaries, decide(summaries)


def decide(summaries: Sequence[ArmSummary]) -> Decision:
    """Apply the pre-registered decision rule mechanically."""
    by_arm = {summary.arm: summary for summary in summaries}
    baseline = by_arm.get("A")
    if baseline is None or baseline.recovery_fraction_common is None:
        return Decision(None, (), (), "Sem decisão: braço A sem dados de recuperação.")
    floor = baseline.recovery_fraction_common - RECOVERY_MARGIN
    eligible = tuple(
        summary.arm
        for summary in summaries
        if summary.recovery_fraction_common is not None
        and summary.recovery_fraction_common >= floor - 1e-12
        and summary.delta_std_common is not None
    )
    if not eligible:
        return Decision(
            None, (), (), "Sem decisão: nenhum braço elegível tem DP do delta."
        )
    best = min(by_arm[arm].delta_std_common or 0.0 for arm in eligible)
    tied = tuple(
        arm
        for arm in eligible
        if (by_arm[arm].delta_std_common or 0.0) <= best + TIE_KCAL_MOL + 1e-12
    )
    winner = min(
        tied,
        key=lambda arm: (
            by_arm[arm].cost_seconds_per_molecule or 0.0,
            ARM_NAMES.index(arm),
        ),
    )
    excluded = [s.arm for s in summaries if s.arm not in eligible]
    explanation = (
        f"Piso de recuperação {floor:.3f} (A {baseline.recovery_fraction_common:.3f} "
        f"- {RECOVERY_MARGIN:.2f}); elegíveis {list(eligible)}"
        + (f", excluídos {excluded}" if excluded else "")
        + f". Menor DP do delta {best:.3f} kcal/mol; dentro de {TIE_KCAL_MOL} "
        f"kcal/mol: {list(tied)}. "
        + (
            f"Empate: vence o mais barato, {winner}."
            if len(tied) > 1
            else f"Vencedor: {winner}."
        )
    )
    return Decision(winner, eligible, tied, explanation)


def _fmt(value: object) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def render_report(summaries: Sequence[ArmSummary], decision: Decision) -> str:
    """Markdown report of the summaries and the decision."""
    columns = [
        ("arm", "Braço"),
        ("molecules", "Moléculas"),
        ("verified", "verified"),
        ("saddle", "saddle"),
        ("unverified", "unverified"),
        ("failed", "failed"),
        ("error", "erro"),
        ("delta_n_common", "n delta (comum)"),
        ("delta_std_common", "DP delta (comum)"),
        ("delta_mad_common", "MAD delta (comum)"),
        ("delta_std_all", "DP delta (todas)"),
        ("recovery_fraction_common", "Recuperação < 0,5 Å"),
        ("final_rmsd_median", "RMSD final (mediana)"),
        ("mean_overlap_with_a", "|X∩A| médio"),
        ("optimizations_total", "Otimizações"),
        ("force_total", "FORCE"),
        ("cost_seconds_per_molecule", "Custo s/molécula"),
    ]
    lines = [
        "# Validação da seleção de conformeros",
        "",
        "Gerado por `scripts/conformer_validation summarize`. Delta = "
        "`H298_cbs − ΔHf_PM7`, só mínimos verificados; DP amostral; MAD = "
        'desvio absoluto mediano. Conjuntos "comum": moléculas com valor em '
        "todos os braços.",
        "",
        "| " + " | ".join(title for _, title in columns) + " |",
        "|" + "---|" * len(columns),
    ]
    for summary in summaries:
        payload = asdict(summary)
        lines.append("| " + " | ".join(_fmt(payload[key]) for key, _ in columns) + " |")
    lines += [
        "",
        "## Decisão",
        "",
        f"**{decision.winner or 'sem decisão'}** — {decision.explanation}",
        "",
    ]
    return "\n".join(lines)


def write_outputs(
    summaries: Sequence[ArmSummary], decision: Decision, out_dir: Path
) -> tuple[Path, Path]:
    """Write ``summary.csv`` and ``report.md`` into ``out_dir``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        fields = list(ArmSummary.__dataclass_fields__)
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for summary in summaries:
            writer.writerow(asdict(summary))
    report_path = out_dir / "report.md"
    report_path.write_text(render_report(summaries, decision), encoding="utf-8")
    return summary_path, report_path


__all__ = [
    "OUTCOME_CLASSES",
    "RECOVERY_MARGIN",
    "TIE_KCAL_MOL",
    "ArmSummary",
    "Decision",
    "decide",
    "latest_rows",
    "outcome_class",
    "render_report",
    "summarize",
    "write_outputs",
]
