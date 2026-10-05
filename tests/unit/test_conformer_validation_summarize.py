"""Per-arm metrics and the pre-registered decision rule, on synthetic rows."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pytest

from tests.unit.conformer_validation_support import load

summarize = load("summarize")
metrics = load("metrics")

MOLECULES = [f"m{i}" for i in range(6)]


def row(
    mol_id: str,
    arm: str,
    *,
    delta: float | None,
    rmsd: float | None = 0.2,
    final_class: str = "verified",
    seconds: float = 1.0,
    error: str | None = None,
) -> dict[str, Any]:
    if error is not None:
        return {"mol_id": mol_id, "arm": arm, "error": error, "error_code": "x"}
    return {
        "mol_id": mol_id,
        "arm": arm,
        "error": None,
        "final_class": final_class,
        "delta_kcal_mol": delta if final_class == "verified" else None,
        "min_rmsd": rmsd,
        "final_rmsd": rmsd,
        "selected_indices": list(range(10)),
        "overlap_with_A": 10 if arm == "A" else 8,
        "n_optimizations": 10,
        "n_force": 10,
        "selection_seconds": seconds / 2,
        "mopac_seconds": seconds / 2,
    }


def arm_rows(
    arm: str,
    deltas: list[float],
    *,
    rmsds: list[float] | None = None,
    seconds: float = 1.0,
) -> list[dict[str, Any]]:
    rmsds = rmsds or [0.2] * len(deltas)
    return [
        row(mol, arm, delta=d, rmsd=r, seconds=seconds)
        for mol, d, r in zip(MOLECULES, deltas, rmsds, strict=True)
    ]


SPREAD_A = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
SPREAD_B = [3.0, 4.0, 5.0, 5.0, 6.0, 7.0]


def by_arm(summaries: list[Any]) -> dict[str, Any]:
    return {summary.arm: summary for summary in summaries}


def test_delta_statistics_use_verified_rows_only() -> None:
    rows = arm_rows("A", SPREAD_A)
    rows[0] = row("m0", "A", delta=None, final_class="saddle")
    rows[1] = row("m1", "A", delta=None, error="boom")

    (summary,), _ = summarize.summarize(rows)

    kept = SPREAD_A[2:]
    assert summary.molecules == 6
    assert (summary.verified, summary.saddle, summary.error) == (4, 1, 1)
    assert summary.delta_n_all == 4
    assert summary.delta_std_all == pytest.approx(metrics.sample_std(kept))
    assert summary.delta_mad_all == pytest.approx(
        metrics.median_absolute_deviation(kept)
    )
    assert summary.optimizations_total == 50


def test_latest_line_of_a_pair_wins() -> None:
    first = row("m0", "A", delta=None, error="boom")
    retried = row("m0", "A", delta=1.0)

    (summary,), _ = summarize.summarize([first, retried])

    assert summary.error == 0 and summary.verified == 1


def test_common_set_pairs_molecules_across_arms() -> None:
    rows = arm_rows("A", SPREAD_A) + arm_rows("B", SPREAD_B)
    rows[6] = row("m0", "B", delta=None, final_class="unverified")

    summaries, _ = summarize.summarize(rows)

    a, b = by_arm(summaries)["A"], by_arm(summaries)["B"]
    assert a.delta_n_common == b.delta_n_common == 5
    assert a.delta_std_common == pytest.approx(metrics.sample_std(SPREAD_A[1:]))
    assert a.delta_n_all == 6


def test_lowest_deviation_wins_when_recovery_holds() -> None:
    rows = arm_rows("A", SPREAD_A) + arm_rows("B", SPREAD_B)

    _, decision = summarize.summarize(rows)

    assert decision.winner == "B"
    assert decision.eligible == ("A", "B")


def test_arm_losing_too_much_recovery_is_not_eligible() -> None:
    # A recovers 6/6; B only 5/6 (-16.7 pp), beyond the 5 pp margin.
    rows = arm_rows("A", SPREAD_A) + arm_rows(
        "B", SPREAD_B, rmsds=[0.2, 0.2, 0.2, 0.2, 0.2, 0.9]
    )

    summaries, decision = summarize.summarize(rows)

    assert by_arm(summaries)["B"].recovery_fraction_common == pytest.approx(5 / 6)
    assert decision.eligible == ("A",)
    assert decision.winner == "A"


def test_recovery_drop_of_exactly_five_points_is_still_eligible() -> None:
    a = [{"min_rmsd": 0.2}] * 20
    b = [{"min_rmsd": 0.2}] * 19 + [{"min_rmsd": 0.9}]
    rows = []
    for i, (ra, rb) in enumerate(zip(a, b, strict=True)):
        rows.append(row(f"x{i}", "A", delta=float(i), rmsd=ra["min_rmsd"]))
        rows.append(row(f"x{i}", "B", delta=float(i) / 2, rmsd=rb["min_rmsd"]))

    _, decision = summarize.summarize(rows)

    assert decision.eligible == ("A", "B")
    assert decision.winner == "B"


def test_tie_within_a_tenth_goes_to_the_cheapest_arm() -> None:
    shifted = [value + 0.05 * (i % 2) for i, value in enumerate(SPREAD_A)]
    rows = (
        arm_rows("A", SPREAD_A, seconds=5.0)
        + arm_rows("B", shifted, seconds=1.0)
        + arm_rows("C", SPREAD_B, seconds=9.0)
    )
    summaries, decision = summarize.summarize(rows)
    stds = {s.arm: s.delta_std_common for s in summaries}
    assert stds["C"] < stds["A"]  # C is clearly best: no tie with it

    assert decision.winner == "C"

    rows = arm_rows("A", SPREAD_A, seconds=5.0) + arm_rows("B", shifted, seconds=1.0)
    summaries, decision = summarize.summarize(rows)
    stds = {s.arm: s.delta_std_common for s in summaries}
    assert abs(stds["A"] - stds["B"]) <= summarize.TIE_KCAL_MOL
    assert set(decision.tied) == {"A", "B"}
    assert decision.winner == "B"  # cheaper, even if its SD is not the lowest


def test_exact_cost_tie_goes_to_the_earlier_arm() -> None:
    rows = arm_rows("A", SPREAD_A) + arm_rows("B", SPREAD_A)

    _, decision = summarize.summarize(rows)

    assert decision.tied == ("A", "B")
    assert decision.winner == "A"


def test_no_decision_without_arm_a() -> None:
    _, decision = summarize.summarize(arm_rows("B", SPREAD_B))

    assert decision.winner is None


def test_outputs_are_written(tmp_path: Path) -> None:
    summaries, decision = summarize.summarize(
        arm_rows("A", SPREAD_A) + arm_rows("B", SPREAD_B)
    )

    summary_path, report_path = summarize.write_outputs(summaries, decision, tmp_path)

    with summary_path.open() as handle:
        lines = list(csv.DictReader(handle))
    assert [line["arm"] for line in lines] == ["A", "B"]
    report = report_path.read_text()
    assert "**B**" in report
    assert "| A |" in report
