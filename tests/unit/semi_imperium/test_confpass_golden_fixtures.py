"""Integrity of the CONFPASS golden files the backend port will be held to.

These tests do not run CONFPASS. They pin the fixtures themselves: the
SDF inputs still hash to what the original was run on, they still read
back through the production adapter in CREST order, and every recorded
outcome is either a full permutation of the ensemble or an explicit
failure.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from semi_imperium.conformers import read_sd_record
from semi_imperium.conformers.confpass import INDEX_FIELD

FIXTURES = Path(__file__).resolve().parents[2] / "fixtures" / "confpass_golden"
METHODS = ("pipe_x_as", "pipe_x", "pipe_as")


def load(name: str) -> dict[str, Any]:
    data: dict[str, Any] = json.loads((FIXTURES / name).read_text())
    return data


INPUTS = load("inputs.json")
GOLDEN = load("golden.json")
CASES = [case["case"] for case in INPUTS["cases"]]


def input_case(name: str) -> dict[str, Any]:
    return next(case for case in INPUTS["cases"] if case["case"] == name)


def golden_case(name: str) -> dict[str, Any]:
    return next(case for case in GOLDEN["cases"] if case["case"] == name)


def records(sdf: str) -> list[str]:
    return [block for block in sdf.split("$$$$\n") if block.strip()]


def test_golden_covers_exactly_the_inputs_and_pins_the_original() -> None:
    assert [case["case"] for case in GOLDEN["cases"]] == CASES
    assert GOLDEN["confpass"]["commit"] == "1b5efb69585ea1f51bedccfed1d9d07133c18a53"
    assert GOLDEN["parameters"] == {"x": 0.8, "x_as": 0.2}


@pytest.mark.parametrize("name", CASES)
def test_sdf_input_is_the_one_confpass_was_run_on(name: str) -> None:
    case = input_case(name)
    digest = hashlib.sha256((FIXTURES / "inputs" / case["sdf"]).read_bytes())

    assert digest.hexdigest() == case["sdf_sha256"]
    assert golden_case(name)["sdf_sha256"] == case["sdf_sha256"]


@pytest.mark.parametrize("name", CASES)
def test_sdf_reads_back_in_crest_order(name: str) -> None:
    case = input_case(name)
    blocks = records((FIXTURES / "inputs" / case["sdf"]).read_text())
    structures = [read_sd_record(block + "$$$$") for block in blocks]

    assert len(structures) == case["conformers"]
    assert [int(s.data[INDEX_FIELD]) for s in structures] == list(
        range(case["conformers"])
    )
    assert {s.geometry.atom_count for s in structures} == {case["atoms"]}
    assert len({s.topology for s in structures}) == 1


@pytest.mark.parametrize("name", CASES)
def test_every_outcome_is_a_full_ranking_or_an_explicit_failure(name: str) -> None:
    conformers = input_case(name)["conformers"]
    priority = golden_case(name)["priority"]

    assert set(priority) == set(METHODS)
    for outcome in priority.values():
        if isinstance(outcome, dict):
            assert set(outcome["error"]) == {"type", "message"}
        else:
            assert sorted(outcome) == list(range(conformers))


def test_ensemble_without_variable_dihedrals_is_recorded_as_a_failure() -> None:
    case = golden_case("methyl_acetate")

    assert case["descriptor_columns"] == []
    assert case["clustering_rows"] == 0
    assert all("error" in outcome for outcome in case["priority"].values())
