"""Fidelity of the CONFPASS PART 1 port to the original's golden files.

Every case in ``tests/fixtures/confpass_golden`` was run through the
original CONFPASS; the port must reproduce its kept dihedral columns,
its clusters at ``x = 0.8`` and the priority lists of all three methods,
or fail with an explicit code where the original crashed.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from semi_imperium.conformers import (
    ConformerBackendError,
    ConfPassCandidate,
    read_sd_record,
)
from semi_imperium.conformers.confpass import INDEX_FIELD
from semi_imperium.conformers.confpass_port import (
    METHODS,
    NO_VARIABLE_DIHEDRAL,
    PortedConfPass,
    analyse,
    dihedral_descriptor,
)

sklearn_cluster = pytest.importorskip("sklearn.cluster")

FIXTURES = Path(__file__).resolve().parents[2] / "fixtures" / "confpass_golden"
GOLDEN: dict[str, Any] = json.loads((FIXTURES / "golden.json").read_text())
CASES = {case["case"]: case for case in GOLDEN["cases"]}
X = GOLDEN["parameters"]["x"]
X_AS = GOLDEN["parameters"]["x_as"]
SUCCEEDED = [name for name, case in CASES.items() if case["clustering_rows"]]


def sd_records(name: str) -> list[str]:
    text = (FIXTURES / "inputs" / f"{name}.sdf").read_text()
    return [block + "$$$$" for block in text.split("$$$$\n") if block.strip()]


def candidates(name: str) -> list[ConfPassCandidate]:
    return [
        ConfPassCandidate(
            index=int(read_sd_record(record).data[INDEX_FIELD]), sd_record=record
        )
        for record in sd_records(name)
    ]


@pytest.mark.parametrize("name", SUCCEEDED)
def test_port_keeps_the_dihedrals_the_original_kept(name: str) -> None:
    analysis = analyse(sd_records(name))

    assert sorted(analysis.descriptor_columns) == CASES[name]["descriptor_columns"]
    assert analysis.conformer_count == CASES[name]["clustering_rows"]


@pytest.mark.parametrize("name", SUCCEEDED)
def test_port_reproduces_the_clusters_at_x(name: str) -> None:
    golden = CASES[name]["clusters_at_x"]
    clusters = analyse(sd_records(name)).clusters_at(X)

    assert len(clusters) == golden["n_clusters"]
    assert [list(cluster) for cluster in clusters] == golden["members"]


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("name", SUCCEEDED)
def test_port_reproduces_every_priority_list(name: str, method: str) -> None:
    analysis = analyse(sd_records(name))

    assert list(analysis.priority(method, x=X, x_as=X_AS)) == (
        CASES[name]["priority"][method]
    )


@pytest.mark.parametrize("name", SUCCEEDED)
def test_one_ward_tree_matches_a_fit_per_cluster_count(name: str) -> None:
    """clustering_dih_v7 refits sklearn once per k; the port cuts one tree."""
    records = sd_records(name)
    analysis = analyse(records)
    _columns, descriptor = dihedral_descriptor(records)

    for k in range(1, len(records) + 1):
        model = sklearn_cluster.AgglomerativeClustering(n_clusters=k).fit(descriptor)
        groups: dict[int, list[int]] = {}
        for position, label in enumerate(model.labels_):
            groups.setdefault(int(label), []).append(position)
        expected = tuple(sorted(tuple(group) for group in groups.values()))
        assert analysis.clusterings[k - 1] == expected


def test_backend_ranks_every_candidate_by_the_default_method() -> None:
    name = "triacetin"
    given = candidates(name)

    rankings = PortedConfPass().prioritize(given)

    by_priority = sorted(rankings, key=lambda ranking: ranking.priority)
    assert [r.index for r in by_priority] == CASES[name]["priority"]["pipe_x_as"]
    assert [r.priority for r in by_priority] == list(range(len(given)))
    assert all(r.pas_completeness_class is None for r in rankings)


def test_backend_maps_positions_back_to_candidate_indices() -> None:
    given = [
        ConfPassCandidate(index=100 + c.index, sd_record=c.sd_record)
        for c in candidates("butanol")
    ]

    rankings = PortedConfPass(method="pipe_x").prioritize(given)

    expected = [100 + i for i in CASES["butanol"]["priority"]["pipe_x"]]
    assert [r.index for r in sorted(rankings, key=lambda r: r.priority)] == expected


def test_ensemble_without_a_variable_dihedral_fails_with_its_own_code() -> None:
    with pytest.raises(ConformerBackendError) as failure:
        PortedConfPass().prioritize(candidates("methyl_acetate"))

    assert failure.value.code == NO_VARIABLE_DIHEDRAL


def test_a_single_conformer_cannot_be_clustered() -> None:
    with pytest.raises(ConformerBackendError) as failure:
        PortedConfPass().prioritize(candidates("butanol")[:1])

    assert failure.value.code == "confpass_too_few_conformers"


def test_unreadable_sd_record_is_reported() -> None:
    broken = ConfPassCandidate(index=0, sd_record="not an sd record")

    with pytest.raises(ConformerBackendError) as failure:
        PortedConfPass().prioritize([broken, broken])

    assert failure.value.code == "sdf_parse_failed"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"method": "random"}, "Unknown CONFPASS method"),
        ({"x": 0.0}, "x must be in"),
        ({"x_as": 1.5}, "x_as must be in"),
    ],
)
def test_backend_rejects_unusable_parameters(
    kwargs: dict[str, Any], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        PortedConfPass(**kwargs)
