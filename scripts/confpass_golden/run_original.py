"""Run the original CONFPASS PART 1 on the golden inputs and record its output.

Runs in the *isolated* CONFPASS environment, never in the project one:
it must stay Python 3.8 compatible and imports CONFPASS from its source
checkout. For every SDF listed in ``inputs.json`` it records:

* the priority list of ``pipe_x_as`` (CONFPASS default, x=0.8, x_as=0.2)
  and of its two components ``pipe_x`` and ``pipe_as``, 0-based in SDF
  order (the CLI prints them 1-based);
* the intermediate the port is easiest to get wrong: the rotatable-bond
  dihedrals CONFPASS kept (column labels, CONFPASS's 1-based atom names,
  sorted: CONFPASS builds them through ``set`` iteration, so their order
  follows ``PYTHONHASHSEED``) and the clusters at
  ``n = round(0.8 * conformers)``;
* any exception CONFPASS raises, by type and message, per method.

The clustering is computed once per case and shared by the three methods,
exactly as ``confpass.conp.get_priority`` builds it.

Usage::

    <isolated-env>/bin/python scripts/confpass_golden/run_original.py \\
        --confpass-src ~/Estagio/confpass-original/src/confpass \\
        --fixtures tests/fixtures/confpass_golden
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from importlib import import_module
from pathlib import Path
from typing import Any

X = 0.8
X_AS = 0.2
N_NTH = 3
METHODS = ("pipe_x_as", "pipe_x", "pipe_as")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plain(value: Any) -> Any:
    """Convert numpy scalars and nested sequences to JSON primitives."""
    if isinstance(value, (list, tuple)):  # noqa: UP038 - must run on Python 3.8
        return [plain(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return value


def descriptor_columns(modules: dict[str, Any], sdf: str) -> list[str]:
    """Repeat the descriptor pipeline of ``get_cluster_df`` up to its columns."""
    selected = modules["isolate"].isolate_dihedral(sdf)
    raw = modules["cal"].dihedral_df(sdf, selected[0], selected[1])
    fixed = modules["param"].remove_fixed_dihedrals(raw, selected[0], selected[1])[2]
    final = modules["param"].remove_fixed_bond_df(raw, fixed)[0]
    corrected = modules["correct"].correction_by_gap(final)[0]
    return sorted(str(column) for column in corrected.columns)


def failure(exc: Exception) -> dict[str, Any]:
    return {"error": {"type": type(exc).__name__, "message": str(exc)}}


def run_case(modules: dict[str, Any], sdf: Path) -> dict[str, Any]:
    """CONFPASS's own calls first, so a recorded failure is CONFPASS's."""
    path = str(sdf)
    try:
        clustering = modules["clustering"].get_cluster_df(path)
    except Exception as exc:  # noqa: BLE001 - the failure itself is the golden output
        return {"clustering": failure(exc)}

    priority: dict[str, Any] = {}
    for method in METHODS:
        try:
            getter = modules["priority"].GetPriority(clustering)
            getter.priority_df(x_=X, x_de_=X_AS, n_=N_NTH, method_=method)
            priority[method] = plain(getter.priority_ls_df["priority_ls"].tolist()[0])
        except Exception as exc:  # noqa: BLE001 - recorded as the golden output
            priority[method] = failure(exc)

    result: dict[str, Any] = {
        "descriptor_columns": descriptor_columns(modules, path),
        "clustering_rows": len(clustering),
        "priority": priority,
    }
    n_x = round(len(clustering) * X)
    at_x = clustering[clustering["n"] == n_x]["clusters"].tolist()
    if at_x:
        result["clusters_at_x"] = {
            "n_clusters": n_x,
            "members": [plain(members) for _label, members in at_x[0]],
        }
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--confpass-src", type=Path, required=True)
    parser.add_argument("--fixtures", type=Path, required=True)
    args = parser.parse_args(argv)

    source: Path = args.confpass_src.expanduser().resolve()
    sys.path.insert(0, str(source))
    modules = {
        "clustering": import_module("clustering_dih_v7"),
        "priority": import_module("GetPriority_v3"),
        "isolate": import_module("isolate_key_dihedral_v5"),
        "cal": import_module("cal_dihedral_v2"),
        "param": import_module("dihedral_parameter_v2"),
        "correct": import_module("correcting_dihedral_v1"),
    }
    versions = {
        name: import_module(name).__version__
        for name in ("numpy", "pandas", "sklearn", "rdkit")
    }

    fixtures: Path = args.fixtures
    inputs = json.loads((fixtures / "inputs.json").read_text())
    cases = []
    for case in inputs["cases"]:
        sdf = fixtures / "inputs" / case["sdf"]
        if sha256_file(sdf) != case["sdf_sha256"]:
            raise SystemExit(f"{sdf} does not match inputs.json; rebuild the inputs")
        result = run_case(modules, sdf)
        cases.append(dict(case=case["case"], sdf_sha256=case["sdf_sha256"], **result))
        errors = [
            name
            for name, value in result.get("priority", result).items()
            if isinstance(value, dict) and "error" in value
        ]
        print(
            f"{case['case']:<18} {'errors: ' + ', '.join(errors) if errors else 'ok'}"
        )

    commit_file = source.parent / "COMMIT"
    golden = {
        "confpass": {
            "repository": "https://github.com/Goodman-lab/CONFPASS",
            "commit": commit_file.read_text().strip() if commit_file.exists() else None,
            "license": "MIT",
            "sources_sha256": {
                path.name: sha256_file(path) for path in sorted(source.glob("*.py"))
            },
        },
        "environment": dict(python=platform.python_version(), **versions),
        "parameters": {"x": X, "x_as": X_AS},
        "indexing": "0-based positions in SDF record order",
        "cases": cases,
    }
    (fixtures / "golden.json").write_text(json.dumps(golden, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
