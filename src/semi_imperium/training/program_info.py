"""Which CREST, xtb and MOPAC produced a number, recorded with the number.

Same version string does not mean same binary (a Fedora RPM, a static build
and a conda-forge build all say ``xtb 6.7.1``), so each program is recorded
by its reported version *and* the sha256 of the file that actually ran.
A program that cannot be found or does not report a version is recorded as
such, never omitted: the run itself decides whether it can proceed.

The binaries' hashes miss the shared libraries they load (``libmopac.so``,
openblas), so the conda environment around them is recorded too, by its
package list (``conda-meta``: name, version and build of each package).
"""

from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from collections.abc import Callable, Mapping
from pathlib import Path

ProgramProbe = Callable[[Mapping[str, str]], dict[str, str]]

UNAVAILABLE = "unavailable"
UNKNOWN = "unknown"

# crest: "Version 3.0.2, Tue, ..."; xtb: "xtb version 6.7.1 (edcfbbe)";
# mopac: "MOPAC version 23.2.2".
_VERSION = re.compile(r"\bversion\s+(\d+(?:\.\d+)*)", re.IGNORECASE)
_VERSION_TIMEOUT_S = 30.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _reported_version(executable: str) -> str:
    try:
        completed = subprocess.run(  # noqa: S603  (executable from run settings)
            [executable, "--version"],
            capture_output=True,
            text=True,
            errors="replace",
            timeout=_VERSION_TIMEOUT_S,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return UNKNOWN
    match = _VERSION.search(completed.stdout + completed.stderr)
    return match.group(1) if match else UNKNOWN


def _environment_sha256(binaries: list[Path]) -> str:
    """Hash of the package list of the conda environment(s) holding ``binaries``.

    Built from file names only (``<name>-<version>-<build>.json``), never from
    their contents, which embed the install prefix: the same lockfile in
    another home folder gives the same hash.
    """
    packages: set[str] = set()
    for binary in binaries:
        meta = binary.parent.parent / "conda-meta"
        if meta.is_dir():
            packages.update(p.name for p in meta.glob("*.json"))
    if not packages:
        return UNAVAILABLE
    return hashlib.sha256("\n".join(sorted(packages)).encode()).hexdigest()


def probe_programs(executables: Mapping[str, str]) -> dict[str, str]:
    """Map each program name to its version, plus ``<name>_sha256`` and ``<name>_path``.

    ``environment_sha256`` names the conda environment of the programs
    (``unavailable`` outside one).
    """
    info: dict[str, str] = {}
    found: list[Path] = []
    for name, executable in executables.items():
        resolved = shutil.which(executable)
        if resolved is None:
            info[name] = UNAVAILABLE
            info[f"{name}_path"] = executable
            continue
        path = Path(resolved)
        found.append(path)
        info[name] = _reported_version(resolved)
        info[f"{name}_sha256"] = _sha256(path)
        info[f"{name}_path"] = str(path)
    info["environment_sha256"] = _environment_sha256(found)
    return info


def program_mismatches(
    central: Mapping[str, str], cell: Mapping[str, str]
) -> list[str]:
    """Why ``cell``'s programs differ from ``central``'s; empty when they match.

    Versions and sha256 must be equal; install paths (``*_path``) may differ.
    A program missing or unreadable on either side never matches, so a
    broken central refuses every cell instead of admitting broken ones.
    """
    differences = []
    for key in sorted(set(central) | set(cell)):
        if key.endswith("_path"):
            continue
        ours, theirs = central.get(key), cell.get(key)
        if ours == theirs and ours not in (UNAVAILABLE, UNKNOWN):
            continue
        differences.append(
            f"{key}: {theirs or 'ausente'} (central {ours or 'ausente'})"
        )
    return differences
