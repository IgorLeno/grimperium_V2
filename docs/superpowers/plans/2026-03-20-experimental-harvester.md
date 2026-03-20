# Experimental Harvester Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone Python CLI tool that collects experimental ΔfH° (gas-phase standard enthalpy of formation) for a list of SMILES from NIST WebBook, NistChemPy, and PubChem, writing curated results to `output/experimental.csv`.

**Architecture:** Three-source cascade (NIST WebBook HTML scrape → NistChemPy library → PubChem REST JSON) orchestrated by a pipeline module; SMILES resolved to CAS/InChIKey via RDKit + PubChem; checkpoint/resume for long-running batches.

**Tech Stack:** Python 3.10+, requests, BeautifulSoup4, pandas, rdkit, nistchempy, tqdm, argparse, unittest.mock (tests only)

---

## File Map

| File | Responsibility |
|------|---------------|
| `experimental_harvester/pyproject.toml` | Package metadata, deps, tool config (mypy, ruff, black) |
| `experimental_harvester/requirements.txt` | Pip-installable deps for quick setup |
| `experimental_harvester/resolver.py` | SMILES → InChIKey (RDKit) → CAS (PubChem REST) |
| `experimental_harvester/sources/__init__.py` | Empty package marker |
| `experimental_harvester/sources/nist_api.py` | Fetch & parse NIST WebBook HTML for ΔfH° |
| `experimental_harvester/sources/nist_chempy.py` | NistChemPy library wrapper for ΔfH° |
| `experimental_harvester/sources/pubchem.py` | PubChem PUG REST for ΔfH° |
| `experimental_harvester/pipeline.py` | `HarvestResult` dataclass + cascading `harvest_one` + batch with checkpoint |
| `experimental_harvester/harvest.py` | CLI entry point (argparse), reads CSV, calls pipeline, writes output |
| `experimental_harvester/tests/__init__.py` | Empty package marker |
| `experimental_harvester/tests/test_resolver.py` | Unit tests for resolver (mocked requests) |
| `experimental_harvester/tests/test_nist_api.py` | Unit tests for nist_api (mocked HTML responses) |
| `experimental_harvester/tests/test_pubchem.py` | Unit tests for pubchem (mocked JSON responses) |
| `experimental_harvester/tests/test_pipeline.py` | Unit tests for pipeline cascading logic |

---

## Task 1: Project Scaffold

**Files:**
- Create: `experimental_harvester/pyproject.toml`
- Create: `experimental_harvester/requirements.txt`
- Create: `experimental_harvester/sources/__init__.py`
- Create: `experimental_harvester/tests/__init__.py`

- [ ] **Step 1: Create project root and pyproject.toml**

```bash
mkdir -p experimental_harvester/sources experimental_harvester/tests experimental_harvester/output
```

`pyproject.toml`:
```toml
[project]
name = "experimental-harvester"
version = "0.1.0"
requires-python = ">=3.10"
dependencies = [
    "requests>=2.31",
    "pandas>=2.0",
    "rdkit>=2023.9",
    "nistchempy>=0.3",
    "tqdm>=4.66",
    "beautifulsoup4>=4.12",
    "lxml>=4.9",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.4",
    "pytest-cov>=4.1",
    "black>=23.0",
    "ruff>=0.1",
    "mypy>=1.5",
    "types-requests>=2.31",
    "types-beautifulsoup4>=4.12",
]

[tool.black]
line-length = 88
target-version = ["py310"]

[tool.ruff]
line-length = 88
target-version = "py310"
select = ["E", "F", "W", "I", "N", "UP", "B", "C4", "SIM"]

[tool.mypy]
python_version = "3.10"
strict = true
ignore_missing_imports = true

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-v --tb=short"
```

- [ ] **Step 2: Create requirements.txt**

```
requests>=2.31
pandas>=2.0
rdkit>=2023.9
nistchempy>=0.3
tqdm>=4.66
beautifulsoup4>=4.12
lxml>=4.9
```

- [ ] **Step 3: Create empty package markers**

`sources/__init__.py` and `tests/__init__.py` — both empty files.

- [ ] **Step 4: Commit scaffold**

```bash
cd experimental_harvester
git init
git add .
git commit -m "chore: init experimental_harvester project scaffold"
```

---

## Task 2: resolver.py — SMILES → InChIKey → CAS

**Files:**
- Create: `experimental_harvester/resolver.py`
- Create: `experimental_harvester/tests/test_resolver.py`

**IMPORTANT:** All tests mock `requests.get` — NEVER make real network calls.

### 2a — Write and run failing tests

- [ ] **Step 1: Write failing tests**

`tests/test_resolver.py`:
```python
"""Tests for resolver.py — SMILES → InChIKey → CAS."""
import unittest
from unittest.mock import MagicMock, patch

import resolver


class TestSmilesToInchikey(unittest.TestCase):
    def test_valid_smiles_returns_inchikey(self):
        # Water: known InChIKey
        result = resolver.smiles_to_inchikey("O")
        self.assertIsNotNone(result)
        self.assertTrue(result.startswith("XLYOFNOQVPJJNP"))

    def test_invalid_smiles_returns_none(self):
        result = resolver.smiles_to_inchikey("NOT_A_SMILES!!!")
        self.assertIsNone(result)

    def test_empty_smiles_returns_none(self):
        result = resolver.smiles_to_inchikey("")
        self.assertIsNone(result)


class TestInchikeyToCas(unittest.TestCase):
    def test_cas_found(self):
        # Mock both PubChem calls: CID lookup + CAS lookup
        cid_response = MagicMock()
        cid_response.status_code = 200
        cid_response.json.return_value = {
            "PropertyTable": {"Properties": [{"CID": 962}]}
        }

        cas_response = MagicMock()
        cas_response.status_code = 200
        cas_response.json.return_value = {
            "Record": {
                "Section": [
                    {
                        "TOCHeading": "Names and Identifiers",
                        "Section": [
                            {
                                "TOCHeading": "Other Identifiers",
                                "Section": [
                                    {
                                        "TOCHeading": "CAS",
                                        "Information": [
                                            {"Value": {"StringWithMarkup": [{"String": "7732-18-5"}]}}
                                        ],
                                    }
                                ],
                            }
                        ],
                    }
                ]
            }
        }

        with patch("requests.get", side_effect=[cid_response, cas_response]):
            result = resolver.inchikey_to_cas("XLYOFNOQVPJJNP-UHFFFAOYSA-N")

        self.assertEqual(result, "7732-18-5")

    def test_cas_not_found_returns_none(self):
        not_found = MagicMock()
        not_found.status_code = 404
        not_found.json.side_effect = Exception("not found")

        with patch("requests.get", return_value=not_found):
            result = resolver.inchikey_to_cas("XXXXXXXXXXXXXXX-UHFFFAOYSA-N")

        self.assertIsNone(result)

    def test_no_cas_section_returns_none(self):
        cid_response = MagicMock()
        cid_response.status_code = 200
        cid_response.json.return_value = {
            "PropertyTable": {"Properties": [{"CID": 999}]}
        }
        empty_cas = MagicMock()
        empty_cas.status_code = 200
        empty_cas.json.return_value = {"Record": {"Section": []}}

        with patch("requests.get", side_effect=[cid_response, empty_cas]):
            result = resolver.inchikey_to_cas("XLYOFNOQVPJJNP-UHFFFAOYSA-N")

        self.assertIsNone(result)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd experimental_harvester
python -m pytest tests/test_resolver.py -v
```
Expected: ImportError or NameError — `resolver` module not found.

### 2b — Implement resolver.py

- [ ] **Step 3: Implement resolver.py**

```python
"""resolver.py — Convert SMILES to InChIKey and CAS number."""
import logging
import time
from typing import Any

import requests
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey, MolToInchi  # type: ignore[import]

logger = logging.getLogger(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
TIMEOUT = 10


def smiles_to_inchikey(smiles: str) -> str | None:
    """Convert SMILES to InChIKey using RDKit (local, no network)."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.debug("Invalid SMILES: %s", smiles)
        return None
    inchi = MolToInchi(mol)
    if inchi is None:
        return None
    return MolToInchiKey(inchi)


def _extract_cas_from_record(record: dict[str, Any]) -> str | None:
    """Recursively search PubChem record for CAS number."""
    for section in record.get("Section", []):
        heading = section.get("TOCHeading", "")
        if heading == "CAS":
            infos = section.get("Information", [])
            if infos:
                val = infos[0].get("Value", {})
                markup = val.get("StringWithMarkup", [])
                if markup:
                    return markup[0].get("String")
        # Recurse into nested sections
        result = _extract_cas_from_record(section)
        if result:
            return result
    return None


def inchikey_to_cas(inchikey: str) -> str | None:
    """Resolve InChIKey to CAS via PubChem REST API."""
    try:
        # Step 1: InChIKey → CID
        url1 = f"{PUBCHEM_BASE}/compound/inchikey/{inchikey}/property/MolecularFormula/JSON"
        r1 = requests.get(url1, timeout=TIMEOUT)
        if r1.status_code != 200:
            logger.debug("CID not found for InChIKey: %s", inchikey)
            return None
        props = r1.json().get("PropertyTable", {}).get("Properties", [])
        if not props:
            return None
        cid = props[0].get("CID")
        if not cid:
            return None

        time.sleep(0.5)

        # Step 2: CID → CAS
        url2 = f"{PUBCHEM_VIEW}/data/compound/{cid}/JSON/?heading=CAS"
        r2 = requests.get(url2, timeout=TIMEOUT)
        if r2.status_code != 200:
            logger.debug("CAS not found for CID: %s", cid)
            return None
        record = r2.json().get("Record", {})
        return _extract_cas_from_record(record)

    except Exception as exc:
        logger.debug("inchikey_to_cas failed: %s", exc)
        return None


def smiles_to_cas(smiles: str) -> str | None:
    """Convert SMILES to CAS number via InChIKey."""
    inchikey = smiles_to_inchikey(smiles)
    if inchikey is None:
        logger.debug("Could not get InChIKey for SMILES: %s", smiles)
        return None
    return inchikey_to_cas(inchikey)
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
python -m pytest tests/test_resolver.py -v
```
Expected: All 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add resolver.py tests/test_resolver.py
git commit -m "feat: add resolver — SMILES to InChIKey and CAS via RDKit + PubChem"
```

---

## Task 3: sources/nist_api.py — NIST WebBook HTML scraping

**Files:**
- Create: `experimental_harvester/sources/nist_api.py`
- Create: `experimental_harvester/tests/test_nist_api.py`

### 3a — Failing tests first

- [ ] **Step 1: Write failing tests**

`tests/test_nist_api.py`:
```python
"""Tests for sources/nist_api.py — NIST WebBook HTML scraping."""
import unittest
from unittest.mock import MagicMock, patch

from sources import nist_api

# Minimal HTML with experimental ΔfH° at 298 K in kJ/mol
NIST_HTML_EXPERIMENTAL_KJ = """
<html><body>
<h2>Gas phase thermochemistry data</h2>
<table>
<tr><th>Delta fH° gas</th><th>kJ/mol</th><th>Method</th><th>Reference</th></tr>
<tr><td>-241.826</td><td>kJ/mol</td><td>Review</td><td>CODATA 1989</td></tr>
</table>
</body></html>
"""

NIST_HTML_EXPERIMENTAL_KCAL = """
<html><body>
<h2>Gas phase thermochemistry data</h2>
<table>
<tr><th>Delta fH° gas</th><th>kcal/mol</th><th>Method</th><th>Reference</th></tr>
<tr><td>-57.798</td><td>kcal/mol</td><td>Calorimetry</td><td>Cox 1970</td></tr>
</table>
</body></html>
"""

NIST_HTML_CALCULATED = """
<html><body>
<h2>Gas phase thermochemistry data</h2>
<table>
<tr><th>Delta fH° gas</th><th>kJ/mol</th><th>Method</th><th>Reference</th></tr>
<tr><td>-241.0</td><td>kJ/mol</td><td>Calculated</td><td>Some DFT</td></tr>
</table>
</body></html>
"""

NIST_HTML_NOT_FOUND = """
<html><body><h1>Not found</h1></body></html>
"""


class TestFetchNist(unittest.TestCase):
    def _mock_response(self, html: str, status: int = 200) -> MagicMock:
        r = MagicMock()
        r.status_code = status
        r.text = html
        return r

    def test_returns_experimental_value_kj_converted(self):
        with patch("requests.get", return_value=self._mock_response(NIST_HTML_EXPERIMENTAL_KJ)):
            result = nist_api.fetch_nist("7732-18-5")
        self.assertIsNotNone(result)
        # -241.826 kJ/mol / 4.184 ≈ -57.80 kcal/mol
        self.assertAlmostEqual(result["hf_gas_kcal_mol"], -241.826 / 4.184, places=2)
        self.assertEqual(result["source"], "NIST_WebBook")

    def test_returns_experimental_value_kcal_unchanged(self):
        with patch("requests.get", return_value=self._mock_response(NIST_HTML_EXPERIMENTAL_KCAL)):
            result = nist_api.fetch_nist("7732-18-5")
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result["hf_gas_kcal_mol"], -57.798, places=3)

    def test_rejects_calculated_value(self):
        with patch("requests.get", return_value=self._mock_response(NIST_HTML_CALCULATED)):
            result = nist_api.fetch_nist("7732-18-5")
        self.assertIsNone(result)

    def test_cas_not_found_returns_none(self):
        with patch("requests.get", return_value=self._mock_response(NIST_HTML_NOT_FOUND, 404)):
            result = nist_api.fetch_nist("0000-00-0")
        self.assertIsNone(result)

    def test_network_error_returns_none(self):
        with patch("requests.get", side_effect=Exception("timeout")):
            result = nist_api.fetch_nist("7732-18-5")
        self.assertIsNone(result)
```

- [ ] **Step 2: Run to confirm failure**

```bash
python -m pytest tests/test_nist_api.py -v
```
Expected: ImportError.

### 3b — Implement nist_api.py

- [ ] **Step 3: Implement sources/nist_api.py**

```python
"""sources/nist_api.py — Fetch ΔfH° from NIST WebBook HTML (Source 1)."""
import logging
import time
from typing import Any

import requests
from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)

NIST_URL = "https://webbook.nist.gov/cgi/cbook.cgi"
HEADERS = {"User-Agent": "experimental_harvester/1.0 (academic research)"}
TIMEOUT = 10
RETRY_MAX = 3

# Keywords that identify non-experimental (computed) entries — reject these.
COMPUTED_KEYWORDS = [
    "Calculated", "Computed", "Estimated", "Predicted",
    "Semi-empirical", "G2", "G3", "G4", "CBS", "PM7", "DFT",
    "B3LYP", "MP2", "HF/",
]


def _is_computed(method_text: str) -> bool:
    return any(kw.lower() in method_text.lower() for kw in COMPUTED_KEYWORDS)


def _parse_nist_html(html: str) -> dict[str, Any] | None:
    """Extract first experimental ΔfH° at 298 K from NIST WebBook HTML."""
    soup = BeautifulSoup(html, "lxml")

    # NIST tables contain 'Delta' or 'fH' in the header text
    for table in soup.find_all("table"):
        headers = [th.get_text(strip=True) for th in table.find_all("th")]
        headers_lower = [h.lower() for h in headers]

        # Check this is a thermochemistry table with enthalpy of formation
        if not any("fh" in h or "delta" in h or "enthalpy" in h for h in headers_lower):
            continue

        # Detect unit column index
        unit = "kj/mol"  # default
        for h in headers_lower:
            if "kcal" in h:
                unit = "kcal/mol"
                break

        # Scan rows for Method column
        method_col = next(
            (i for i, h in enumerate(headers_lower) if "method" in h), None
        )
        ref_col = next(
            (i for i, h in enumerate(headers_lower) if "ref" in h), None
        )

        for row in table.find_all("tr"):
            cells = row.find_all("td")
            if not cells:
                continue

            raw_value = cells[0].get_text(strip=True)
            # Detect unit override from this cell's neighbor
            for cell in cells[1:3]:
                txt = cell.get_text(strip=True).lower()
                if "kcal" in txt:
                    unit = "kcal/mol"
                elif "kj" in txt:
                    unit = "kj/mol"

            method_text = ""
            if method_col is not None and method_col < len(cells):
                method_text = cells[method_col].get_text(strip=True)

            if _is_computed(method_text):
                logger.debug("Skipping computed entry: method=%s", method_text)
                continue

            try:
                value = float(raw_value.replace(",", ""))
            except ValueError:
                continue

            # Convert to kcal/mol
            if unit == "kj/mol":
                value = value / 4.184
            elif "j/mol" in unit and "k" not in unit:
                value = value / 4184.0

            notes = ""
            if ref_col is not None and ref_col < len(cells):
                notes = cells[ref_col].get_text(strip=True)

            return {
                "hf_gas_kcal_mol": value,
                "uncertainty": None,
                "source": "NIST_WebBook",
                "notes": notes,
            }
    return None


def fetch_nist(cas: str) -> dict[str, Any] | None:
    """Fetch ΔfH° from NIST WebBook for a given CAS number."""
    params = {"ID": cas, "Mask": "1", "Type": "JANAFG", "Table": "on"}

    for attempt in range(RETRY_MAX):
        try:
            response = requests.get(
                NIST_URL, params=params, headers=HEADERS, timeout=TIMEOUT
            )
            if response.status_code != 200:
                logger.debug("NIST returned %d for CAS %s", response.status_code, cas)
                return None
            time.sleep(1.5)
            return _parse_nist_html(response.text)
        except Exception as exc:
            wait = 2 ** attempt
            logger.debug("NIST attempt %d failed (%s); retrying in %ds", attempt + 1, exc, wait)
            if attempt < RETRY_MAX - 1:
                time.sleep(wait)

    return None
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_nist_api.py -v
```
Expected: All 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add sources/nist_api.py tests/test_nist_api.py
git commit -m "feat: add NIST WebBook HTML scraper (source 1)"
```

---

## Task 4: sources/nist_chempy.py — NistChemPy wrapper (Source 2)

**Files:**
- Create: `experimental_harvester/sources/nist_chempy.py`
- No dedicated test file — covered by pipeline tests via mocking.

- [ ] **Step 1: Implement sources/nist_chempy.py**

```python
"""sources/nist_chempy.py — NistChemPy wrapper for ΔfH° (Source 2, fallback)."""
import logging
import time
from typing import Any

logger = logging.getLogger(__name__)

COMPUTED_KEYWORDS = [
    "calculated", "computed", "estimated", "predicted",
    "semi-empirical", "g2", "g3", "g4", "cbs", "pm7", "dft",
]


def fetch_nistchempy(inchikey: str) -> dict[str, Any] | None:
    """Fetch ΔfH° via nistchempy library using InChIKey."""
    try:
        import nistchempy as nc  # type: ignore[import]

        compound = nc.get_compound(inchikey)
        if compound is None:
            logger.debug("NistChemPy: no compound for InChIKey %s", inchikey)
            return None

        thermo = getattr(compound, "thermo_gas", None)
        if thermo is None:
            logger.debug("NistChemPy: no gas-phase thermo for %s", inchikey)
            return None

        hf298 = getattr(thermo, "Hf298", None)
        if hf298 is None:
            return None

        # Check method for computed flags
        method = str(getattr(thermo, "method", "")).lower()
        if any(kw in method for kw in COMPUTED_KEYWORDS):
            logger.debug("NistChemPy: rejecting computed value (method=%s)", method)
            return None

        # nistchempy returns values in kJ/mol
        hf_kcal = float(hf298) / 4.184

        time.sleep(1.0)
        return {
            "hf_gas_kcal_mol": hf_kcal,
            "uncertainty": None,
            "source": "NistChemPy",
            "notes": f"method={method}",
        }

    except Exception as exc:
        logger.debug("NistChemPy error for %s: %s", inchikey, exc)
        return None
```

- [ ] **Step 2: Commit**

```bash
git add sources/nist_chempy.py
git commit -m "feat: add NistChemPy wrapper (source 2, fallback)"
```

---

## Task 5: sources/pubchem.py — PubChem REST (Source 3)

**Files:**
- Create: `experimental_harvester/sources/pubchem.py`
- Create: `experimental_harvester/tests/test_pubchem.py`

### 5a — Write failing tests first

- [ ] **Step 1: Write failing tests**

`tests/test_pubchem.py`:
```python
"""Tests for sources/pubchem.py — PubChem PUG REST."""
import unittest
from unittest.mock import MagicMock, patch

from sources import pubchem

CID_JSON = {
    "IdentifierList": {"CID": [962]}
}

PUBCHEM_VIEW_EXPERIMENTAL = {
    "Record": {
        "Section": [
            {
                "TOCHeading": "Chemical and Physical Properties",
                "Section": [
                    {
                        "TOCHeading": "Experimental Properties",
                        "Section": [
                            {
                                "TOCHeading": "Enthalpy of Formation",
                                "Information": [
                                    {
                                        "Value": {
                                            "StringWithMarkup": [
                                                {"String": "-57.8 kcal/mol"}
                                            ]
                                        },
                                        "Reference": ["NIST 2023"],
                                    }
                                ],
                            }
                        ],
                    }
                ],
            }
        ]
    }
}

PUBCHEM_VIEW_COMPUTED_ONLY = {
    "Record": {
        "Section": [
            {
                "TOCHeading": "Chemical and Physical Properties",
                "Section": [
                    {
                        "TOCHeading": "Computed Properties",
                        "Section": [
                            {
                                "TOCHeading": "Enthalpy of Formation",
                                "Information": [
                                    {
                                        "Value": {"StringWithMarkup": [{"String": "-240 kJ/mol"}]},
                                        "Reference": ["DFT 2020"],
                                    }
                                ],
                            }
                        ],
                    }
                ],
            }
        ]
    }
}


class TestFetchPubchem(unittest.TestCase):
    def _mock(self, json_data: dict, status: int = 200) -> MagicMock:
        m = MagicMock()
        m.status_code = status
        m.json.return_value = json_data
        return m

    def test_returns_experimental_value(self):
        with patch(
            "requests.get",
            side_effect=[self._mock(CID_JSON), self._mock(PUBCHEM_VIEW_EXPERIMENTAL)],
        ):
            result = pubchem.fetch_pubchem("O")
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result["hf_gas_kcal_mol"], -57.8, places=1)
        self.assertEqual(result["source"], "PubChem")

    def test_ignores_computed_section(self):
        with patch(
            "requests.get",
            side_effect=[self._mock(CID_JSON), self._mock(PUBCHEM_VIEW_COMPUTED_ONLY)],
        ):
            result = pubchem.fetch_pubchem("O")
        self.assertIsNone(result)

    def test_smiles_not_found_returns_none(self):
        with patch("requests.get", return_value=self._mock({}, 404)):
            result = pubchem.fetch_pubchem("NOT_VALID")
        self.assertIsNone(result)

    def test_kj_conversion(self):
        view_kj = {
            "Record": {
                "Section": [
                    {
                        "TOCHeading": "Chemical and Physical Properties",
                        "Section": [
                            {
                                "TOCHeading": "Experimental Properties",
                                "Section": [
                                    {
                                        "TOCHeading": "Enthalpy of Formation",
                                        "Information": [
                                            {
                                                "Value": {
                                                    "StringWithMarkup": [
                                                        {"String": "-241.826 kJ/mol"}
                                                    ]
                                                },
                                                "Reference": ["CODATA"],
                                            }
                                        ],
                                    }
                                ],
                            }
                        ],
                    }
                ]
            }
        }
        with patch(
            "requests.get",
            side_effect=[self._mock(CID_JSON), self._mock(view_kj)],
        ):
            result = pubchem.fetch_pubchem("O")
        self.assertIsNotNone(result)
        self.assertAlmostEqual(result["hf_gas_kcal_mol"], -241.826 / 4.184, places=2)
```

- [ ] **Step 2: Run to confirm failure**

```bash
python -m pytest tests/test_pubchem.py -v
```
Expected: ImportError.

### 5b — Implement pubchem.py

- [ ] **Step 3: Implement sources/pubchem.py**

```python
"""sources/pubchem.py — PubChem PUG REST for ΔfH° (Source 3, last resort)."""
import logging
import time
import urllib.parse
from typing import Any

import requests

logger = logging.getLogger(__name__)

PUBCHEM_PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"
TIMEOUT = 15
RETRY_MAX = 2

# Only sections with these headings contain experimental (not computed) data.
COMPUTED_HEADINGS = {"computed properties", "computed", "predicted properties"}


def _parse_value(text: str) -> tuple[float, str] | None:
    """Parse '−57.8 kcal/mol' → (value, unit). Returns None if unparseable."""
    text = text.strip().replace("−", "-").replace("\u2212", "-")
    unit = "kcal/mol"
    lower = text.lower()
    if "kj/mol" in lower:
        unit = "kj/mol"
        text = lower.replace("kj/mol", "").strip()
    elif "j/mol" in lower:
        unit = "j/mol"
        text = lower.replace("j/mol", "").strip()
    else:
        text = lower.replace("kcal/mol", "").strip()
    try:
        return float(text), unit
    except ValueError:
        return None


def _find_experimental_hf(record: dict[str, Any]) -> tuple[float, str, str] | None:
    """Recursively search record for experimental ΔfH° value.

    Returns (raw_value, unit, notes) or None.
    Skips entire 'Computed' / 'Predicted' section branches.
    """
    heading = record.get("TOCHeading", "").lower()

    # Skip entire computed branches without recursing into them
    if heading in COMPUTED_HEADINGS:
        return None

    in_hf_section = "enthalpy" in heading and "formation" in heading

    if in_hf_section:
        for info in record.get("Information", []):
            val_obj = info.get("Value", {})
            for markup in val_obj.get("StringWithMarkup", []):
                parsed = _parse_value(markup.get("String", ""))
                if parsed is not None:
                    raw_value, unit = parsed
                    ref = ""
                    refs = info.get("Reference", [])
                    if refs:
                        ref = refs[0] if isinstance(refs[0], str) else str(refs[0])
                    return raw_value, unit, ref
        return None

    for sub in record.get("Section", []):
        result = _find_experimental_hf(sub)
        if result:
            return result

    return None


def fetch_pubchem(smiles: str) -> dict[str, Any] | None:
    """Fetch ΔfH° from PubChem using SMILES."""
    for attempt in range(RETRY_MAX + 1):
        try:
            # Step 1: SMILES → CID
            encoded = urllib.parse.quote(smiles, safe="")
            r1 = requests.get(
                f"{PUBCHEM_PUG}/compound/smiles/{encoded}/cids/JSON",
                timeout=TIMEOUT,
            )
            if r1.status_code != 200:
                logger.debug("PubChem CID not found for SMILES: %s", smiles)
                return None

            cids = r1.json().get("IdentifierList", {}).get("CID", [])
            if not cids:
                return None
            cid = cids[0]
            time.sleep(1.0)

            # Step 2: CID → thermochemical data
            r2 = requests.get(
                f"{PUBCHEM_VIEW}/data/compound/{cid}/JSON/?heading=Enthalpy+of+Formation",
                timeout=TIMEOUT,
            )
            if r2.status_code != 200:
                return None

            record = r2.json().get("Record", {})
            found = _find_experimental_hf(record)
            if found is None:
                return None

            raw_value, unit, notes = found
            if unit == "kj/mol":
                hf_kcal = raw_value / 4.184
            elif unit == "j/mol":
                hf_kcal = raw_value / 4184.0
            else:
                hf_kcal = raw_value

            return {
                "hf_gas_kcal_mol": hf_kcal,
                "uncertainty_kcal_mol": None,
                "source": "PubChem",
                "notes": notes,
            }

        except Exception as exc:
            logger.debug("PubChem attempt %d failed: %s", attempt + 1, exc)
            if attempt < RETRY_MAX:
                time.sleep(1.0)

    return None
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_pubchem.py -v
```
Expected: All 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add sources/pubchem.py tests/test_pubchem.py
git commit -m "feat: add PubChem REST source (source 3, last resort)"
```

---

## Task 6: pipeline.py — Cascading orchestrator with checkpoint

**Files:**
- Create: `experimental_harvester/pipeline.py`
- Create: `experimental_harvester/tests/test_pipeline.py`

### 6a — Failing tests first

- [ ] **Step 1: Write failing tests**

`tests/test_pipeline.py`:
```python
"""Tests for pipeline.py — cascade orchestration and batch processing."""
import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path
from unittest.mock import MagicMock, patch

import pipeline
from pipeline import HarvestResult


class TestHarvestOne(unittest.TestCase):
    def _nist_result(self) -> dict:
        return {
            "hf_gas_kcal_mol": -57.8,
            "uncertainty": 0.1,
            "source": "NIST_WebBook",
            "notes": "CODATA",
        }

    def _pubchem_result(self) -> dict:
        return {
            "hf_gas_kcal_mol": -57.5,
            "uncertainty_kcal_mol": None,
            "source": "PubChem",
            "notes": "ref",
        }

    def test_uses_nist_first(self):
        with (
            patch("resolver.smiles_to_inchikey", return_value="INCHI-KEY"),
            patch("resolver.smiles_to_cas", return_value="7732-18-5"),
            patch("sources.nist_api.fetch_nist", return_value=self._nist_result()),
        ):
            result = pipeline.harvest_one("O")

        self.assertEqual(result.source, "NIST_WebBook")
        self.assertAlmostEqual(result.hf_gas_kcal_mol, -57.8)

    def test_fallback_to_nistchempy_when_nist_fails(self):
        nc_result = {
            "hf_gas_kcal_mol": -57.7,
            "uncertainty": None,
            "source": "NistChemPy",
            "notes": "",
        }
        with (
            patch("resolver.smiles_to_inchikey", return_value="INCHI-KEY"),
            patch("resolver.smiles_to_cas", return_value="7732-18-5"),
            patch("sources.nist_api.fetch_nist", return_value=None),
            patch("sources.nist_chempy.fetch_nistchempy", return_value=nc_result),
        ):
            result = pipeline.harvest_one("O")
        self.assertEqual(result.source, "NistChemPy")

    def test_fallback_to_pubchem_when_both_nist_fail(self):
        with (
            patch("resolver.smiles_to_inchikey", return_value="INCHI-KEY"),
            patch("resolver.smiles_to_cas", return_value=None),
            patch("sources.nist_api.fetch_nist", return_value=None),
            patch("sources.nist_chempy.fetch_nistchempy", return_value=None),
            patch("sources.pubchem.fetch_pubchem", return_value=self._pubchem_result()),
        ):
            result = pipeline.harvest_one("O")
        self.assertEqual(result.source, "PubChem")

    def test_all_sources_fail_returns_not_found(self):
        with (
            patch("resolver.smiles_to_inchikey", return_value=None),
            patch("resolver.smiles_to_cas", return_value=None),
            patch("sources.nist_api.fetch_nist", return_value=None),
            patch("sources.nist_chempy.fetch_nistchempy", return_value=None),
            patch("sources.pubchem.fetch_pubchem", return_value=None),
        ):
            result = pipeline.harvest_one("UNKNOWN")
        self.assertEqual(result.source, "not_found")
        self.assertIsNone(result.hf_gas_kcal_mol)

    def test_large_value_logs_warning_but_is_kept(self):
        big_result = {
            "hf_gas_kcal_mol": 600.0,
            "uncertainty": None,
            "source": "NIST_WebBook",
            "notes": "",
        }
        with (
            patch("resolver.smiles_to_inchikey", return_value="INCHI-KEY"),
            patch("resolver.smiles_to_cas", return_value="1234-56-7"),
            patch("sources.nist_api.fetch_nist", return_value=big_result),
        ):
            import logging
            with self.assertLogs("pipeline", level="WARNING"):
                result = pipeline.harvest_one("C" * 20)
        # Value kept despite warning
        self.assertIsNotNone(result.hf_gas_kcal_mol)
        self.assertEqual(result.hf_gas_kcal_mol, 600.0)


class TestHarvestBatch(unittest.TestCase):
    def test_skips_already_processed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint = Path(tmpdir) / "checkpoint.csv"
            # Pre-populate checkpoint with one SMILES
            import pandas as pd
            pd.DataFrame([{
                "smiles": "O",
                "hf_gas_kcal_mol": -57.8,
                "source": "NIST_WebBook",
                "uncertainty_kcal_mol": None,
                "notes": "",
            }]).to_csv(checkpoint, index=False)

            call_count = {"n": 0}

            def fake_harvest_one(smi: str) -> HarvestResult:
                call_count["n"] += 1
                return HarvestResult(smiles=smi, hf_gas_kcal_mol=None,
                                     source="not_found", uncertainty_kcal_mol=None, notes="")

            with patch("pipeline.harvest_one", side_effect=fake_harvest_one):
                results = pipeline.harvest_batch(["O", "CC"], checkpoint)

        # "O" was already in checkpoint — harvest_one called only once (for "CC")
        self.assertEqual(call_count["n"], 1)
        self.assertEqual(len(results), 2)

    def test_resumes_from_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint = Path(tmpdir) / "checkpoint.csv"
            import pandas as pd
            pd.DataFrame([{
                "smiles": "O",
                "hf_gas_kcal_mol": -57.8,
                "source": "NIST_WebBook",
                "uncertainty_kcal_mol": None,
                "notes": "",
            }]).to_csv(checkpoint, index=False)

            with patch("pipeline.harvest_one") as mock_ho:
                mock_ho.return_value = HarvestResult(
                    smiles="CC", hf_gas_kcal_mol=-20.0,
                    source="NIST_WebBook", uncertainty_kcal_mol=None, notes=""
                )
                results = pipeline.harvest_batch(["O", "CC"], checkpoint)

            mock_ho.assert_called_once_with("CC")
            sources = {r.smiles: r.source for r in results}
            self.assertEqual(sources["O"], "NIST_WebBook")
```

- [ ] **Step 2: Run to confirm failure**

```bash
python -m pytest tests/test_pipeline.py -v
```
Expected: ImportError.

### 6b — Implement pipeline.py

- [ ] **Step 3: Implement pipeline.py**

```python
"""pipeline.py — Cascade harvesting orchestrator with checkpoint/resume."""
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
from tqdm import tqdm

import resolver
import sources.nist_api as nist_api
import sources.nist_chempy as nist_chempy
import sources.pubchem as pubchem

logger = logging.getLogger(__name__)

CHECKPOINT_INTERVAL = 50
HF_RANGE_WARN = 500.0  # kcal/mol — log warning if |hf| exceeds this

CSV_COLUMNS = ["smiles", "hf_gas_kcal_mol", "source", "uncertainty_kcal_mol", "notes"]


@dataclass
class HarvestResult:
    smiles: str
    hf_gas_kcal_mol: float | None
    source: str  # "NIST_WebBook" | "NistChemPy" | "PubChem" | "not_found"
    uncertainty_kcal_mol: float | None
    notes: str


def _result_from_dict(smiles: str, data: dict[str, Any]) -> HarvestResult:
    return HarvestResult(
        smiles=smiles,
        hf_gas_kcal_mol=data.get("hf_gas_kcal_mol"),
        source=data.get("source", "unknown"),
        uncertainty_kcal_mol=data.get("uncertainty") or data.get("uncertainty_kcal_mol"),
        notes=data.get("notes", ""),
    )


def harvest_one(smiles: str) -> HarvestResult:
    """Cascade through sources to find experimental ΔfH° for a single SMILES."""
    inchikey = resolver.smiles_to_inchikey(smiles)
    cas = resolver.smiles_to_cas(smiles)

    # Source 1: NIST WebBook (requires CAS)
    if cas:
        data = nist_api.fetch_nist(cas)
        if data:
            result = _result_from_dict(smiles, data)
            _check_range(result)
            return result

    # Source 2: NistChemPy (requires InChIKey)
    if inchikey:
        data = nist_chempy.fetch_nistchempy(inchikey)
        if data:
            result = _result_from_dict(smiles, data)
            _check_range(result)
            return result

    # Source 3: PubChem (SMILES directly)
    data = pubchem.fetch_pubchem(smiles)
    if data:
        result = _result_from_dict(smiles, data)
        _check_range(result)
        return result

    return HarvestResult(
        smiles=smiles,
        hf_gas_kcal_mol=None,
        source="not_found",
        uncertainty_kcal_mol=None,
        notes="",
    )


def _check_range(result: HarvestResult) -> None:
    if result.hf_gas_kcal_mol is not None and abs(result.hf_gas_kcal_mol) > HF_RANGE_WARN:
        logger.warning(
            "Large |ΔfH°| = %.1f kcal/mol for SMILES %s (source: %s) — verify manually",
            result.hf_gas_kcal_mol, result.smiles, result.source,
        )


def _load_checkpoint(checkpoint_path: Path) -> dict[str, HarvestResult]:
    """Load previously processed results from checkpoint CSV."""
    if not checkpoint_path.exists():
        return {}
    df = pd.read_csv(checkpoint_path)
    results: dict[str, HarvestResult] = {}
    for _, row in df.iterrows():
        r = HarvestResult(
            smiles=str(row["smiles"]),
            hf_gas_kcal_mol=None if pd.isna(row["hf_gas_kcal_mol"]) else float(row["hf_gas_kcal_mol"]),
            source=str(row["source"]),
            uncertainty_kcal_mol=None if pd.isna(row.get("uncertainty_kcal_mol", float("nan"))) else float(row["uncertainty_kcal_mol"]),
            notes=str(row.get("notes", "")),
        )
        results[r.smiles] = r
    return results


def _save_checkpoint(results: dict[str, HarvestResult], checkpoint_path: Path) -> None:
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "smiles": r.smiles,
            "hf_gas_kcal_mol": r.hf_gas_kcal_mol,
            "source": r.source,
            "uncertainty_kcal_mol": r.uncertainty_kcal_mol,
            "notes": r.notes,
        }
        for r in results.values()
    ]
    pd.DataFrame(rows, columns=CSV_COLUMNS).to_csv(checkpoint_path, index=False)


def harvest_batch(
    smiles_list: list[str], checkpoint_path: Path
) -> list[HarvestResult]:
    """Process a list of SMILES with checkpoint/resume support."""
    done = _load_checkpoint(checkpoint_path)
    pending = [s for s in smiles_list if s not in done]

    logger.info(
        "Batch: %d total, %d already done, %d to process",
        len(smiles_list), len(done), len(pending),
    )

    for i, smiles in enumerate(tqdm(pending, desc="Harvesting")):
        result = harvest_one(smiles)
        done[smiles] = result

        if (i + 1) % CHECKPOINT_INTERVAL == 0:
            _save_checkpoint(done, checkpoint_path)
            logger.info("Checkpoint saved at %d molecules", i + 1)

    # Final checkpoint
    _save_checkpoint(done, checkpoint_path)

    # Return in original smiles_list order
    return [done[s] for s in smiles_list if s in done]
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_pipeline.py -v
```
Expected: All 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add pipeline.py tests/test_pipeline.py
git commit -m "feat: add pipeline — cascade harvesting with checkpoint/resume"
```

---

## Task 7: harvest.py — CLI entry point

**Files:**
- Create: `experimental_harvester/harvest.py`

- [ ] **Step 1: Implement harvest.py**

```python
"""harvest.py — CLI entry point for experimental ΔfH° harvesting."""
import argparse
import logging
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

import pipeline

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(name)s — %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("harvest")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Harvest experimental ΔfH° for SMILES from NIST/PubChem."
    )
    parser.add_argument("--input", required=True, help="Path to input CSV with SMILES")
    parser.add_argument(
        "--output", default="output/experimental.csv", help="Output CSV path"
    )
    parser.add_argument(
        "--smiles-col", default="smiles", help="Name of the SMILES column (default: smiles)"
    )
    parser.add_argument(
        "--resume", action="store_true", help="Resume from existing checkpoint"
    )
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    args = parser.parse_args()

    logging.getLogger().setLevel(args.log_level)

    # Read input
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error("Input file not found: %s", input_path)
        sys.exit(1)

    df = pd.read_csv(input_path)
    if args.smiles_col not in df.columns:
        logger.error(
            "Column '%s' not found. Available: %s", args.smiles_col, list(df.columns)
        )
        sys.exit(1)

    smiles_list = df[args.smiles_col].dropna().unique().tolist()
    logger.info("Found %d unique SMILES. Starting harvest...", len(smiles_list))

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output_path.parent / "checkpoint.csv"

    if not args.resume and checkpoint_path.exists():
        checkpoint_path.unlink()
        logger.info("Cleared existing checkpoint (use --resume to continue)")

    results = pipeline.harvest_batch(smiles_list, checkpoint_path)

    # Write output
    rows = [
        {
            "smiles": r.smiles,
            "hf_gas_kcal_mol": r.hf_gas_kcal_mol,
            "source": r.source,
            "uncertainty_kcal_mol": r.uncertainty_kcal_mol,
            "notes": r.notes,
        }
        for r in results
    ]
    out_df = pd.DataFrame(rows)
    out_df.to_csv(output_path, index=False)
    logger.info("Output written to %s", output_path)

    # Summary
    total = len(results)
    found = sum(1 for r in results if r.hf_gas_kcal_mol is not None)
    not_found = total - found
    source_counts = Counter(r.source for r in results if r.hf_gas_kcal_mol is not None)

    print(
        f"\nCompleted: {total} total | found: {found} | not_found: {not_found} | "
        f"sources: NIST={source_counts.get('NIST_WebBook', 0)}, "
        f"NistChemPy={source_counts.get('NistChemPy', 0)}, "
        f"PubChem={source_counts.get('PubChem', 0)}"
    )

    # Write harvest log
    log_path = output_path.parent / "harvest.log"
    with log_path.open("w") as f:
        for r in results:
            f.write(
                f"{r.smiles}\t{r.source}\t{r.hf_gas_kcal_mol}\t{r.notes}\n"
            )
    logger.info("Harvest log written to %s", log_path)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke test (no network)**

```bash
# Quick structural check — will fail fast without actual input file
python harvest.py --help
```
Expected: Prints usage without error.

- [ ] **Step 3: Commit**

```bash
git add harvest.py
git commit -m "feat: add harvest.py CLI entry point"
```

---

## Task 8: Full test suite and quality gates

- [ ] **Step 1: Run all tests**

```bash
cd experimental_harvester
python -m pytest tests/ -v --tb=short
```
Expected: All tests PASS.

- [ ] **Step 2: Run ruff**

```bash
ruff check .
```
Fix any reported issues.

- [ ] **Step 3: Run black**

```bash
black --check .
```
Fix formatting if needed: `black .`

- [ ] **Step 4: Run mypy**

```bash
mypy resolver.py pipeline.py harvest.py sources/ --strict
```
Fix any type errors reported.

- [ ] **Step 5: Final commit**

```bash
git add -A
git commit -m "chore: pass all quality gates (ruff, black, mypy, pytest)"
```

---

## Data Quality Reference

These rules are enforced at parse time in each source module:

| Rule | Where enforced |
|------|---------------|
| Reject Calculated/DFT/CBS/etc. methods | `nist_api._is_computed`, `nist_chempy`, `pubchem._find_experimental_hf` |
| kJ/mol → kcal/mol: ÷ 4.184 | All three source modules |
| J/mol → kcal/mol: ÷ 4184 | All three source modules |
| Gas phase only | NIST filter on table type (JANAFG = gas) |
| Reference temperature 298 K | NIST Mask=1 + JANAFG type parameter |
| Log reference in notes | All sources populate `notes` field |
| Warn if \|hf\| > 500 kcal/mol | `pipeline._check_range` |

---

## Input/Output Contract

**Input:** `data/thermo_pm7.csv` — confirmed column name: `smiles` (lowercase)

**Output:** `output/experimental.csv`

```
smiles,hf_gas_kcal_mol,source,uncertainty_kcal_mol,notes
O,-57.8,NIST_WebBook,,CODATA 1989
CC,-20.0,PubChem,,NIST 2023
C1CCCCC1,,not_found,,
```
