# Grimperium Testing Guide
## TDD-First Culture + Phase C Quality Standards

**Last Updated:** 2026-01-17
**Status:** ✅ 145 tests passing, 82% coverage (target 85%)

---

## Testing Philosophy

1. **Write tests FIRST** (TDD red-green-refactor)
2. **Coverage ≥85%** (mutation testing recommended later)
3. **Clear test names** (test_[function]_[scenario]_[expected])
4. **Fixtures for reuse** (CSV data, mocked PM7 results)
5. **Edge cases always** (None, empty, invalid inputs)

---

## Current Test Suite

### Status (2026-01-17)
- Total tests: 145
- Passing: 145 (100%)
- Coverage: 82% (target 85%+)
- Tools: pytest, pytest-cov, pytest-mock

### Test Categories
| Category | Count | Coverage |
|----------|-------|----------|
| CLI views | ~40 | 85% |
| Data loading | ~25 | 80% |
| Calculations | ~30 | 88% |
| Utilities | ~20 | 75% |
| Integration | ~30 | 78% |

---

## Test Structure

### Single Test Example
```python
# tests/unit/test_databases_view.py

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

from grimperium.cli.views.databases_view import DatabasesView


def test_load_real_phase_a_results_from_json(tmp_path: Path) -> None:
    """
    Bug #2: Test loading molecule count from phase_a_results.json.

    When JSON exists with n_molecules, that value should be used.
    """
    json_file = tmp_path / "phase_a_results.json"
    json_data = {
        "n_molecules": 29568,
        "results": [
            {"smiles": "CCO", "H298_pm7": -100.5, "timestamp": "2026-01-15T10:30:00"}
        ],
    }
    json_file.write_text(json.dumps(json_data), encoding="utf-8")

    with patch("grimperium.cli.views.databases_view.PHASE_A_RESULTS_FILE", json_file):
        result = DatabasesView.load_real_phase_a_results()

    assert result is not None
    assert result.molecules == 29568
    assert result.name == "CREST PM7"
    assert result.status == "ready"  # Because molecules > 0
```

---

## Running Tests

### Full Suite
```bash
pytest tests/ -v
→ 145 passed in X.XXs
```

### With Coverage
```bash
pytest tests/ --cov=src/ --cov-report=html
→ Coverage report in htmlcov/index.html
```

### Specific Module
```bash
pytest tests/unit/test_databases_view.py -v
→ Tests for database_view only
```

### TDD Workflow
```bash
1. Write failing test
   pytest tests/test_new_feature.py -v
   → RED (expected to fail)

2. Implement minimum code to pass
   [edit src/...]

3. Run test again
   pytest tests/test_new_feature.py -v
   → GREEN (passing)

4. Refactor (maintain tests passing)
   [improve src/...]

5. All tests pass
   pytest tests/ -v
   → 146 passed (new test + existing)
```

---

## Fixtures (Reusable Test Data)

### conftest.py Pattern
```python
# tests/conftest.py

import pytest
from src.data.loader import load_molecules

@pytest.fixture
def sample_molecules():
    """Sample 3-molecule dataset for testing"""
    return {
        'mol_001': {'smiles': 'CCO', 'pm7_energy': -156.234, ...},
        'mol_002': {'smiles': 'CC(C)O', 'pm7_energy': -198.456, ...},
        'mol_003': {'smiles': 'C1=CC=CC=C1', 'pm7_energy': -412.789, ...},
    }

@pytest.fixture
def mock_csv_file(tmp_path):
    """Create temporary CSV for testing data loading"""
    csv_path = tmp_path / "test_molecules.csv"
    csv_path.write_text("mol_id,smiles,pm7_energy\n...")
    return csv_path
```

---

## BATCH 12: Test Coverage Goals

### Current (82%) → Target (85%+)
- Focus on CLI views (database_view, settings_view)
- Add edge case tests for calculations
- Integration tests for full workflows

### New Tests Expected
- ~15-20 new tests
- Each bug fix includes test coverage
- Coverage increase: +3-5% points

---

## Continuous Integration (Future)

### GitHub Actions Template
```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: actions/setup-python@v2
      - run: pip install -r requirements.txt
      - run: pytest tests/ --cov=src/
      - run: mypy --strict src/
      - run: ruff check src/
```

---

**Next:** Add 15-20 tests during BATCH 12 to reach 85% coverage target.
