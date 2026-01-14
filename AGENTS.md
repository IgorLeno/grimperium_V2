# AGENTS.md - Grimperium

ML framework for molecular thermodynamic property prediction with delta-learning. Python 3.10+.

## Commands
- **Test all:** `pytest tests/ -v --cov=src/grimperium`
- **Single test:** `pytest tests/unit/test_file.py::test_name -v`
- **Lint:** `ruff check src/ tests/` (zero warnings required)
- **Format:** `black src/ tests/` (line-length 88)
- **Type check:** `mypy src/ --strict`
- **Pre-commit:** `pre-commit run --all-files`

## Architecture
- `src/grimperium/crest_pm7/` - Phase A: conformer generation, threshold monitoring, timeout prediction
- `src/grimperium/core/` - Delta learning + metrics
- `src/grimperium/models/` - ML models (XGBoost, KernelRidge, ensemble)
- `src/grimperium/data/` - Data loading and fusion
- `tests/{unit,integration}/` - Tests with fixtures in `conftest.py`

## Code Style
- **Type hints:** 100% required on all functions
- **Docstrings:** Google-style with Args/Returns/Example sections
- **Naming:** `snake_case` (functions), `PascalCase` (classes), `UPPER_CASE` (constants)
- **Coverage:** Minimum 85%, markers: `@pytest.mark.unit`, `@pytest.mark.integration`
- **No TODOs, placeholders, or skipped tests**

See `docs/CLAUDE.md` for detailed project guidelines.
