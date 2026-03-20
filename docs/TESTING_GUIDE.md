# Grimperium Testing Guide

This is a supporting testing reference. The main local behavior contract is
[`/CLAUDE.md`](/home/igor/Projetos/grimperium/CLAUDE.md).

## Current Test Layout

- `tests/unit/`
- `tests/integration/`
- `tests/experiments/`
- `tests/cli/`
- `tests/ml/`
- shared fixtures in `tests/conftest.py` and local `conftest.py` files

`pyproject.toml` is the source of truth for pytest options and markers.

## Default Testing Loop

1. Run the closest relevant test first.
2. Expand to the nearest directory or marker slice if the change spans files.
3. Run the full suite before handoff when the change is substantial.

## Useful Commands

```bash
pytest tests/unit/test_file.py::test_name -v
pytest -k "pattern" -v
pytest tests/unit/ -v
pytest tests/integration/ -v
pytest tests/ -v --cov=src/grimperium
pytest --collect-only -q
```

## Notes

- Use `--collect-only` when you need a current test inventory instead of relying
  on frozen prose counts.
- For CLI, parser, schema, and pipeline work, do not stop at one isolated unit
  test if downstream behavior could change.
- When reporting verification, name the exact test command you ran.
