---
name: grimperium-tests
description: Estratégia de execução de testes para o layout atual do Grimperium
tools: [bash]
context: fork
user-invocable: true
allowed-tools:
  - Bash(poetry *)
  - Bash(pytest *)
---

# Skill: Grimperium Test Execution

Use this skill for deciding which tests to run and in what order.

## Current Test Layout

- `tests/unit/`
- `tests/integration/`
- `tests/experiments/`
- `tests/cli/`
- `tests/ml/`

`pyproject.toml` defines the active pytest configuration and markers. Trust that
over stale prose counts.

## Recommended Order

1. Run the closest test file to the change.
2. Run the closest package or marker slice if the behavior crosses files.
3. Run `pytest tests/ -v` or `pytest tests/ -v --cov=src/grimperium` before
   handoff when the change is substantial.

## Useful Commands

Single test:

```bash
pytest tests/unit/test_file.py::test_name -v
```

Pattern match:

```bash
pytest -k "pattern" -v
```

Directory slice:

```bash
pytest tests/unit/ -v
pytest tests/integration/ -v
pytest tests/ml/ -v
```

Full suite with coverage:

```bash
pytest tests/ -v --cov=src/grimperium
```

## Notes

- Prefer targeted reruns during iteration.
- If a failure crosses CLI, data, and pipeline code, expand validation beyond a
  single module.
- Report the exact command run instead of saying "tests passed" generically.
