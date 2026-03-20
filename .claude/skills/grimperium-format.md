---
name: grimperium-format
description: Workflow de formatação e lint alinhado ao pyproject atual
tools: [bash]
context: fork
user-invocable: true
allowed-tools:
  - Bash(black *)
  - Bash(ruff *)
---

# Skill: Grimperium Format And Lint

Use this skill when the task is formatting or lint cleanup.

## Commands

Format:

```bash
black src/ tests/
```

Lint:

```bash
ruff check src/ tests/
```

## Usage Notes

- `black` is the formatter of record for Python files.
- `ruff` is the lint gate; keep it clean before handoff.
- If you only want to inspect formatting drift, use `black --check src/ tests/`.
- Do not claim full verification from a formatting pass alone; pair it with the
  relevant type or test command when behavior changed.
