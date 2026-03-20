---
name: grimperium-ci-fix
description: Reproduz e triageia falhas de CI usando os comandos nativos do repositório
tools: [bash, file]
context: fork
user-invocable: true
allowed-tools:
  - Bash(poetry *)
  - Bash(black *)
  - Bash(ruff *)
  - Bash(mypy *)
  - Bash(git *)
---

# Skill: Grimperium CI Triage

Use this skill when a CI report, failing workflow log, or reproduced local gate
needs to be turned into a concrete debugging loop.

`CLAUDE.md` remains the main authority. This skill is only a focused workflow.

## Inputs

- A CI error summary or workflow log, or
- A local failing command that matches the CI symptom

## Reproduction Order

1. Identify the failing job and exact command from the log.
2. Reproduce the smallest failing check locally.
3. Fix the smallest real cause, not the broadest symptom.
4. Rerun the same failing check.
5. Run the broader related gate before handoff.

## Repo-Native Commands

Formatting and lint:

```bash
black src/ tests/
ruff check src/ tests/
```

Type checking:

```bash
mypy src/ --strict
```

Tests:

```bash
pytest tests/ -v --tb=short
pytest path/to/test_file.py::test_name -v
pytest -k "pattern" -v
```

## Triage Rules

- Do not auto-commit or auto-push.
- Do not claim CI is fixed until the relevant local gate reruns cleanly.
- Prefer targeted reruns while iterating, then a broader validating rerun.
- For dataset, schema, or parser failures, confirm the current repo state before
  trusting the log's narrative.
