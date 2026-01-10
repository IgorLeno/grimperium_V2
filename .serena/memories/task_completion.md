# Grimperium - Task Completion Checklist

## Before Completing Any Task

### 1. Run Linting

```bash
ruff check . --fix
black .
```

Ensure no linting errors remain.

### 2. Run Type Checking

```bash
mypy src/ --ignore-missing-imports
```

Fix any type errors before proceeding.

### 3. Run Tests

```bash
pytest -v
```

All tests must pass. If modifying core functionality, run:

```bash
pytest --cov=src/grimperium --cov-report=term-missing
```


### 4. Pre-commit Hooks

```bash
pre-commit run --all-files
```

This runs all quality checks automatically.

## For Code Changes


### Added New Functionality

- [ ] Added appropriate type hints
- [ ] Added docstring with Parameters/Returns/Example
- [ ] Added unit tests in `tests/unit/`
- [ ] Updated relevant README.md

### Modified Existing Code
- [ ] Existing tests still pass
- [ ] No regression in functionality
- [ ] Updated docstrings if API changed

### Bug Fixes
- [ ] Added test case that reproduces the bug
- [ ] Verified fix resolves the issue
- [ ] No new regressions introduced

## Commit Guidelines


### Commit Message Format

```
type: brief description

- Detailed bullet points if needed
- Reference issues: Fixes #123
```

### Types
- `feat`: New feature
- `fix`: Bug fix
- `refactor`: Code restructuring
- `test`: Adding/updating tests
- `docs`: Documentation changes
- `chore`: Maintenance tasks

## Final Verification Sequence

```bash
# Full verification before committing
ruff check . --fix && black . && mypy src/ && pytest -v

# Or use pre-commit
pre-commit run --all-files
```
