# Grimperium - Suggested Commands

## Development Environment

### Install Dependencies
```bash
poetry install                    # Install all dependencies
poetry install --with dev         # Include dev dependencies
poetry install --extras docs      # Include documentation deps
```

### Activate Virtual Environment
```bash
poetry shell                      # Activate poetry venv
# or
source $(poetry env info --path)/bin/activate
```

## Testing

### Run All Tests
```bash
pytest                            # Run all tests
pytest -v                         # Verbose output
pytest --tb=short                 # Short traceback
```

### Run Specific Test Files
```bash
pytest tests/unit/                # Unit tests only
pytest tests/integration/         # Integration tests only
pytest tests/unit/test_delta_learning.py  # Single file
```

### Run with Coverage
```bash
pytest --cov=src/grimperium --cov-report=term-missing
pytest --cov=src/grimperium --cov-report=html  # HTML report
```

### Run Parallel Tests
```bash
pytest -n auto                    # Use all CPUs
pytest -n 4                       # Use 4 workers
```

### Skip Slow Tests
```bash
pytest -m "not slow"              # Skip slow tests
pytest -m "not integration"       # Skip integration tests
```

## Linting & Formatting

### Ruff (Linting)
```bash
ruff check .                      # Check for issues
ruff check . --fix                # Auto-fix issues
```

### Black (Formatting)
```bash
black .                           # Format all files
black --check .                   # Check without modifying
black --diff .                    # Show what would change
```

### MyPy (Type Checking)
```bash
mypy src/                         # Check source files
mypy src/ --ignore-missing-imports
```

### All Pre-commit Hooks
```bash
pre-commit run --all-files        # Run all hooks
pre-commit install                # Install git hooks
```

## Build & Package

### Build Package
```bash
poetry build                      # Build wheel and sdist
```

### Documentation
```bash
pdoc --html src/grimperium -o docs/pdoc
sphinx-build docs/source docs/build/html
```

## Tox (Multi-Python Testing)
```bash
tox                               # Run all environments
tox -e py310                      # Python 3.10 only
tox -e lint                       # Lint only
```

## Git Workflow
```bash
git status                        # Check status
git diff                          # View changes
git add -p                        # Interactive staging
git commit -m "type: description" # Commit
```

## Utility Commands
```bash
ls -la                            # List files
grep -r "pattern" src/            # Search in source
find . -name "*.py" -type f       # Find Python files
cat pyproject.toml                # View config
```
