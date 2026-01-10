# Grimperium - Code Style & Conventions

## Python Style

### Formatting
- **Line Length**: 88 characters (Black default)
- **Formatter**: Black
- **Import Sorting**: Ruff with isort rules
- **Target Python**: 3.9+

### Type Hints
- **Required**: All public functions must have type hints
- **mypy**: Strict mode enabled (`disallow_untyped_defs = true`)
- **Style**: Use `Optional`, `Union` from typing module
- **Arrays**: Use `np.ndarray` for NumPy arrays

```python
def fit(
    self, X: np.ndarray, y_cbs: np.ndarray, y_pm7: np.ndarray
) -> "DeltaLearner":
```

### Docstrings
- **Format**: NumPy-style docstrings
- **Required**: All public classes and methods
- **Sections**: Description, Parameters, Returns, Example

```python
def predict(self, X: np.ndarray, y_pm7: np.ndarray) -> np.ndarray:
    """
    Predict H298_CBS using delta-learning.

    Parameters:
        X: Features (n_samples, n_features)
        y_pm7: PM7 approximation (n_samples,)

    Returns:
        Predicted H298_CBS (n_samples,)

    Example:
        y_pred = learner.predict(X_test, y_pm7_test)
    """
```

## Naming Conventions

### Variables
- **snake_case**: For variables and functions
- **UPPER_CASE**: For constants
- **_private**: Single underscore for internal use
- **Prefixes**: `y_` for target arrays, `X` for feature matrices

### Classes
- **PascalCase**: Class names
- **Base prefix**: `BaseModel` for abstract classes
- **Descriptive**: `DeltaLearner`, `KernelRidgeModel`

### Files
- **snake_case**: All Python files
- **test_ prefix**: Test files (e.g., `test_delta_learning.py`)

## Code Patterns

### Class Structure
```python
class ClassName:
    """Class docstring."""
    
    def __init__(self, ...):
        """Initialize with parameters."""
        self.attr: Type = value
        self._private: bool = False
        
    def __repr__(self) -> str:
        return f"ClassName(...)"
    
    def public_method(self) -> ReturnType:
        """Method docstring."""
        pass
```

### Error Handling
- Use explicit checks with descriptive messages
- Raise appropriate exceptions
- Validate inputs at method boundaries

```python
if not self.is_fitted:
    raise ValueError("DeltaLearner not fitted. Call fit() first.")
```

### Assertions
- Use `assert` for development-time sanity checks
- Include descriptive messages

```python
assert X.shape[0] == len(y_cbs), "X and y_cbs must have same n_samples"
```

## Ruff Rules (Enabled)
- E, W: pycodestyle errors/warnings
- F: Pyflakes
- I: isort
- B: flake8-bugbear
- C4: flake8-comprehensions
- UP: pyupgrade
- ARG: unused arguments
- SIM: simplify

## Testing Conventions
- Use pytest
- Test files: `tests/unit/test_*.py`
- Test functions: `def test_*():`
- Use fixtures from `conftest.py`
- Mark slow tests with `@pytest.mark.slow`
