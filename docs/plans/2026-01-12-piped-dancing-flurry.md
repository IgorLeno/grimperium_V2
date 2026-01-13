# Grimperium CLI Implementation Plan

## Overview
Build a professional, interactive CLI application for GRIMPERIUM using Rich library with 6 main modules: CALC, DATABASES, MODELS, RESULTS, SETTINGS, and ABOUT.

---

## Phase 1: Foundation & Architecture

### 1.1 File Structure
```
src/grimperium/cli/
├── __init__.py           # CLI package exports
├── __main__.py           # Entry point: python -m grimperium.cli
├── app.py                # Main application class
├── controller.py         # Navigation and state management
├── menu.py               # Menu rendering and selection
├── styles.py             # Colors, themes, ASCII art
├── mock_data.py          # Mock data for MVP
└── views/
    ├── __init__.py
    ├── base_view.py      # Abstract base view class
    ├── calc_view.py      # CALC module (predictions)
    ├── databases_view.py # DATABASES module
    ├── models_view.py    # MODELS module
    ├── results_view.py   # RESULTS module
    ├── settings_view.py  # SETTINGS module [IN DEVELOPMENT]
    └── about_view.py     # ABOUT module
```

### 1.2 Dependencies
```toml
# Add to pyproject.toml [tool.poetry.dependencies]
rich = ">=13.0"
questionary = ">=2.0.0"  # Cross-platform arrow key menu navigation
```

### 1.3 Entry Point Update
```toml
# pyproject.toml
[tool.poetry.scripts]
grimperium = "grimperium.api:main"
grimperium-cli = "grimperium.cli:main"  # New CLI entry
```

---

## Phase 2: Core Components

### 2.1 styles.py - Theme & Colors
```python
# Color scheme
COLORS = {
    "calc": "#00D9FF",        # Cyan
    "databases": "#0080FF",   # Blue
    "models": "#FF00FF",      # Magenta
    "results": "#00FF00",     # Green
    "settings": "#FFFF00",    # Yellow
    "about": "#F5F5F5",       # White
    "error": "#FF4444",       # Red
    "success": "#44FF44",     # Green
    "warning": "#FFAA00",     # Orange
    "in_dev": "#888888",      # Gray
}

# ASCII art banner for GRIMPERIUM
BANNER = """..."""
```

### 2.2 controller.py - State Management
```python
class CliController:
    """Manages navigation state and view transitions."""

    def __init__(self):
        self.history: list[str] = []  # Breadcrumb navigation
        self.current_view: str = "main"
        self.current_model: str = "DeltaRF_v1.0"
        self.console = Console()

    def navigate_to(self, view: str) -> None
    def go_back(self) -> bool
    def run(self) -> int  # Main loop, returns exit code
```

### 2.3 menu.py - Interactive Menu
```python
def create_menu(
    options: list[MenuOption],
    title: str = "",
    show_status: bool = True
) -> int | None:
    """
    Display interactive menu with arrow key navigation.
    Returns selected index or None if cancelled.
    """
```

### 2.4 base_view.py - Abstract Base
```python
class BaseView(ABC):
    """Base class for all view modules."""

    def __init__(self, controller: CliController):
        self.controller = controller
        self.console = controller.console

    @abstractmethod
    def render(self) -> None:
        """Render the view content."""

    @abstractmethod
    def get_menu_options(self) -> list[MenuOption]:
        """Return menu options for this view."""

    def show_in_development(self, feature: str) -> None:
        """Show [IN DEVELOPMENT] message."""
```

---

## Phase 3: View Implementations

### 3.1 CALC View (calc_view.py)
**Features:**
- Single molecule SMILES input
- Mock prediction display (HOF value, confidence, model)
- Option to predict another or go back

**UI Flow:**
```
[CALC] Molecular Property Prediction
─────────────────────────────────────
Enter SMILES: CCO

┌─ Prediction Result ─────────────────────┐
│ SMILES:      CCO (Ethanol)              │
│ Property:    Heat of Formation (HOF)    │
│ Predicted:   -56.21 kcal/mol            │
│ Confidence:  95.2%                      │
│ Model:       DeltaRF_v1.0               │
└─────────────────────────────────────────┘

❯ Predict another molecule
  Back to main menu
```

### 3.2 DATABASES View (databases_view.py)
**Features:**
- List available databases with status
- Show summary for each (molecules count, last updated)
- [IN DEVELOPMENT] for: Calculate, Edit, Delete, Add new

**Mock Data:**
```python
DATABASES = [
    {"name": "CBS Reference", "molecules": 30026, "status": "ready", "updated": "2026-01-11"},
    {"name": "CREST PM7", "molecules": 30026, "status": "ready", "updated": "2026-01-12"},
    {"name": "NIST Experimental", "molecules": 0, "status": "in_development"},
]
```

### 3.3 MODELS View (models_view.py)
**Features:**
- List trained models with metrics
- Show details: algorithm, MAE, R², training date, hyperparameters
- [IN DEVELOPMENT] for: Train, Test, Compare

**Mock Data:**
```python
MODELS = [
    {"name": "DeltaRF_v1.0", "algorithm": "Random Forest", "mae": 2.34, "r2": 0.945, "date": "2026-01-12"},
    {"name": "DeltaXGB_v1.0", "algorithm": "XGBoost", "mae": 2.12, "r2": 0.952, "date": "2026-01-12"},
    {"name": "DeltaNN_v1.0", "algorithm": "Neural Network", "status": "in_development"},
]
```

### 3.4 RESULTS View (results_view.py)
**Features:**
- Model performance comparison table
- CBS vs PM7 divergence analysis
- Accuracy metrics by severity bands

**Sections:**
1. Model Comparison Table (MAE, R², training time)
2. Divergence Distribution (LOW/MEDIUM/HIGH/CRITICAL)
3. Accuracy by Severity

### 3.5 SETTINGS View (settings_view.py)
**Status:** [IN DEVELOPMENT]
- Show message and return to main menu

### 3.6 ABOUT View (about_view.py)
**Features:**
- Version: 1.0.0-beta
- Description of GRIMPERIUM
- System status (databases, models ready)
- Links (documentation, GitHub)

---

## Phase 4: Integration & Testing

### 4.1 Integration with Existing Code
```python
# calc_view.py - Future integration point
from grimperium.data import ChemperiumLoader
from grimperium.models import DeltaLearningEnsemble
from grimperium.config import GrimperiumConfig

# For MVP: Use mock predictions
# For production: Load actual model and data
```

### 4.2 Test Strategy
- Unit tests for each view class
- Integration test for navigation flow
- Manual testing for keyboard interaction

---

## Phase 5: Implementation Order

### Step 1: Foundation (styles.py, controller.py, menu.py)
- [ ] Create `src/grimperium/cli/` directory
- [ ] Implement `styles.py` with colors and ASCII banner
- [ ] Implement `controller.py` with state management
- [ ] Implement `menu.py` with arrow key navigation
- [ ] Create `__init__.py` and `__main__.py`

### Step 2: Base Infrastructure (base_view.py, mock_data.py)
- [ ] Implement `BaseView` abstract class
- [ ] Create `mock_data.py` with all mock datasets
- [ ] Create `views/__init__.py`

### Step 3: Views Implementation
- [ ] Implement `about_view.py` (simplest, good for testing)
- [ ] Implement `databases_view.py` (tables, lists)
- [ ] Implement `models_view.py` (similar structure)
- [ ] Implement `results_view.py` (complex tables)
- [ ] Implement `calc_view.py` (user input handling)
- [ ] Implement `settings_view.py` ([IN DEVELOPMENT])

### Step 4: Main Application (app.py)
- [ ] Implement welcome screen with ASCII art
- [ ] Implement main menu loop
- [ ] Wire all views together
- [ ] Handle Ctrl+C gracefully

### Step 5: Entry Point & Testing
- [ ] Update `pyproject.toml` with new entry point
- [ ] Create `__main__.py` for `python -m grimperium.cli`
- [ ] Write unit tests for views
- [ ] Manual E2E testing

---

## Key Design Decisions

### Navigation Library Choice
**Decision:** `questionary` (cross-platform, built on prompt_toolkit)

**Rationale:**
- Works on Windows, Linux, and macOS
- Rich integration via prompt_toolkit styling
- Well-maintained, good documentation
- Supports: select, checkbox, confirm, text input

### Mock Data Strategy
- All mock data centralized in `mock_data.py`
- Easy to swap for real data later
- Clear separation of concerns

### Error Handling
- Graceful Ctrl+C handling at all levels
- [IN DEVELOPMENT] features show message, don't crash
- Invalid input re-prompts user

---

## Verification Checklist

- [ ] App starts with welcome screen
- [ ] ASCII banner renders correctly
- [ ] Main menu shows all 6 modules
- [ ] Arrow key navigation works
- [ ] Enter selects menu item
- [ ] Ctrl+C returns to previous menu or exits
- [ ] CALC accepts SMILES input and shows result
- [ ] DATABASES lists all databases with stats
- [ ] MODELS lists all models with metrics
- [ ] RESULTS shows comparison tables
- [ ] SETTINGS shows [IN DEVELOPMENT] message
- [ ] ABOUT shows version and info
- [ ] No unhandled exceptions
- [ ] Runs with `python -m grimperium.cli`

---

## Files to Create

| File | Lines (est.) | Priority |
|------|--------------|----------|
| `cli/__init__.py` | 10 | P1 |
| `cli/__main__.py` | 15 | P1 |
| `cli/styles.py` | 80 | P1 |
| `cli/controller.py` | 100 | P1 |
| `cli/menu.py` | 80 | P1 |
| `cli/mock_data.py` | 100 | P1 |
| `cli/app.py` | 150 | P1 |
| `cli/views/__init__.py` | 10 | P2 |
| `cli/views/base_view.py` | 60 | P2 |
| `cli/views/about_view.py` | 80 | P2 |
| `cli/views/databases_view.py` | 120 | P2 |
| `cli/views/models_view.py` | 120 | P2 |
| `cli/views/results_view.py` | 150 | P2 |
| `cli/views/calc_view.py` | 130 | P2 |
| `cli/views/settings_view.py` | 40 | P3 |

**Total estimated:** ~1,245 lines of code

---

## Dependencies Summary

```toml
# pyproject.toml additions
[tool.poetry.dependencies]
rich = ">=13.0"
questionary = ">=2.0.0"

[tool.poetry.scripts]
grimperium-cli = "grimperium.cli:main"
```
