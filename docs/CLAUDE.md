# Claude Code Configuration - Grimperium Project
## Phase C: CLI Interactive Application - BATCH 12 Focus

---

## 🎯 PROJECT CONTEXT

**Project Name:** Grimperium Delta-Learning Framework  
**Current Phase:** C - CLI Interactive Application  
**Current BATCH:** 12 (Critical Fixes)  
**Deadline:** BATCH 12 smoke test completion (3 molecules)  
**Stack:** Python 3.10+, Rich, Questionary, Pytest, Ruff, Black, MyPy  

### Quick Stats

```
Lines of Code: ~3,500 (src/)
Test Coverage Target: 85%+
Type Hints: 100% required
CLI Modules: 7 (CALC, DATABASES, MODELS, RESULTS, SETTINGS, ABOUT, BATCH)
Open Bugs: 11 (BATCH 12)
Estimated Fix Time: 90-120 minutes
```

---

## 📊 GRIMPERIUM-SPECIFIC AGENTS

These agents activate automatically in Grimperium context:

### 1. **Python Development Agent** (Core)

**Activated When:** Any Python file mentioned

**Grimperium-Specific Knowledge:**
- Type hints 100% requirement (mypy --strict)
- Async patterns in CLI event loops
- Rich library integration patterns
- Questionary interactive flows
- Pandas/NumPy for molecular data

**Example Activation:**
```
"Implementa csv_manager.py com type hints 100%"
→ Agent knows: Must pass mypy --strict, use Pandas properly
```

### 2. **Testing Expert Agent** (Pytest Focus)

**Activated When:** Tests, coverage, or pytest mentioned

**Grimperium-Specific Knowledge:**
- Pytest fixtures for molecular data
- Mock CREST results
- Coverage target: 85%+
- Test structure mirrors src/ structure
- Molecular property assertions

**Example Activation:**
```
"Escreve 15 testes para database_view.py com cov > 85%"
→ Agent knows: Must mock CSV loading, verify counts, test update flow
```

### 3. **Code Review Expert Agent** (Quality Gate)

**Activated When:** Code review requested

**Grimperium-Specific Quality Criteria:**
1. **Type Hints** (100%)
2. **Test Coverage** (85%+)
3. **Linting** (Ruff clean)
4. **Formatting** (Black aligned)
5. **Molecular Correctness** (math validated)
6. **CLI UX** (user-friendly)

**Example Activation:**
```
"/superpowers:code-review src/cli/views/database_view.py"
→ Reviews all 6 criteria above
```

### 4. **CLI Development Agent** (New - Specialized)

**Activated When:** CLI, menu, view, or user-interaction mentioned

**Grimperium-Specific Knowledge:**
- Rich library (panels, tables, console output)
- Questionary patterns (choices, prompts)
- Navigation flow (main menu → submenu → action → back)
- Error handling for user input
- Progress indicators for long operations

**Example Activation:**
```
"Adiciona botão Back to Main Menu em Settings"
→ Agent knows: Rich patterns, Questionary navigation, proper state management
```

### 5. **Superpowers Suite** (Always Available)

**Commands:**
```
/superpowers:brainstorm    → Design CLI features, UX flows
/superpowers:write-plan    → Create BATCH plans with time estimates
/superpowers:execute-plan  → Run multi-step implementation
/superpowers:debug         → Debug CLI behavior
/superpowers:code-review   → Quality gate before commit
```

---

## 🎨 GRIMPERIUM-SPECIFIC SLASH COMMANDS

Add these to your Claude Code config (custom commands):

### BATCH Management

```
/grimperium:batch-start
  Input: BATCH number, title, bug list
  Output: Full BATCH context, time estimates, code-before/after
  Example: /grimperium:batch-start 12 "CLI Critical Fixes"

/grimperium:batch-audit
  Input: Bugs reported by user
  Output: Complete audit, categorized, prioritized
  Example: "Testei e encontrei bugs X, Y, Z" → /grimperium:batch-audit

/grimperium:batch-execute
  Input: BATCH number
  Output: Step-by-step execution with validation
  Example: /grimperium:batch-execute 12
```

### Testing & Quality

```
/grimperium:test-coverage
  Input: Module name (optional)
  Output: Coverage report, missing tests highlighted
  Example: /grimperium:test-coverage database_view

/grimperium:test-all
  Input: None
  Output: Full test suite + coverage report
  Example: /grimperium:test-all

/grimperium:quality-gate
  Input: File or module
  Output: Passes: type hints, linting, formatting, coverage
  Example: /grimperium:quality-gate src/cli/views/
```

### CLI Testing & Validation

```
/grimperium:cli-test
  Input: Module to test (CALC, DATABASES, MODELS, RESULTS, SETTINGS, ABOUT)
  Output: Manual testing flow with expected outputs
  Example: /grimperium:cli-test DATABASES

/grimperium:cli-flow
  Input: None
  Output: Full CLI user flow diagram (ASCII)
  Example: /grimperium:cli-flow

/grimperium:validate-molecules
  Input: Count or range
  Output: Loads test molecules, validates counts, outputs summary
  Example: /grimperium:validate-molecules 3
```

### Documentation

```
/grimperium:update-docs
  Input: Section to update (PHASE-A, PHASE-B, CLAUDE, architecture)
  Output: Updated documentation file
  Example: /grimperium:update-docs CLAUDE

/grimperium:docs-sync
  Input: None
  Output: Checks docs against code, updates mismatches
  Example: /grimperium:docs-sync
```

### Debug & Fix

```
/grimperium:bug-details
  Input: Bug ID or description
  Output: Bug details, file location, code context
  Example: /grimperium:bug-details "CBS Original database"

/grimperium:fix-suggest
  Input: Bug description
  Output: Multiple fix options, pros/cons
  Example: /grimperium:fix-suggest "Update button not working"
```

---

## 🔄 GRIMPERIUM WORKFLOW (Phase C - BATCH Model)

### Pre-BATCH Workflow (Planning)

```
1. User Tests CLI
   "Testei e encontrei bugs X, Y, Z"
   ↓
2. Auditor (Human) Makes Complete Audit
   → Categorizes bugs (CRITICAL, HIGH, MEDIUM)
   → Estimates time per bug
   → Proposes ordering
   ↓
3. /superpowers:write-plan
   → Generate BATCH X with specifications
   → Code BEFORE/AFTER for each fix
   → Success criteria
   ↓
4. User Validates Plan
   → Asks questions
   → Approves or adjusts
   ↓
5. /grimperium:batch-start X
   → Initialize BATCH with full context
```

### Per-BATCH Workflow (Implementation)

```
1. git checkout -b fix/batch-X
   ↓
2. /grimperium:batch-execute X
   → Step-by-step implementation
   → Follows code BEFORE/AFTER specifications
   ↓
3. /grimperium:test-coverage
   → Run tests after each change
   → Ensure coverage maintained (85%+)
   ↓
4. /grimperium:quality-gate
   → Type hints ✅
   → Linting ✅
   → Formatting ✅
   → Tests ✅
   ↓
5. /grimperium:cli-test [MODULE]
   → Manual testing of changed module
   ↓
6. /superpowers:code-review
   → Final 6-aspect review
   ↓
7. git commit -m "fix(batch-X): [description]"
   ↓
8. git push origin fix/batch-X
```

### Post-BATCH Workflow (Validation)

```
1. Manual CLI smoke test (3 molecules)
   ↓
2. /grimperium:validate-molecules 3
   → Confirms molecules loaded correctly
   ↓
3. Merge branch: git merge --ff-only fix/batch-X
   ↓
4. Delete branch: git branch -d fix/batch-X
   ↓
5. Report results to user
   ↓
6. Plan next BATCH (if needed)
```

---

## 📁 CRITICAL FILE REFERENCES

These files are context for agents (mention them for better understanding):

### Core CLI Structure
```
src/grimperium/cli/
├─ __init__.py
├─ __main__.py                   (Module entry point)
├─ app.py                        (CLI application)
├─ controller.py                (Main orchestrator)
├─ menu.py                      (Menu logic)
├─ styles.py                    (Rich styling)
├─ settings_manager.py          (Settings persistence)
└─ views/
   ├─ base_view.py             (Abstract view)
   ├─ calc_view.py             (Molecular calculation)
   ├─ databases_view.py        ← BUGS: CBS, count, flow, update
   ├─ models_view.py           (ML model management)
   ├─ results_view.py          (Analytics display)
   ├─ settings_view.py         ← BUGS: headers, back button
   ├─ about_view.py            (Info display)
   └─ batch_view.py            (Batch processing management)

Data Management
src/grimperium/crest_pm7/
├─ batch/
│  ├─ csv_manager.py           ← BUG: wrong count
│  └─ processor_adapter.py     ← BUG: wrong flow

Test Structure (85%+ coverage required)
tests/
├─ unit/
│  ├─ test_calc_view.py
│  ├─ test_databases_view.py    ← Needs updates for bugs
│  ├─ test_settings_view.py     ← Needs updates for bugs
│  └─ test_settings_manager.py
```

---

## 🐛 BATCH 12: THE 11 BUGS (Quick Reference)

For context when discussing bugs:

| # | Severity | Issue | File | Status |
|---|----------|-------|------|--------|
| 1 | 🔴 CRITICAL | CBS Original (deleted file) | database_view.py | 👀 Ready |
| 2 | 🔴 CRITICAL | CREST PM7 wrong count | csv_manager.py | 👀 Ready |
| 3 | 🔴 BLOCKER | Calc flow broken | database_view.py | 👀 Ready |
| 4 | 🟡 HIGH | Update button doesn't work | database_view.py | 👀 Ready |
| 5 | 🟡 HIGH | Duplicate headers in Settings | settings_view.py | 👀 Ready |
| 6 | 🔴 BLOCKER | No back button in Settings | settings_view.py | 👀 Ready |
| 7 | 🟢 MEDIUM | MOPAC title misaligned | settings_manager.py | 👀 Ready |
| 8 | 🟡 HIGH | xTB should be inside CREST | settings_manager.py | 👀 Ready |
| 9 | 🟢 MEDIUM | CREST toggles not intuitive | settings_manager.py | 👀 Ready |
| 10 | 🟢 MEDIUM | MOPAC settings not intuitive | settings_manager.py | 👀 Ready |
| 11 | 🟢 MEDIUM | Missing visual feedback | General | 👀 Ready |

---

## ✅ QUALITY REQUIREMENTS (Non-Negotiable)

Every change in Grimperium must meet:

```
1. TYPE HINTS (100%)
   → mypy --strict src/ (zero errors)
   
2. LINTING (Ruff clean)
   → ruff check src/ (zero errors)
   
3. FORMATTING (Black)
   → black --check src/ (zero changes)
   
4. TESTS (85%+ coverage)
   → pytest tests/ --cov=src/ --cov-report=term
   → Minimum 85% line coverage
   
5. CODE REVIEW (6 aspects)
   → /superpowers:code-review [file]
   → Readability, Performance, Security, Maintainability, etc.
   
6. MOLECULAR CORRECTNESS
   → Calculations validated against research
   → CSV counts accurate
   → Data flows correctly
```

---

## 🎯 SUCCESS CRITERIA (BATCH 12)

BATCH 12 is done when:

- ✅ All 11 bugs fixed and tested
- ✅ CLI smoke test passes (3 molecules load correctly)
- ✅ Settings menu works: back button, proper layout, no dupes
- ✅ CREST PM7 shows correct count (3, not wrong number)
- ✅ Update button functional
- ✅ xTB inside CREST section
- ✅ All tests pass (85%+ coverage maintained)
- ✅ Code review: all 6 aspects pass
- ✅ Zero mypy/ruff/black errors
- ✅ Documentation updated (CLAUDE.md, PHASE-A-START-HERE)

---

## 🚀 QUICK START (Per Session)

When starting work on Grimperium:

```
1. Check current BATCH
   /serena:list-memories | grep "current-batch"

2. Read BATCH specifications
   /grimperium:batch-start [number]

3. Create feature branch
   git checkout -b fix/batch-X

4. Implement changes (following BEFORE/AFTER code)
   /grimperium:batch-execute X

5. Test thoroughly
   /grimperium:test-coverage
   /grimperium:quality-gate

6. Review your code
   /superpowers:code-review

7. Commit and push
   git commit -m "fix(batch-X): [description]"
   git push origin fix/batch-X

8. Validate manually
   /grimperium:cli-test [module]
   /grimperium:validate-molecules 3

9. Plan next BATCH (if applicable)
```

---

## 📚 REFERENCE DOCUMENTATION

| Document | Purpose | Location |
|----------|---------|----------|
| **architecture.md** | System design, data flow, Phase A/B/C | /docs/ |
| **PHASE-A-START-HERE.md** | Phase A validation, success metrics | /docs/ |
| **CHANGELOG.md** | Version history, what changed | /root/ |
| **README.md** | Project overview, quick start | /root/ |
| **CLAUDE.md (Global)** | Universal best practices | /root/ |

---

## 🔐 IMPORTANT RULES FOR GRIMPERIUM

1. **Always TDD** - Write tests first, implementation second
2. **Type Hints First** - 100% before implementation
3. **No Vague Prompts** - Always mention file names and line numbers
4. **Test Before Commit** - Every commit must pass `/grimperium:quality-gate`
5. **Molecular Correctness** - Math validated against research
6. **CLI Flow** - User can always return to main menu (back button)
7. **Documentation Sync** - Update CLAUDE.md if workflow changes
8. **Batch Model** - Work in BATCHes, not ad-hoc changes

---

## 🎓 WHEN TO CALL THE ORIENTADOR (Human)

These situations require human decision-making:

- New BATCH discovered during testing (call for audit)
- Design question (CLI structure, workflow)
- Trade-off decision (speed vs. quality, scope vs. time)
- Phase transition (A → B, B → C)
- Technology decision (add new library, change architecture)
- Priority reordering (reorder bugs in BATCH)

**How to call:**
```
"Preciso de uma auditoria em [module].
Encontrei: [X bug], [Y issue], [Z question]"

→ Orientador responds with:
   - Complete audit
   - Categorization
   - Time estimates
   - Prioritization
   - Recomendações
```

---

**Version:** 1.0 (Phase C - BATCH 12)
**Last Updated:** 2026-01-17  
