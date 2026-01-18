# Grimperium Development Workflow v2.0
## awesome-claude-code Integration + Phase C Ready

**Last Updated:** 2026-01-17  
**Status:** ✅ Phase C (CLI Interactive) - Ready to start BATCH 12

---

## 🎯 Overview

Grimperium segue padrões de **awesome-claude-code** para máxima eficiência:

1. **Specification-Driven Development** - ANTES/DEPOIS contracts
2. **Progressive Disclosure** - Agentes carregam sob demanda
3. **Session Continuity** - Serena memory persiste contexto
4. **Ralph Wiggum Loop** - Execução autônoma com validação
5. **TDD-First Culture** - Testes antes do código
6. **Continuous Code Review** - 6 aspectos validados
7. **Documentation Sync** - Docs sempre alinhadas com código

---

## 📋 WORKFLOW: BATCH Execution

### Phase 1: Audit & Planning
```
/grimperium:batch-audit 12
→ Identifica 11 bugs críticos
→ Categoriza por prioridade
→ Estima tempo por bug
```

### Phase 2: Specification First
Para **cada bug**:
```markdown
## FIX: [Bug Title]

### ANTES (Current - Broken)
[código atual]

### DEPOIS (Target - Fixed)
[código corrigido]

### Tests
[testes de validação]

### Success Criteria
- [ ] Code matches DEPOIS
- [ ] All tests pass
- [ ] Coverage ≥ 85%
- [ ] Review passes 6/6
```

### Phase 3: Implementation (Ralph Wiggum Loop)
```
1. Implement DEPOIS
2. Run tests → RED?
   - Fix → GREEN
3. Code review → Fail?
   - Fix → Pass
4. Update memory
5. Next bug
```

### Phase 4: Quality Gates
```
/grimperium:quality-gate
→ mypy --strict src/: 0 errors ✅
→ ruff check src/: 0 errors ✅
→ black --check src/: OK ✅
→ pytest --cov: ≥85% ✅
→ CLI smoke test: OK ✅
```

### Phase 5: Docs & Commit
```
/grimperium:docs-sync
/grimperium:changelog-entry "Phase C: resolve 11 critical CLI bugs"
git commit -m "fix(batch-12): [message]"
```

---

## 🧠 Progressive Agent Loading

### Always Available (Core)
- `/superpowers:brainstorm` - Design & exploration
- `/superpowers:write-plan` - Planning
- `/superpowers:execute-plan` - Implementation
- `/superpowers:debug` - Debugging
- `/superpowers:code-review` - Quality check

### Load on Context
| Context Keyword | Agent | Cost |
|-----------------|-------|------|
| "test", "pytest" | Testing Expert | 180 tokens |
| "type", "async" | Python Expert | 150 tokens |
| "review", "quality" | Review Expert | 200 tokens |
| "cli", "menu" | CLI Expert | 160 tokens |
| "docs", "changelog" | Docs Expert | 140 tokens |

**Token Optimization:** ~70% savings vs. monolithic loading

---

## 💾 Serena Memory Pattern

### Per-BATCH Memory
```bash
/serena:write-memory "batch-12-status: {
  total: 11,
  fixed: 0,
  current: 'audit-phase',
  start: '2026-01-17T18:51',
  target-completion: '2026-01-17T21:30'
}"
```

### Update After Each Bug
```bash
/serena:edit-memory "batch-12-status: {
  ...,
  fixed: 1,
  remaining: [2,3,4,5,6,7,8,9,10,11],
  last-fix: 'Bug #1 (CBS Original)',
  time-spent: '10 min'
}"
```

---

## 🚀 Day-to-Day Commands

### Session Start
```bash
cd /home/igor/Projetos/grimperium
/serena:list-memories | grep batch-12
→ See exact status, resume from there
```

### Start BATCH
```bash
/grimperium:batch-execute 12
→ Runs Auditor → Planner → Executor pipeline
```

### Test Individual Module
```bash
/grimperium:test-module "database_view"
→ Runs tests for specific module only
```

### Code Review
```bash
/superpowers:code-review
→ 6-aspect review: type hints, tests, linting, format, perf, correctness
```

### Final Validation
```bash
/grimperium:quality-gate
→ All checks before merge
```

---

## 📊 Quality Standards (Non-Negotiable)

### Type Hints
```bash
mypy --strict src/
→ 0 errors required
```

### Tests
```bash
pytest tests/ --cov=src/ --cov-fail-under=85
→ ≥85% coverage required
```

### Linting
```bash
ruff check src/
→ 0 errors required
```

### Formatting
```bash
black --check src/
→ 0 changes required
```

### Code Review (6 Aspects)
1. **Type Hints** - 100% coverage
2. **Tests** - Comprehensive, ≥85% coverage
3. **Linting** - Ruff clean
4. **Formatting** - Black aligned
5. **Performance** - Optimized
6. **Correctness** - Validated

---

## 🔄 Git Hooks

### Pre-commit (Automatic)
Runs BEFORE commit:
```bash
✅ mypy --strict src/
✅ ruff check src/
✅ black --check src/
✅ pytest tests/
```

Blocks commit se qualquer check falha (use `--no-verify` só em emergências).

---

## 📚 Key Files

| File | Purpose | Last Updated |
|------|---------|--------------|
| `CLAUDE.md` v2.0 | Behavioral guide for Claude Code | 2026-01-17 |
| `.claude/settings.json` | Claude Code config (minimal, hooks ready) | 2026-01-17 |
| `.git/hooks/pre-commit` | Quality gate automático | 2026-01-17 |
| `VALIDATION_CHECKLIST.md` | Setup validation (antes de BATCH 12) | 2026-01-17 |
| `CHANGELOG.md` | Unreleased: Phase C setup | 2026-01-17 |
| `README.md` | Phase C overview | 2026-01-17 |
| `architecture.md` | v2.0 com CLI redesign | 2026-01-17 |

---

## 🎯 Success = Combination Of

1. **Clear SPEC** (ANTES/DEPOIS) ✅
2. **Tests First** (TDD) ✅
3. **Code Review** (6 aspects) ✅
4. **Quality Gates** (mypy, ruff, black, pytest) ✅
5. **Memory Persistence** (Serena) ✅
6. **Documentation Sync** (auto-update) ✅

**Follow this and Phase C completion = guaranteed.** 🚀

---

**Version:** 2.0 (awesome-claude-code integrated)  
**Status:** Ready for production  
**Next:** Execute VALIDATION_CHECKLIST.md, then start BATCH 12
