# Grimperium Workflow Reference

This is a supporting workflow document. The local authority for agent behavior
is [`/CLAUDE.md`](/home/igor/Projetos/grimperium/CLAUDE.md).

## Default Flow

1. Clarify scope and plan non-trivial work.
2. Read the relevant code, tests, and adjacent docs.
3. Implement the smallest coherent change.
4. Run targeted validation while iterating.
5. Run the broadest relevant gate before handoff.
6. Update `CHANGELOG.md` and affected docs when workflow-visible behavior
   changed.

## Repo-Native Verification Commands

```bash
black src/ tests/
ruff check src/ tests/
mypy src/ --strict
pytest tests/ -v --cov=src/grimperium
```

## Batch Habit

- Use batch language for coherent non-trivial work, but do not assume a fixed
  historical batch number.
- Keep `[Unreleased]` in `CHANGELOG.md` current for completed batches or
  meaningful workflow-visible changes.

## Related Files

- Main contract:
  [`CLAUDE.md`](/home/igor/Projetos/grimperium/CLAUDE.md)
- Execution details:
  [`.claude/DEVELOPMENT_PROCEDURES.md`](/home/igor/Projetos/grimperium/.claude/DEVELOPMENT_PROCEDURES.md)
- Testing reference:
  [`docs/TESTING_GUIDE.md`](/home/igor/Projetos/grimperium/docs/TESTING_GUIDE.md)
- Historical transition note:
  [`docs/.archive/deprecated/2026-03-20-instruction-system-transition.md`](/home/igor/Projetos/grimperium/docs/.archive/deprecated/2026-03-20-instruction-system-transition.md)
