# CHANGELOG Entry Template

Use this template when updating `CHANGELOG.md` under `[Unreleased]`.

## Minimal Template

```markdown
### Changed
- **<scope or batch label>** (YYYY-MM-DD)
  - `path/to/file.py`: concise summary of what changed
  - Verification: `pytest path/to/test.py -v`

### Fixed
- **<scope or batch label>** (YYYY-MM-DD)
  - Fixed <bug or regression summary>
  - `path/to/file.py`: note the corrected behavior
```

## Expanded Template

```markdown
### Added
- **<scope or batch label>** (YYYY-MM-DD)
  - `path/to/new_file.py`: what was added
  - `tests/path_to_test.py`: coverage added for the new behavior

### Changed
- **<scope or batch label>** (YYYY-MM-DD)
  - `path/to/file.py`: what changed and why

### Fixed
- **<scope or batch label>** (YYYY-MM-DD)
  - Fixed <user-visible or workflow-visible issue>
  - Verification: `pytest ...`, `mypy ...`, `ruff ...`

### Deprecated
- **<scope or batch label>** (YYYY-MM-DD)
  - `old_api_or_file`: reason
  - Replacement: `new_api_or_file`

### Removed
- **<scope or batch label>** (YYYY-MM-DD)
  - Removed <obsolete file, path, or behavior>
```

## Checklist

- Entry is under `[Unreleased]`
- Date uses `(YYYY-MM-DD)`
- Section type matches the change
- Paths and code identifiers use backticks
- Verification is included when it helps the reader
- No speculative version number is added
