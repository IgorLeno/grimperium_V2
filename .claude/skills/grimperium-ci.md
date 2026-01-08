---
name: grimperium-ci-fix
description: Analisa CI Error Summary Report e corrige todos os erros automaticamente
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

# 🔧 Skill: Fix Grimperium CI Errors

**Propósito:** Automatizar correção de erros de CI/CD

**Quando usar:** Quando você recebe um CI Error Summary Report do GitHub Actions

## Limitações / intervenção manual

- Se houver erros que exigem decisão humana, eu aviso
- Type errors podem exigir edição manual de signature de função
- Testes podem requerer mudança de lógica (não apenas formatting)

## O que esta skill faz

1. ✅ Parse do CI Error Report
2. ✅ Identifica tipo de erro (Lint, Type, Tests)
3. ✅ Aplica fixes automáticos principalmente para Lint e erros simples (ex.: formatting)
4. ✅ Para Type e Tests, pode apenas sugerir fixes e/ou aplicar patches diagnósticos; eu mostro os resultados de validação e posso pausar para decisão do dev e edições manuais antes de seguir (ver nota em Type e seção "Notas")
5. ✅ Commit + Push somente após revisão humana quando houver Type/Tests que exigem decisão/edição manual (posso parar antes do commit/push para você ajustar e confirmar)

## Tipos de Erros Suportados

### Lint Errors (Black Format)

Arquivos desformatados são rewritados automaticamente:

```bash
poetry run black src/ tests/
```

### Type Errors (Mypy)

Quando o job `typecheck` falha, siga um fluxo **iterativo**: isole o erro, aplique um fix pequeno, re-execute até ficar limpo.

- **Rodar mypy no projeto (baseline do CI)**:

```bash
poetry run mypy src/grimperium --ignore-missing-imports
```

- **Rodar mypy em um único arquivo (para iterar rápido)**:

```bash
poetry run mypy src/grimperium/caminho/do_arquivo.py --ignore-missing-imports
```

- **Padrões comuns e como corrigir**
  - **Missing return types / return implícito**: mypy reclama de `Function is missing a return type annotation` ou detecta `Any`.
    - **Fix**: adicione anotação explícita e garanta `return` em todos os caminhos.
  - **Union/Optional misuse**: `Item "None" of "X | None" has no attribute ...` (ou acesso a atributo/método sem checar `None`).
    - **Fix**: use `Optional[T]` (ou `T | None`) + guard clause (`if x is None: ...`) antes do uso.
  - **Incorrect imports / stubs ausentes**: bibliotecas sem type hints geram `Skipping analyzing ...`/`Library stubs not installed`.
    - **Fix**: instale stubs quando existirem, ajuste imports, ou isole o uso com `typing.cast(...)` / `# type: ignore[<code>]` **apenas** no ponto mínimo necessário (com comentário curto do motivo).
  - **Wrong attribute / container types**: `Incompatible types in assignment` / `Argument 1 ... has incompatible type ...`.
    - **Fix**: corrija assinatura/retorno (fonte da verdade), ajuste tipos de atributos/estruturas, e evite “forçar” com cast sem necessidade.

- **Exemplos de fixes (padrões frequentes)**
- **Adicionar return annotation explícita**:

```python
def parse_user_id(raw: str) -> int:
    return int(raw)
```

- **Optional + checagem antes do uso**:

```python
from typing import Optional

def normalize_name(name: Optional[str]) -> str:
    if name is None:
        return ""
    return name.strip()
```

- **Ajustar assinatura/retorno para bater com o uso real**:

```python
from collections.abc import Sequence

def get_ids(values: Sequence[str]) -> list[int]:
    return [int(v) for v in values]
```

- **Cast/ignore localizado (quando inevitável)**:

```python
from typing import cast

raw = get_untyped_value()
value = cast(str, raw)  # cast localizado: API externa não tipada
```

- **Workflow recomendado**
  - Rode mypy **no arquivo com erro**, aplique o menor patch possível, re-rote até zerar.
  - Ao final, rode novamente o mypy do pacote (`src/grimperium`) para garantir que não surgiram efeitos colaterais.

### Test Errors (Pytest)

Quando o job `tests` falha, o objetivo é **reproduzir localmente**, entender o traceback e corrigir de forma incremental.

- **Rodar a suíte (baseline do CI)**:

```bash
poetry run pytest tests/ -v --tb=short
```

- **Reproduzir um teste específico (rápido e determinístico)**:

```bash
poetry run pytest -k "<test_name>" -q
```

- **Inspecionar o traceback**:
  - Comece com `--tb=short` (mais legível).
  - Se precisar de mais contexto, remova `--tb=short` ou use `-vv`.

- **Debug interativo com PDB (quando o erro não está óbvio)**:

```bash
poetry run pytest -k "<test_name>" -q --pdb
```

- **Fixes comuns**
  - **Ajustar assertions frágeis**: normalize dados (ordem, timezone, floats), ou valide apenas o que importa.
  - **Mockar dependências externas** (rede/FS/tempo/UUID): use `monkeypatch`, fixtures e fakes para evitar flakiness.
  - **Setup/teardown incorreto em `conftest.py`**: revise fixtures com escopo errado (ex.: `session` vs `function`), estado global vazando, ou cleanup faltando.

- **Workflow recomendado (iterativo)**
  - Edite **o teste ou o source**, rode **um único teste** com `-k`, repita até passar.
  - Depois rode a suíte completa (`tests/`) para garantir que não quebrou nada.
  - Se ficar incerto, adicione um **TODO focado** (ex.: “TODO: investigar flakiness em {condição} / mock de `API`”) e abra uma issue com o traceback e passos de reprodução.

## Como Usar

### Opção 1: Colar CI Error Summary

```text
@claude /grimperium-ci-fix

Aqui está o CI Error Report:

# CI/CD Error Summary Report

Generated: 2026-01-07 09:38:30 UTC
Commit: `738992bf4f69`
Branch: `main`
Run: #15

***

## Overall Status

❌ FAILURES DETECTED

***

## Run Metadata

- Workflow: `CI`
- Job(s): `lint`, `typecheck`, `tests`
- Trigger: `push` (branch `main`)
- Actor: `igor`
- Runner: `ubuntu-latest`
- Python: `3.12`
- Poetry: `1.8.x`

***

## Summary

- Overall: ❌ FAILED
- Failed jobs: 1 / 3
- Failed steps: 1
- Total duration: 06m 42s

***

## Per-job Breakdown

### Job: `lint`

- Status: ✅ PASSED
- Duration: 00m 41s

### Job: `typecheck`

- Status: ✅ PASSED
- Duration: 01m 12s

### Job: `tests`

- Status: ❌ FAILED
- Duration: 04m 26s

#### Failed step: `Run pytest`

- Command:
  - `poetry run pytest tests/ -q --maxfail=1 --disable-warnings`
- Error (excerpt):

```text
=================================== FAILURES ===================================
_____________________ test_full_pipeline_real_dataset __________________________
tests/integration/test_full_pipeline.py:84: in test_full_pipeline_real_dataset
    result = pipeline.run(dataset="real")
src/grimperium/api.py:132: in run
    data = load_dataset(dataset=dataset, **kwargs)
src/grimperium/data/loader.py:211: in load_dataset
    raise FileNotFoundError(f"Arquivo de dataset não encontrado: {path}")
E   FileNotFoundError: Arquivo de dataset não encontrado: /home/runner/work/grimperium/grimperium/data/real/thermo_cbs_opt.csv
=========================== short test summary info ============================
FAILED tests/integration/test_full_pipeline.py::test_full_pipeline_real_dataset
1 failed, 87 passed in 31.42s
ERROR: Process completed with exit code 1.
```

#### Logs relevantes (tail)

```text
[grimperium] dataset=real
[grimperium] resolved path=/home/runner/work/grimperium/grimperium/data/real/thermo_cbs_opt.csv
[grimperium] hint: configure GRIMPERIUM_DATA_DIR or use fixtures/real_data.py
```

***

## Failed Tests Summary

- Suite: `tests/integration`
- Failing tests (1):
  - `tests/integration/test_full_pipeline.py::test_full_pipeline_real_dataset`

### Stack trace (short)

```text
FileNotFoundError: Arquivo de dataset não encontrado: /home/runner/work/grimperium/grimperium/data/real/thermo_cbs_opt.csv
  at src/grimperium/data/loader.py:211 in load_dataset
  called by src/grimperium/api.py:132 in run
  called by tests/integration/test_full_pipeline.py:84 in test_full_pipeline_real_dataset
```

***

## Artifacts / Links

- GitHub Run: `https://github.com/<org>/<repo>/actions/runs/123456789`
- Job logs:
  - `tests`: `https://github.com/<org>/<repo>/actions/runs/123456789/job/987654321`
- Artifacts (expected):
  - `pytest-report.xml` (missing)
  - `coverage.xml` (present)
  - `ci-error-summary.md` (present)

***

## Recommended Next Steps

- Verificar se o CI está baixando/coplando o dataset real no job `tests` (ou se o teste deveria usar fixtures mockadas).
- Confirmar o caminho esperado do arquivo (`thermo_cbs_opt.csv`) e alinhar com `GRIMPERIUM_DATA_DIR` / `src/grimperium/data/loader.py`.
- Reproduzir localmente:
  - `poetry install`
  - `poetry run pytest tests/integration/test_full_pipeline.py -v --tb=short`
- Se o dataset real não deve rodar no CI, marcar o teste com skip condicional e garantir cobertura via fixtures.

### Opção 2: Ler arquivo diretamente

```text
@claude /grimperium-ci-fix CI_ERROR_SUMMARY.md
```

## Processo Passo a Passo

```text
1️⃣ Parse Error Report
   ├─ Lint errors (Black)? → black src/ tests/
   ├─ Type errors (Mypy)? → Listar e corrigir
   └─ Test errors? → Analisar traceback

2️⃣ Execute Fixes
   ├─ poetry run black src/ tests/
   ├─ poetry run mypy src/ --ignore-missing-imports
   └─ poetry run pytest tests/ -v

3️⃣ Validate
   ├─ Ruff: ruff check src/ tests/
   ├─ Mypy: mypy src/grimperium
   └─ Pytest: 88 passed, 0 errors

4️⃣ Commit + Push
   ├─ git add -A
   ├─ git commit -m "Fix: Resolve CI errors - format + type + tests"
   └─ git push origin main
```

## Output

- ✅ Lista de todos os erros corrigidos
- ✅ Comando de cada fix executado
- ✅ Resultado final (pass/fail)
- ✅ Commit hash se bem-sucedido

## Notas

Ver **Limitações / intervenção manual** (acima).

## Next Steps Quando Falhar

Se a skill não conseguir corrigir:

```bash
# Você executa manualmente
poetry run black src/ tests/
poetry run mypy src/grimperium --ignore-missing-imports

# Então pede ajuda
@claude analyze my errors [copy error messages]
```
