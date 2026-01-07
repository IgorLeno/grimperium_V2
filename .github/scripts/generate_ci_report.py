#!/usr/bin/env python3
"""
Gera um relatório consolidado do CI a partir dos logs baixados como artifacts.

Este script foi pensado para ser resiliente: caso algum log não exista
(job não rodou, artifact não foi gerado, etc), ele não deve quebrar o CI —
apenas reportar o estado.
"""
import os
import re
from pathlib import Path
from datetime import datetime, timezone

# Limites de exibição no relatório (para evitar reports enormes no GitHub UI).
MAX_LINT_ERRORS_DISPLAYED = 30
MAX_TYPE_ERRORS_DISPLAYED = 30
MAX_TEST_FAILURES_DISPLAYED = 20
MAX_TEST_ERRORS_DISPLAYED = 20

# Define paths
WORKSPACE = Path(os.environ.get("GITHUB_WORKSPACE", "."))
LOGS_DIR = WORKSPACE / "logs"
LINT_LOG = LOGS_DIR / "lint-log" / "lint.log"
TYPE_LOG = LOGS_DIR / "type-log" / "type.log"
TEST_LOGS = list(LOGS_DIR.glob("test-log-*/test-*.log"))
REPORT_FILE = WORKSPACE / "CI_ERROR_SUMMARY.md"
GITHUB_STEP_SUMMARY = Path(
    os.environ.get("GITHUB_STEP_SUMMARY", "/tmp/summary.md")
)

# Collect metadata
COMMIT = os.environ.get("GITHUB_SHA", "unknown")[:12]
BRANCH = os.environ.get("GITHUB_REF_NAME", "unknown")
RUN_NUMBER = os.environ.get("GITHUB_RUN_NUMBER", "unknown")
RUN_ID = os.environ.get("GITHUB_RUN_ID", "unknown")
REPO = os.environ.get("GITHUB_REPOSITORY", "unknown")
SERVER_URL = os.environ.get("GITHUB_SERVER_URL", "https://github.com")
TIMESTAMP = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def read_log(path):
    """
    Lê um arquivo de log de forma segura.

    Não deve explodir em caso de IO/encoding.
    """
    if path.exists():
        try:
            return path.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            return f"Error reading log: {e}"
    return "Log file not found"


def sanitize_log_excerpt(
    text: str, *, max_lines: int = 12, max_chars: int = 400
) -> str:
    """
    Retorna um trecho curto e sanitizado de um log para exibição em relatórios.

    Objetivo: evitar vazamento acidental de segredos (tokens/chaves) e
    paths locais.
    """
    if not text:
        return ""

    # Linhas que tendem a conter segredos ou material sensível
    # devem ser removidas.
    drop_line_patterns = [
        re.compile(r"(?i)-----BEGIN (?:RSA|OPENSSH|EC|PGP) PRIVATE KEY-----"),
        re.compile(r"(?i)\baws_secret_access_key\b"),
    ]

    # Substituições pontuais (redação inline) para preservar contexto do erro.
    redact_substitutions: list[tuple[re.Pattern[str], str]] = [
        # Tokens comuns (GitHub / genéricos)
        (re.compile(r"\bgh[pousr]_[A-Za-z0-9_]{20,}\b"), "<redacted-token>"),
        (re.compile(r"\bgithub_pat_[A-Za-z0-9_]{20,}\b"), "<redacted-token>"),
        # JWT (3 segmentos base64url)
        (
            re.compile(
                r"\beyJ[a-zA-Z0-9_-]{10,}\."
                r"[a-zA-Z0-9_-]{10,}\."
                r"[a-zA-Z0-9_-]{10,}\b"
            ),
            "<redacted-jwt>",
        ),
        # AWS Access Key Id
        (re.compile(r"\bAKIA[0-9A-Z]{16}\b"), "<redacted-aws-key-id>"),
        # Authorization headers
        (
            re.compile(r"(?i)(authorization:\s*bearer\s+)(\S+)"),
            r"\1<redacted>",
        ),
        (
            re.compile(r"(?i)(authorization:\s*token\s+)(\S+)"),
            r"\1<redacted>",
        ),
        # Assignments de variáveis sensíveis
        (
            re.compile(
                r"(?i)\b(token|api[_-]?key|secret|password|passwd)\b"
                r"(\s*[:=]\s*)"
                r"(?P<q>['\"])?"
                r"(?P<val>(?(q).*?(?=(?P=q))|\S+))"
                r"(?(q)(?P=q))"
            ),
            r"\1\2\g<q><redacted>\g<q>",
        ),
        # Paths absolutos comuns (Linux/macOS runners)
        (re.compile(r"(/(?:home|Users|runner|root)/[^\s:]+)"), "<path>"),
    ]

    sanitized_lines: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if not line.strip():
            continue

        if any(p.search(line) for p in drop_line_patterns):
            continue

        for pattern, replacement in redact_substitutions:
            line = pattern.sub(replacement, line)

        sanitized_lines.append(line)
        if len(sanitized_lines) >= max_lines:
            break

    excerpt = "\n".join(sanitized_lines).strip()
    if len(excerpt) > max_chars:
        excerpt = excerpt[: max_chars - 1].rstrip() + "…"
    return excerpt


def should_include_log_excerpt(status: str, text: str) -> bool:
    """
    Heurística simples: só mostra trecho quando houver indício claro de erro.
    """
    if status == "failed":
        return True
    if not text:
        return False
    return bool(re.search(r"\bERROR\b", text))


def parse_lint_log(text):
    """Parse ruff/black output"""
    # Padroniza ausência de log como "not_run"
    # (job não executou / artifact ausente).
    if "Log file not found" in text or "Error reading" in text:
        return [], "not_run"

    # Extract ruff errors
    ruff_errors = re.findall(
        r'^.+?:\d+:\d+: [A-Z]\d+.+$',
        text,
        re.MULTILINE
    )

    # Extract black formatting issues
    black_errors = re.findall(
        r'^would reformat .+$',
        text,
        re.MULTILINE
    )

    all_errors = ruff_errors + black_errors

    # Determine status
    if not all_errors:
        status = "passed"
    else:
        status = "failed"

    return all_errors, status


def parse_type_log(text):
    """Parse mypy output"""
    # Padroniza ausência de log como "not_run"
    # (job não executou / artifact ausente).
    if "Log file not found" in text or "Error reading" in text:
        return [], "not_run"

    # Extract mypy errors
    errors = re.findall(
        r'^.+?:\d+: error:.+$',
        text,
        re.MULTILINE
    )

    # Check for success message
    if "Success: no issues found" in text:
        status = "passed"
    elif errors:
        status = "failed"
    else:
        status = "not_run"

    return errors, status


def parse_test_log(text, python_version=""):
    """Parse pytest output"""
    # Padroniza ausência de log como "not_run"
    # (job não executou / artifact ausente).
    if "Log file not found" in text or "Error reading" in text:
        return [], "Log file not found", [], "not_run"

    # Extract summary line
    summary_match = re.search(
        r'=+ (.+ passed.+|.+ failed.+|.+ error.+) =+',
        text
    )
    summary = summary_match.group(1) if summary_match else "No summary found"

    # Extract failure details
    failures = re.findall(
        r'^FAILED .+$',
        text,
        re.MULTILINE
    )

    # Extract error details
    errors = re.findall(
        r'^ERROR .+$',
        text,
        re.MULTILINE
    )

    # Determine status
    if "passed" in text and not failures and not errors:
        status = "passed"
    elif failures or errors:
        status = "failed"
    else:
        status = "not found"

    return failures, summary, errors, status


# Read logs
print(f"📂 Reading logs from {LOGS_DIR}...")
lint_log = read_log(LINT_LOG)
type_log = read_log(TYPE_LOG)

# Parse
print("🔍 Parsing lint log...")
lint_errors, lint_status = parse_lint_log(lint_log)

print("🔍 Parsing type check log...")
type_errors, type_status = parse_type_log(type_log)

# Parse test logs for all Python versions
test_results = []
print("🔍 Parsing test logs...")
if TEST_LOGS:
    for test_log_path in TEST_LOGS:
        # Extract Python version from path (e.g., test-3.11.log)
        version_match = re.search(r'test-(\d+\.\d+)\.log', test_log_path.name)
        py_version = version_match.group(1) if version_match else "unknown"

        test_log = read_log(test_log_path)
        failures, summary, errors, status = parse_test_log(
            test_log, py_version
        )

        test_results.append({
            "version": py_version,
            "failures": failures,
            "errors": errors,
            "summary": summary,
            "status": status
        })
else:
    print("⚠️  No test logs found")
    test_results.append({
        "version": "unknown",
        "failures": [],
        "errors": [],
        "summary": "No test logs found",
        "status": "not_run"
    })

# Determine overall status (agregado a partir dos componentes)
#
# Status suportados:
# - passed: check executou e passou
# - failed: check executou e falhou
# - not_run: check não executou / artifact ausente
#
# Regras:
# - Se qualquer check falhou -> FAILURES DETECTED
# - Se todos estão not_run -> CHECKS NOT RUN
# - Se todos passaram -> ALL PASSED
# - Caso misto (pass + not_run) -> INCOMPLETE — SOME CHECKS NOT RUN


def normalize_status(status: str) -> str:
    if status in {"passed", "failed", "not_run"}:
        return status
    return "not_run"


lint_status = normalize_status(lint_status)
type_status = normalize_status(type_status)
for t in test_results:
    t["status"] = normalize_status(t.get("status", "not_run"))

component_statuses = (
    [lint_status, type_status] + [t["status"] for t in test_results]
)

any_failed = any(s == "failed" for s in component_statuses)
all_not_run = all(s == "not_run" for s in component_statuses)
all_passed = all(s == "passed" for s in component_statuses)

if any_failed:
    overall_status = "❌ FAILURES DETECTED"
elif all_not_run:
    overall_status = "⚠️ CHECKS NOT RUN"
elif all_passed:
    overall_status = "✅ ALL PASSED"
else:
    overall_status = "⚠️ INCOMPLETE — SOME CHECKS NOT RUN"

# Format status badges


def status_badge(status):
    if status == "passed":
        return "✅ PASSED"
    elif status == "failed":
        return "❌ FAILED"
    else:
        return "⚠️  NOT RUN"


lint_badge = status_badge(lint_status)
type_badge = status_badge(type_status)

# Generate report
print("📝 Generating report...")
report = f"""# CI/CD Error Summary Report

**Generated:** {TIMESTAMP}
**Commit:** `{COMMIT}`
**Branch:** `{BRANCH}`
**Run:** [#{RUN_NUMBER}]({SERVER_URL}/{REPO}/actions/runs/{RUN_ID})

---

## Overall Status

{overall_status}

---

## 1️⃣ Lint Check (Ruff + Black)

**Status:** {lint_badge}

"""

if lint_status == "failed" and lint_errors:
    report += f"**Errors Found:** {len(lint_errors)}\n\n"
    report += "<details>\n<summary>Show lint errors</summary>\n\n```\n"
    for error in lint_errors[:MAX_LINT_ERRORS_DISPLAYED]:
        report += f"{error}\n"
    if len(lint_errors) > MAX_LINT_ERRORS_DISPLAYED:
        more_lint_errors = len(lint_errors) - MAX_LINT_ERRORS_DISPLAYED
        report += f"\n... and {more_lint_errors} more errors\n"
    report += "```\n</details>\n\n"
elif lint_status == "passed":
    report += "No lint errors found. ✨\n\n"
else:
    # Evita vazar conteúdo bruto de log
    # (pode conter paths locais, tokens, etc).
    report += "Lint output omitted for safety.\n\n"
    if should_include_log_excerpt(lint_status, lint_log):
        excerpt = sanitize_log_excerpt(lint_log)
        if excerpt:
            report += (
                "<details>\n"
                "<summary>Show sanitized lint excerpt</summary>\n\n"
                "```\n"
            )
            report += f"{excerpt}\n"
            report += "```\n</details>\n\n"

report += """---

## 2️⃣ Type Check (Mypy)

**Status:** """ + type_badge + "\n\n"

if type_status == "failed" and type_errors:
    report += f"**Errors Found:** {len(type_errors)}\n\n"
    report += "<details>\n<summary>Show type errors</summary>\n\n```\n"
    for error in type_errors[:MAX_TYPE_ERRORS_DISPLAYED]:
        report += f"{error}\n"
    if len(type_errors) > MAX_TYPE_ERRORS_DISPLAYED:
        more_type_errors = len(type_errors) - MAX_TYPE_ERRORS_DISPLAYED
        report += f"\n... and {more_type_errors} more errors\n"
    report += "```\n</details>\n\n"
elif type_status == "passed":
    report += "No type errors found. ✨\n\n"
else:
    # Evita vazar conteúdo bruto de log
    # (pode conter paths locais, tokens, etc).
    report += "Type check output omitted for safety.\n\n"
    if should_include_log_excerpt(type_status, type_log):
        excerpt = sanitize_log_excerpt(type_log)
        if excerpt:
            report += (
                "<details>\n"
                "<summary>Show sanitized type-check excerpt</summary>\n\n"
                "```\n"
            )
            report += f"{excerpt}\n"
            report += "```\n</details>\n\n"

report += """---

## 3️⃣ Tests (Pytest)

"""

for i, test_result in enumerate(test_results):
    test_badge = status_badge(test_result["status"])

    report += f"""### Python {test_result['version']}

**Status:** {test_badge}

"""

    if test_result["summary"]:
        report += f"**Summary:** {test_result['summary']}\n\n"

    if test_result["failures"]:
        report += f"**Failures:** {len(test_result['failures'])}\n\n"
        report += "<details>\n<summary>Show test failures</summary>\n\n```\n"
        for failure in test_result["failures"][:MAX_TEST_FAILURES_DISPLAYED]:
            report += f"{failure}\n"
        if len(test_result["failures"]) > MAX_TEST_FAILURES_DISPLAYED:
            more_failures = (
                len(test_result["failures"]) - MAX_TEST_FAILURES_DISPLAYED
            )
            report += f"\n... and {more_failures} more failures\n"
        report += "```\n</details>\n\n"

    if test_result["errors"]:
        report += f"**Errors:** {len(test_result['errors'])}\n\n"
        report += "<details>\n<summary>Show test errors</summary>\n\n```\n"
        errors_to_show = test_result["errors"][:MAX_TEST_ERRORS_DISPLAYED]
        for error in errors_to_show:
            report += f"{error}\n"
        if len(test_result["errors"]) > MAX_TEST_ERRORS_DISPLAYED:
            more_errors = (
                len(test_result["errors"]) - MAX_TEST_ERRORS_DISPLAYED
            )
            report += f"\n... and {more_errors} more errors\n"
        report += "```\n</details>\n\n"

    if test_result["status"] == "passed":
        report += "All tests passed! ✨\n\n"

report += """---

## 📊 Summary Table

| Component | Status | Details |
|-----------|--------|---------|
"""
report += (
    "| Lint (Ruff + Black) | "
    + lint_badge
    + " | "
    + (
        "Not run"
        if lint_status == "not_run"
        else (f"{len(lint_errors)} errors" if lint_errors else "Clean")
    )
    + " |\n"
)
report += (
    "| Type Check (Mypy) | "
    + type_badge
    + " | "
    + (
        "Not run"
        if type_status == "not_run"
        else (f"{len(type_errors)} errors" if type_errors else "Clean")
    )
    + " |\n"
)

for test_result in test_results:
    test_badge = status_badge(test_result["status"])
    if test_result["status"] == "not_run":
        details = "Not run"
    elif test_result["errors"]:
        details = f"{len(test_result['errors'])} errors"
    elif test_result["failures"]:
        details = f"{len(test_result['failures'])} failures"
    else:
        details = "All passed"
    report += (
        f"| Tests (Python {test_result['version']}) | {test_badge} | "
        f"{details} |\n"
    )

report += f"| **Overall** | **{overall_status}** | - |\n"

report += """
---

## 🔧 Next Steps

"""

if lint_errors:
    report += """
### Fix Lint Errors
Run locally:
```bash
poetry run ruff check src/ tests/
poetry run black src/ tests/
```
"""

if type_errors:
    report += """
### Fix Type Errors
Run locally:
```bash
poetry run mypy src/grimperium --ignore-missing-imports
```
"""

if any(t["failures"] or t["errors"] for t in test_results):
    report += """
### Fix Test Failures
Run locally:
```bash
poetry run pytest tests/ -v
```
"""

if all_passed:
    report += """
### ✅ All Checks Passed!
Great work! All CI checks are passing.
"""

# Write report file
print(f"💾 Writing report to {REPORT_FILE}...")
REPORT_FILE.write_text(report)
print(f"✅ Report written to {REPORT_FILE}")

# Append to Job Summary (visible in GitHub UI)
print("📤 Updating GitHub Step Summary...")
try:
    GITHUB_STEP_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with open(GITHUB_STEP_SUMMARY, "a") as f:
        f.write(report)
        f.write("\n\n---\n\n")
        f.write(
            "📥 **Full report available in artifacts: "
            "`CI-Error-Summary-Report`**\n"
        )
    print("✅ Job Summary updated")
except Exception as e:
    print(f"⚠️  Could not update Job Summary: {e}")

print("\n" + "="*60)
print("REPORT GENERATION COMPLETE")
print("="*60)
