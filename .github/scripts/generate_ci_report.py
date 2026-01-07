#!/usr/bin/env python3
"""
Generate comprehensive CI error report from captured logs
"""
import os
import re
from pathlib import Path
from datetime import datetime

# Define paths
WORKSPACE = Path(os.environ.get("GITHUB_WORKSPACE", "."))
LOGS_DIR = WORKSPACE / "logs"
LINT_LOG = LOGS_DIR / "lint-log" / "lint.log"
TYPE_LOG = LOGS_DIR / "type-log" / "type.log"
TEST_LOGS = list(LOGS_DIR.glob("test-log-*/test-*.log"))
REPORT_FILE = WORKSPACE / "CI_ERROR_SUMMARY.md"
GITHUB_STEP_SUMMARY = Path(os.environ.get("GITHUB_STEP_SUMMARY", "/tmp/summary.md"))

# Collect metadata
COMMIT = os.environ.get("GITHUB_SHA", "unknown")[:12]
BRANCH = os.environ.get("GITHUB_REF_NAME", "unknown")
RUN_NUMBER = os.environ.get("GITHUB_RUN_NUMBER", "unknown")
RUN_ID = os.environ.get("GITHUB_RUN_ID", "unknown")
REPO = os.environ.get("GITHUB_REPOSITORY", "unknown")
SERVER_URL = os.environ.get("GITHUB_SERVER_URL", "https://github.com")
TIMESTAMP = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")


def read_log(path):
    """Read log file safely"""
    if path.exists():
        try:
            return path.read_text()
        except Exception as e:
            return f"Error reading log: {e}"
    return "Log file not found"


def parse_lint_log(text):
    """Parse ruff/black output"""
    if "Log file not found" in text or "Error reading" in text:
        return [], text

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
    status = "not found" if "not found" in text else ("passed" if not all_errors else "failed")

    return all_errors, status


def parse_type_log(text):
    """Parse mypy output"""
    if "Log file not found" in text or "Error reading" in text:
        return [], text

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
        status = "not found"

    return errors, status


def parse_test_log(text, python_version=""):
    """Parse pytest output"""
    if "Log file not found" in text or "Error reading" in text:
        return [], None, [], text

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
        failures, summary, errors, status = parse_test_log(test_log, py_version)

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
        "status": "not found"
    })

# Determine overall status
all_passed = (
    lint_status in ["passed", "not found"] and
    type_status in ["passed", "not found"] and
    all(t["status"] in ["passed", "not found"] for t in test_results)
)

overall_status = "✅ ALL PASSED" if all_passed else "❌ FAILURES DETECTED"

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
    for error in lint_errors[:30]:  # Limit to first 30
        report += f"{error}\n"
    if len(lint_errors) > 30:
        report += f"\n... and {len(lint_errors) - 30} more errors\n"
    report += "```\n</details>\n\n"
elif lint_status == "passed":
    report += "No lint errors found. ✨\n\n"
else:
    report += f"{lint_log[:200]}\n\n"

report += """---

## 2️⃣ Type Check (Mypy)

**Status:** """ + type_badge + "\n\n"

if type_status == "failed" and type_errors:
    report += f"**Errors Found:** {len(type_errors)}\n\n"
    report += "<details>\n<summary>Show type errors</summary>\n\n```\n"
    for error in type_errors[:30]:
        report += f"{error}\n"
    if len(type_errors) > 30:
        report += f"\n... and {len(type_errors) - 30} more errors\n"
    report += "```\n</details>\n\n"
elif type_status == "passed":
    report += "No type errors found. ✨\n\n"
else:
    report += f"{type_log[:200]}\n\n"

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
        for failure in test_result["failures"][:20]:
            report += f"{failure}\n"
        if len(test_result["failures"]) > 20:
            report += f"\n... and {len(test_result['failures']) - 20} more failures\n"
        report += "```\n</details>\n\n"

    if test_result["errors"]:
        report += f"**Errors:** {len(test_result['errors'])}\n\n"
        report += "<details>\n<summary>Show test errors</summary>\n\n```\n"
        for error in test_result["errors"][:20]:
            report += f"{error}\n"
        if len(test_result["errors"]) > 20:
            report += f"\n... and {len(test_result['errors']) - 20} more errors\n"
        report += "```\n</details>\n\n"

    if test_result["status"] == "passed":
        report += "All tests passed! ✨\n\n"

report += """---

## 📊 Summary Table

| Component | Status | Details |
|-----------|--------|---------|
| Lint (Ruff + Black) | """ + lint_badge + " | " + (f"{len(lint_errors)} errors" if lint_errors else "Clean") + """ |
| Type Check (Mypy) | """ + type_badge + " | " + (f"{len(type_errors)} errors" if type_errors else "Clean") + """ |
"""

for test_result in test_results:
    test_badge = status_badge(test_result["status"])
    details = f"{len(test_result['failures'])} failures" if test_result['failures'] else "All passed"
    report += f"| Tests (Python {test_result['version']}) | {test_badge} | {details} |\n"

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
print(f"📤 Updating GitHub Step Summary...")
try:
    GITHUB_STEP_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with open(GITHUB_STEP_SUMMARY, "a") as f:
        f.write(report)
        f.write("\n\n---\n\n")
        f.write("📥 **Full report available in artifacts: `CI-Error-Summary-Report`**\n")
    print("✅ Job Summary updated")
except Exception as e:
    print(f"⚠️  Could not update Job Summary: {e}")

print("\n" + "="*60)
print("REPORT GENERATION COMPLETE")
print("="*60)
