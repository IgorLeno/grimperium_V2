#!/usr/bin/env bash
# setup_worker.sh — One-shot Grimperium worker setup for Linux / WSL2
#
# Usage:
#   bash setup_worker.sh [REPO_URL]
#
# Default REPO_URL: https://github.com/IgorLeno/grimperium_V2.git
#
# What this script does:
#   1. Verify platform (Linux native or WSL2)
#   2. Install system build dependencies via apt
#   3. Install Poetry (if absent)
#   4. Ensure $HOME/.local/bin is in PATH
#   5. Clone the repository
#   6. poetry install --extras "worker"
#   7. Install CREST + MOPAC (via scripts/install_tools.sh)
#   8. Verify the worker can import correctly
#   9. Print a summary of what succeeded and what failed

set -uo pipefail

REPO_URL="${1:-https://github.com/IgorLeno/grimperium_V2.git}"
REPO_DIR="$HOME/grimperium_V2"

# ── Colour helpers ────────────────────────────────────────────────────────────
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'
ok()   { echo -e "${GREEN}✓ $*${NC}"; }
fail() { echo -e "${RED}✗ $*${NC}"; }
warn() { echo -e "${YELLOW}! $*${NC}"; }

# Tracking for final summary
declare -a PASSED=()
declare -a FAILED=()

# ── Step 1: Platform check ────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Grimperium Worker Setup"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [[ "$(uname -s)" != "Linux" ]]; then
    fail "This script only runs on Linux (native or WSL2)."
    echo "  For Windows, open Ubuntu in WSL2 and run this script from there."
    exit 1
fi

if uname -r | grep -qi microsoft; then
    ok "Platform: WSL2 detected"
else
    ok "Platform: Linux native"
fi

# ── Step 2: System build dependencies (apt only) ─────────────────────────────
echo ""
echo "[ 2/9 ] Installing system build dependencies..."
if command -v apt-get &>/dev/null; then
    sudo apt-get update -qq && \
    sudo apt-get install -y -qq \
        git curl wget build-essential \
        python3 python3-pip python3-venv \
        2>/dev/null
    ok "System dependencies installed"
    PASSED+=("System deps (apt)")
else
    warn "apt-get not found — skipping system dep install."
    warn "Make sure git, curl, wget, python3 are available."
    FAILED+=("System deps (no apt)")
fi

# ── Step 3: Poetry ───────────────────────────────────────────────────────────
echo ""
echo "[ 3/9 ] Checking Poetry..."
if ! command -v poetry &>/dev/null; then
    echo "  Installing Poetry via official installer..."
    curl -sSL https://install.python-poetry.org | python3 -
    ok "Poetry installed"
    PASSED+=("Poetry install")
else
    ok "Poetry already installed ($(poetry --version 2>/dev/null | head -1))"
    PASSED+=("Poetry (pre-existing)")
fi

# ── Step 4: PATH ─────────────────────────────────────────────────────────────
echo ""
echo "[ 4/9 ] Ensuring \$HOME/.local/bin is in PATH..."
if ! grep -q 'HOME/.local/bin' ~/.bashrc 2>/dev/null; then
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
    ok "Added \$HOME/.local/bin to ~/.bashrc"
else
    ok "\$HOME/.local/bin already in ~/.bashrc"
fi
export PATH="$HOME/.local/bin:$PATH"

# ── Step 5: Clone repository ─────────────────────────────────────────────────
echo ""
echo "[ 5/9 ] Cloning repository..."
if [[ -d "$REPO_DIR/.git" ]]; then
    ok "Repository already exists at $REPO_DIR — skipping clone"
    PASSED+=("Repo (pre-existing)")
else
    if git clone "$REPO_URL" "$REPO_DIR"; then
        ok "Cloned to $REPO_DIR"
        PASSED+=("Git clone")
    else
        fail "git clone failed. Check URL and network access."
        FAILED+=("Git clone")
        exit 1
    fi
fi

# ── Step 6: poetry install ───────────────────────────────────────────────────
echo ""
echo "[ 6/9 ] Running poetry install --extras worker ..."
cd "$REPO_DIR"
if poetry install --extras "worker" --no-interaction; then
    ok "Python dependencies installed"
    PASSED+=("poetry install --extras worker")
else
    fail "poetry install failed."
    FAILED+=("poetry install")
    exit 1
fi

# ── Step 7: CREST + MOPAC ────────────────────────────────────────────────────
echo ""
echo "[ 7/9 ] Installing CREST and MOPAC..."
if bash scripts/install_tools.sh; then
    ok "install_tools.sh completed"
    PASSED+=("install_tools.sh")
else
    warn "install_tools.sh exited non-zero — check output above"
    FAILED+=("install_tools.sh")
fi

# ── Step 7b: Configure LD_LIBRARY_PATH for current session ───────────────────
if [ -d "$HOME/.local/lib" ]; then
    export LD_LIBRARY_PATH="$HOME/.local/lib:${LD_LIBRARY_PATH:-}"
    ok "LD_LIBRARY_PATH configured for current session"
fi

# ── Step 8: Verification ─────────────────────────────────────────────────────
echo ""
echo "[ 8/9 ] Verifying installation..."

if command -v crest &>/dev/null; then
    ok "crest: $(crest --version 2>/dev/null | head -1 || echo 'found')"
    PASSED+=("crest binary")
else
    fail "crest not found in PATH"
    FAILED+=("crest binary")
fi

if command -v mopac &>/dev/null && mopac --version 2>&1 | grep -qi "^MOPAC version"; then
    ok "mopac: $(mopac --version 2>&1 | head -1)"
    PASSED+=("mopac binary")
else
    fail "mopac not found or Qt installer detected in PATH"
    FAILED+=("mopac binary")
fi

if poetry run python -c "from grimperium.worker import WorkerRunner" 2>/dev/null; then
    ok "grimperium.worker.WorkerRunner imports successfully"
    PASSED+=("Worker Python import")
else
    fail "grimperium.worker.WorkerRunner import failed"
    FAILED+=("Worker Python import")
fi

# ── Step 9: Summary ──────────────────────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Setup Summary"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

for item in "${PASSED[@]:-}"; do
    [[ -n "$item" ]] && echo -e "  ${GREEN}✓${NC} $item"
done

for item in "${FAILED[@]:-}"; do
    [[ -n "$item" ]] && echo -e "  ${RED}✗${NC} $item"
done

echo ""
if [[ ${#FAILED[@]} -eq 0 ]]; then
    echo -e "${GREEN}All steps completed successfully.${NC}"
else
    echo -e "${YELLOW}Setup completed with ${#FAILED[@]} issue(s). See above.${NC}"
fi

echo ""
echo "Next steps:"
echo "  1. Copy your .env config file into: $REPO_DIR/.env"
echo "     (Ask the coordinator for the SERVER_URL, WORKER_ID, and API_TOKEN values)"
echo ""
echo "  2. Launch the worker:"
echo "     cd $REPO_DIR"
echo "     poetry run python -m grimperium.worker.runner"
echo ""
echo "  3. If tools are not found after reboot, run:"
echo "     source ~/.bashrc"
echo ""
