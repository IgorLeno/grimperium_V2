#!/usr/bin/env bash
set -euo pipefail

TOOLS_DIR="$HOME/.local/bin"
mkdir -p "$TOOLS_DIR"

OS=$(uname -s)    # Linux | Darwin
ARCH=$(uname -m)  # x86_64 | aarch64 | arm64

# ── Add ~/.local/bin to PATH in shell rc files ───────────────────────────────
for rc in ~/.bashrc ~/.zshrc ~/.profile; do
    if [ -f "$rc" ] && ! grep -q 'HOME/.local/bin' "$rc"; then
        echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$rc"
    fi
done
export PATH="$HOME/.local/bin:$PATH"

# ── Install CREST + xTB ──────────────────────────────────────────────────────
install_xtb_crest() {
    if command -v conda &>/dev/null; then
        conda install -c conda-forge crest xtb -y
        return
    fi

    if command -v mamba &>/dev/null; then
        mamba install -c conda-forge crest xtb -y
        return
    fi

    if command -v dnf &>/dev/null; then
        sudo dnf install xtb -y || true
    fi

    # Fallback — download crest binary directly from GitHub Releases
    URL="https://github.com/crest-lab/crest/releases/latest/download/crest-gnu-12-amd64-static.tar.xz"
    wget -q --show-progress "$URL" -O /tmp/crest.tar.xz
    tar -xf /tmp/crest.tar.xz -C "$TOOLS_DIR"
    chmod +x "$TOOLS_DIR/crest"
}

# ── Install MOPAC ────────────────────────────────────────────────────────────
install_mopac() {
    if command -v dnf &>/dev/null; then
        sudo dnf install mopac -y || true
    fi

    if ! command -v mopac &>/dev/null; then
        URL="https://github.com/openmopac/mopac/releases/latest/download/mopac-linux.tar.gz"
        wget -q --show-progress "$URL" -O /tmp/mopac.tar.gz
        tar -xf /tmp/mopac.tar.gz -C "$TOOLS_DIR"
        chmod +x "$TOOLS_DIR/mopac"
    fi
}

install_xtb_crest
install_mopac

echo ""
echo "Installation complete. Installed tools:"
command -v crest && crest --version 2>/dev/null | head -1 || echo "  crest: NOT FOUND"
command -v xtb   && xtb --version   2>/dev/null | head -1 || echo "  xtb:   NOT FOUND"
command -v mopac && mopac -v        2>/dev/null | head -1 || echo "  mopac: NOT FOUND"
echo ""
echo "If tools are still not found, restart the terminal and run grimperium again."
