#!/usr/bin/env bash
set -uo pipefail

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

    # Fallback — download crest binary from continuous release (tag: latest)
    # NOTE: The stable release v3.0.2 does NOT have pre-built binaries.
    # The continuous release (tag "latest") has: crest-gnu-12-ubuntu-latest.tar.xz
    CREST_URL="https://github.com/crest-lab/crest/releases/download/latest/crest-gnu-12-ubuntu-latest.tar.xz"
    echo ">>> Downloading CREST (continuous release)..."
    echo "    $CREST_URL"
    mkdir -p /tmp/crest_extracted
    if wget -q --show-progress "$CREST_URL" -O /tmp/crest.tar.xz; then
        tar -xf /tmp/crest.tar.xz -C /tmp/crest_extracted/
        CREST_BIN=$(find /tmp/crest_extracted -name "crest" -type f | head -1)
        if [ -n "$CREST_BIN" ]; then
            cp "$CREST_BIN" "$TOOLS_DIR/crest"
            chmod +x "$TOOLS_DIR/crest"
            echo ">>> CREST installed to $TOOLS_DIR/crest"
        else
            echo "!!! crest binary not found inside archive."
        fi
    else
        echo "!!! CREST download failed."
        echo "    Install manually: https://crest-lab.github.io/crest-docs/page/installation"
    fi
}

# ── Install MOPAC ────────────────────────────────────────────────────────────
install_mopac() {
    if command -v dnf &>/dev/null; then
        sudo dnf install mopac -y || true
    fi

    if ! command -v mopac &>/dev/null; then
        MOPAC_URL=$(curl -s https://api.github.com/repos/openmopac/mopac/releases/latest \
            | grep "browser_download_url" \
            | grep -i "linux" \
            | grep -v '\.sha256\|\.asc\|\.sig' \
            | head -1 \
            | cut -d '"' -f 4)

        if [ -z "$MOPAC_URL" ]; then
            echo "!!! Could not resolve MOPAC URL. Install manually: https://github.com/openmopac/mopac/releases"
        else
            wget -q --show-progress "$MOPAC_URL" -O /tmp/mopac_pkg \
                && tar -xf /tmp/mopac_pkg -C "$TOOLS_DIR" 2>/dev/null \
                || cp /tmp/mopac_pkg "$TOOLS_DIR/mopac" \
                && chmod +x "$TOOLS_DIR"/mopac* 2>/dev/null \
                || echo "!!! MOPAC install failed — install manually."
        fi
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
