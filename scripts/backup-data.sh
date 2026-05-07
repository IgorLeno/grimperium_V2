#!/usr/bin/env bash
# backup-data.sh — Push Grimperium data to GitHub (private) and Google Drive
#
# Usage:
#   ./scripts/backup-data.sh           # both targets
#   ./scripts/backup-data.sh github    # GitHub only
#   ./scripts/backup-data.sh gdrive    # Google Drive only
#
# Requirements:
#   - git + gh CLI authenticated (for GitHub)
#   - rclone configured with remote "gdrive" (for Google Drive)

set -euo pipefail

PROJ_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="$PROJ_ROOT/data"
DATA_REPO="$HOME/Projetos/grimperium-data"
GDRIVE_REMOTE="gdrive:grimperium-backups"
TIMESTAMP="$(date +%Y-%m-%d_%H%M%S)"

FILES=(
    "thermo_pm7.csv"
    "thermo_pm7.csv.bak"
    "thermo_cbs_chon.csv"
)

target="${1:-both}"

log() { echo "[backup] $(date +%H:%M:%S) $*"; }

# ─── GitHub push ───────────────────────────────────────────────
backup_github() {
    log "Starting GitHub backup..."

    if [[ ! -d "$DATA_REPO/.git" ]]; then
        echo "ERROR: Data repo not found at $DATA_REPO" >&2
        echo "  Clone it: git clone https://github.com/IgorLeno/grimperium-data.git $DATA_REPO" >&2
        return 1
    fi

    local changed=false
    for f in "${FILES[@]}"; do
        src="$DATA_DIR/$f"
        dst="$DATA_REPO/$f"
        if [[ -f "$src" ]]; then
            if ! cmp -s "$src" "$dst" 2>/dev/null; then
                cp "$src" "$dst"
                changed=true
                log "  Updated: $f"
            else
                log "  Unchanged: $f"
            fi
        else
            log "  SKIP (not found): $f"
        fi
    done

    if $changed; then
        cd "$DATA_REPO"
        git add -A
        git commit -m "backup: $TIMESTAMP"
        git push
        log "GitHub push complete."
    else
        log "No changes to push to GitHub."
    fi
}

# ─── Google Drive sync ─────────────────────────────────────────
backup_gdrive() {
    log "Starting Google Drive backup..."

    if ! command -v rclone &>/dev/null; then
        echo "ERROR: rclone not installed. Install with: sudo apt install rclone" >&2
        return 1
    fi

    if ! rclone listremotes 2>/dev/null | grep -q "^gdrive:"; then
        echo "ERROR: rclone remote 'gdrive' not configured." >&2
        echo "  Run: rclone config  (choose Google Drive, name it 'gdrive')" >&2
        return 1
    fi

    for f in "${FILES[@]}"; do
        src="$DATA_DIR/$f"
        if [[ -f "$src" ]]; then
            rclone copy "$src" "$GDRIVE_REMOTE/" --progress
            log "  Uploaded: $f"
        else
            log "  SKIP (not found): $f"
        fi
    done

    log "Google Drive backup complete."
}

# ─── Main ──────────────────────────────────────────────────────
case "$target" in
    github)  backup_github ;;
    gdrive)  backup_gdrive ;;
    both)
        backup_github
        echo ""
        backup_gdrive
        ;;
    *)
        echo "Usage: $0 [github|gdrive|both]" >&2
        exit 1
        ;;
esac

log "Done."
