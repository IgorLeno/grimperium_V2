"""Styles and constants for grimperium_mini TUI."""

from rich.console import Console

MINI_BANNER = r"""
   ____ ____  ___ __  __ ____  _____ ____  ___ _   _ __  __   __  __ ___ _   _ ___
  / ___|  _ \|_ _|  \/  |  _ \| ____|  _ \|_ _| | | |  \/  | |  \/  |_ _| \ | |_ _|
 | |  _| |_) || || |\/| | |_) |  _| | |_) || || | | | |\/| | | |\/| || ||  \| || |
 | |_| |  _ < | || |  | |  __/| |___|  _ < | || |_| | |  | | | |  | || || |\  || |
  \____|_| \_\___|_|  |_|_|   |_____|_| \_\___|\___/|_|  |_| |_|  |_|___|_| \_|___|
"""

MINI_SUBTITLE = "[bold cyan]TCC Pipeline · CREST + MOPAC · AM1 / PM3 / PM7[/bold cyan]"
MINI_VERSION = "[dim]v0.1.0[/dim]"

COLORS = {
    "primary": "cyan",
    "success": "green",
    "warning": "yellow",
    "error": "red",
    "muted": "dim",
    "border": "bright_black",
}

console = Console()
