"""
Styles and theming for the GRIMPERIUM CLI.

This module defines colors, themes, and ASCII art for the CLI interface.
"""

from rich.style import Style
from rich.theme import Theme

# Color palette
COLORS = {
    "calc": "#00D9FF",  # Cyan - predictions
    "databases": "#0080FF",  # Blue - data management
    "models": "#FF00FF",  # Magenta - ML models
    "results": "#00FF00",  # Green - analytics
    "settings": "#FFFF00",  # Yellow - configuration
    "about": "#F5F5F5",  # White - info
    "error": "#FF4444",  # Red
    "success": "#44FF44",  # Green
    "warning": "#FFAA00",  # Orange
    "in_dev": "#888888",  # Gray - in development
    "muted": "#666666",  # Dark gray
    "highlight": "#FFFFFF",  # White
    "border": "#444444",  # Border color
}

# Rich styles
STYLES = {
    "calc": Style(color=COLORS["calc"], bold=True),
    "databases": Style(color=COLORS["databases"], bold=True),
    "models": Style(color=COLORS["models"], bold=True),
    "results": Style(color=COLORS["results"], bold=True),
    "settings": Style(color=COLORS["settings"], bold=True),
    "about": Style(color=COLORS["about"], bold=True),
    "error": Style(color=COLORS["error"], bold=True),
    "success": Style(color=COLORS["success"], bold=True),
    "warning": Style(color=COLORS["warning"], bold=True),
    "in_dev": Style(color=COLORS["in_dev"], italic=True),
    "muted": Style(color=COLORS["muted"]),
    "title": Style(color=COLORS["highlight"], bold=True),
    "subtitle": Style(color=COLORS["muted"], italic=True),
}

# Rich theme for console
CLI_THEME = Theme(
    {
        "info": COLORS["about"],
        "warning": COLORS["warning"],
        "error": COLORS["error"],
        "success": COLORS["success"],
        "calc": COLORS["calc"],
        "databases": COLORS["databases"],
        "models": COLORS["models"],
        "results": COLORS["results"],
        "settings": COLORS["settings"],
        "about": COLORS["about"],
        "in_dev": COLORS["in_dev"],
        "muted": COLORS["muted"],
    }
)

# ASCII art banner
BANNER = r"""
[bold cyan]
 ██████╗ ██████╗ ██╗███╗   ███╗██████╗ ███████╗██████╗ ██╗██╗   ██╗███╗   ███╗
██╔════╝ ██╔══██╗██║████╗ ████║██╔══██╗██╔════╝██╔══██╗██║██║   ██║████╗ ████║
██║  ███╗██████╔╝██║██╔████╔██║██████╔╝█████╗  ██████╔╝██║██║   ██║██╔████╔██║
██║   ██║██╔══██╗██║██║╚██╔╝██║██╔═══╝ ██╔══╝  ██╔══██╗██║██║   ██║██║╚██╔╝██║
╚██████╔╝██║  ██║██║██║ ╚═╝ ██║██║     ███████╗██║  ██║██║╚██████╔╝██║ ╚═╝ ██║
 ╚═════╝ ╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝ ╚═════╝ ╚═╝     ╚═╝
[/bold cyan]
"""

# Compact banner for smaller terminals
BANNER_COMPACT = r"""
[bold cyan]
  _____ _____ _____ ___ _ _____ _____ _____ _____ _____ _____
 |   __| __  |     |   | |  _  |   __| __  |     |  |  |     |
 |  |  |    -|-   -| | | |   __|   __|    -|-   -|  |  | | | |
 |_____|__|__|_____|_|___|__|  |_____|__|__|_____|_____|_|_|_|
[/bold cyan]
"""

# Welcome screen subtitle
SUBTITLE = "[muted]ML-Enhanced PM7 Molecular Property Predictor[/muted]"
VERSION_LINE = "[muted]Delta-Learning Framework v1.0.0-beta[/muted]"

# Menu icons
ICONS = {
    "calc": "🧪",
    "databases": "📊",
    "models": "🤖",
    "results": "📈",
    "settings": "⚙️",
    "about": "ℹ️",
    "back": "◀",
    "exit": "✖",
    "success": "✓",
    "error": "✗",
    "warning": "⚠",
    "in_dev": "🚧",
    "arrow": "❯",
}

# Status indicators
STATUS_READY = f"[success]{ICONS['success']} Ready[/success]"
STATUS_IN_DEV = f"[in_dev]{ICONS['in_dev']} In Development[/in_dev]"
STATUS_ERROR = f"[error]{ICONS['error']} Error[/error]"

# Footer text
FOOTER_NAVIGATION = "[muted]↑↓ Navigate  ↵ Select  Ctrl+C Exit[/muted]"
FOOTER_INPUT = "[muted]Type your input and press Enter  Ctrl+C Cancel[/muted]"
