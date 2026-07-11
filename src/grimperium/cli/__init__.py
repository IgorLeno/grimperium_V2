"""
GRIMPERIUM CLI - Interactive Terminal Interface.

This package provides a professional CLI application for
molecular property prediction using Rich and questionary.

Usage:
    python -m grimperium.cli
    # or
    grimperium-cli

Modules:
    app: Main application class
    controller: Navigation and state management
    menu: Interactive menu system
    styles: Colors, themes, and ASCII art
    viewmodels: UI-facing result structures
    views: View modules for each section
"""

from grimperium.cli.app import GrimperiumCLI, main

__all__ = ["GrimperiumCLI", "main"]
