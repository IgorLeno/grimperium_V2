"""
Entry point for running GRIMPERIUM CLI as a module.

Usage:
    python -m grimperium.cli
"""

import sys

from grimperium.cli.app import main

if __name__ == "__main__":
    sys.exit(main())
