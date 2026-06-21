"""Executable module for python -m grimperium_mini."""

import sys


def main() -> int:
    if "--no-tui" in sys.argv:
        from .cli import main as cli_main

        sys.argv.remove("--no-tui")
        return cli_main()
    from .app import MiniApp

    return MiniApp().run()


if __name__ == "__main__":
    sys.exit(main())
