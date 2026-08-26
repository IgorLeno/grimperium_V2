"""Focused Semi-Imperium application package."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("grimperium")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ["__version__"]
