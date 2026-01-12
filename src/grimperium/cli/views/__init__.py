"""
View modules for GRIMPERIUM CLI.

Each view handles a specific section of the application.
"""

from grimperium.cli.views.about_view import AboutView
from grimperium.cli.views.base_view import BaseView
from grimperium.cli.views.calc_view import CalcView
from grimperium.cli.views.databases_view import DatabasesView
from grimperium.cli.views.models_view import ModelsView
from grimperium.cli.views.results_view import ResultsView
from grimperium.cli.views.settings_view import SettingsView

__all__ = [
    "BaseView",
    "CalcView",
    "DatabasesView",
    "ModelsView",
    "ResultsView",
    "SettingsView",
    "AboutView",
]
