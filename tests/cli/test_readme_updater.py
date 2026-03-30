"""Tests for the README updater CLI workflow."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from grimperium.cli.readme_updater import ReadmeUpdater

README_FIXTURE = """# Grimperium

Introdução.

## Resultados Atuais

Conteúdo antigo que será substituído.

---

## Outra Seção

Texto que deve permanecer idêntico.
"""


def test_collect_pm7_stats_ok(tmp_path: Path) -> None:
    """collect_pm7_stats returns filtered absolute-error metrics."""
    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        {
            "status": ["OK", "OK", "FAILED", "FAILED", "OK"],
            "H298_pm7": [12.0, 18.0, 30.0, 40.0, -100.0],
            "H298_cbs": [10.0, 20.0, 30.0, 40.0, 100.0],
            "cbs_quality_flag": ["OK", "OK", "OK", "OK", "SUSPECT"],
        }
    ).to_csv(csv_path, index=False)

    updater = ReadmeUpdater(pm7_csv_path=csv_path)
    stats = updater.collect_pm7_stats()

    assert stats["total_ok"] == 3
    assert stats["total_rows"] == 5
    assert stats["r2"] == pytest.approx(0.84)
    assert stats["mare"] == pytest.approx(2.0)
    assert stats["p50"] == pytest.approx(2.0)
    assert stats["p90"] == pytest.approx(2.0)
    assert stats["bias"] == pytest.approx(0.0)


def test_collect_pm7_stats_file_not_found(tmp_path: Path) -> None:
    """Missing CSV returns zero counts and null metrics."""
    updater = ReadmeUpdater(pm7_csv_path=tmp_path / "missing.csv")

    stats = updater.collect_pm7_stats()

    assert stats["total_ok"] == 0
    assert stats["total_rows"] == 0
    assert stats["r2"] is None


def test_collect_model_stats_no_model(tmp_path: Path) -> None:
    """Empty model directory returns model_found=False."""
    models_dir = tmp_path / "models"
    models_dir.mkdir()

    updater = ReadmeUpdater(model_bundle_path=models_dir)
    stats = updater.collect_model_stats()

    assert stats["model_found"] is False
    assert stats["r2"] is None


def test_update_readme_dry_run(tmp_path: Path) -> None:
    """Dry run returns updated content without writing the README file."""
    readme_path = tmp_path / "README.md"
    readme_path.write_text(README_FIXTURE, encoding="utf-8")

    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        {
            "status": ["OK", "OK"],
            "H298_pm7": [10.0, 20.0],
            "H298_cbs": [10.0, 20.0],
            "cbs_quality_flag": ["OK", "OK"],
        }
    ).to_csv(csv_path, index=False)

    updater = ReadmeUpdater(
        readme_path=readme_path,
        pm7_csv_path=csv_path,
        model_bundle_path=tmp_path / "models",
    )

    success, new_content = updater.update_readme(dry_run=True)

    assert success is True
    assert "| Moléculas calculadas (status OK) | 2 |" in new_content
    assert "| Correlação PM7 vs CBS-QB3 (R²) | 1.0000 |" in new_content
    assert "Nenhum modelo ML treinado ainda" in new_content
    assert readme_path.read_text(encoding="utf-8") == README_FIXTURE


def test_update_readme_preserves_other_sections(tmp_path: Path) -> None:
    """Only the dynamic README section is replaced."""
    readme_path = tmp_path / "README.md"
    readme_path.write_text(README_FIXTURE, encoding="utf-8")

    csv_path = tmp_path / "thermo_pm7.csv"
    pd.DataFrame(
        {
            "status": ["OK"],
            "H298_pm7": [5.0],
            "H298_cbs": [4.0],
            "cbs_quality_flag": ["OK"],
        }
    ).to_csv(csv_path, index=False)

    updater = ReadmeUpdater(
        readme_path=readme_path,
        pm7_csv_path=csv_path,
        model_bundle_path=tmp_path / "models",
    )

    success, new_content = updater.update_readme(dry_run=True)

    assert success is True
    assert new_content.startswith("# Grimperium\n\nIntrodução.\n\n")
    assert new_content.endswith(
        "\n---\n\n## Outra Seção\n\nTexto que deve permanecer idêntico.\n"
    )
