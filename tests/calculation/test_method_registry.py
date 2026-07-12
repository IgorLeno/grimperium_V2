"""Tests for calculation method registry and YAML definitions."""

from __future__ import annotations

import pytest

from grimperium.calculation.methods.registry import (
    CalculationMethodDefinition,
    discover_calculation_methods,
    get_calculation_method,
    list_calculation_methods,
    list_calculation_properties,
    parse_method_definition,
)


def test_parse_method_definition_rejects_unknown_feature_schema() -> None:
    payload = {
        "method_id": "bad_method",
        "version": "0.1.0",
        "display_name": "Bad Method",
        "property_id": "standard_enthalpy_of_formation",
        "property_name": "Standard enthalpy of formation",
        "conformer_selection": {"strategy": "lowest_crest_energy", "settings": {}},
        "model_requirement": {
            "model_required": True,
            "model_binding": "session_selected",
        },
        "compatibility": {
            "property": "standard_enthalpy_of_formation",
            "baseline_program": "mopac",
            "baseline_hamiltonian": "PM7",
            "feature_schema": "unknown_schema",
        },
        "xtb": {"optional": True, "enabled_by_default": True},
        "output": {"units": ["kcal/mol", "kJ/mol", "both"]},
    }

    with pytest.raises(ValueError, match="Unknown feature schema"):
        parse_method_definition(payload)


def test_parse_method_definition_rejects_absolute_paths() -> None:
    payload = {
        "method_id": "bad_method",
        "version": "0.1.0",
        "display_name": "Bad Method",
        "property_id": "standard_enthalpy_of_formation",
        "property_name": "Standard enthalpy of formation",
        "conformer_selection": {
            "strategy": "lowest_crest_energy",
            "settings": {"scratch_path": "/tmp/grimperium"},
        },
        "model_requirement": {"model_required": False, "model_binding": None},
        "compatibility": {
            "property": "standard_enthalpy_of_formation",
            "baseline_program": "mopac",
            "baseline_hamiltonian": None,
            "feature_schema": None,
        },
        "xtb": {"optional": True, "enabled_by_default": True},
        "output": {"units": ["kcal/mol", "kJ/mol", "both"]},
    }

    with pytest.raises(ValueError, match="absolute path"):
        parse_method_definition(payload)


def test_list_standard_enthalpy_methods_returns_method_a_and_b() -> None:
    methods = list_calculation_methods("standard_enthalpy_of_formation")

    assert [method.method_id for method in methods] == [
        "crest_pm7",
        "pm7_delta_learning",
        "semiempirical_am1_pm3_pm7",
    ]
    assert all(isinstance(method, CalculationMethodDefinition) for method in methods)


def test_list_calculation_properties_exposes_standard_enthalpy() -> None:
    properties = list_calculation_properties()

    assert [prop.property_id for prop in properties] == [
        "standard_enthalpy_of_formation"
    ]


def test_method_a_definition_matches_pr4_contract() -> None:
    method = get_calculation_method(
        "semiempirical_am1_pm3_pm7",
        property_id="standard_enthalpy_of_formation",
    )

    assert method.property_id == "standard_enthalpy_of_formation"
    assert method.conformer_selection.strategy == "lowest_crest_energy"
    assert method.conformer_selection.settings == {"max_selected": 1}
    assert method.model_requirement.model_required is False
    assert method.model_requirement.model_binding is None
    assert method.xtb.optional is True
    assert method.xtb.enabled_by_default is True


def test_method_b_definition_matches_pr4_contract() -> None:
    method = get_calculation_method(
        "pm7_delta_learning",
        property_id="standard_enthalpy_of_formation",
    )

    assert method.conformer_selection.strategy == "lowest_pm7_hof_within_crest_subset"
    assert method.conformer_selection.settings == {"crest_subset_size": 3}
    assert method.model_requirement.model_required is True
    assert method.model_requirement.model_binding == "session_selected"
    assert method.compatibility.baseline_program == "mopac"
    assert method.compatibility.baseline_hamiltonian == "PM7"
    assert method.compatibility.feature_schema == "grimperium_dhf_v1"


def test_discovery_rejects_invalid_yaml(tmp_path, monkeypatch) -> None:
    package_dir = tmp_path / "badpkg"
    package_dir.mkdir()
    (package_dir / "__init__.py").write_text("", encoding="utf-8")
    (package_dir / "bad.yaml").write_text(
        "method_id: bad\nversion: 1.0\n",
        encoding="utf-8",
    )
    monkeypatch.syspath_prepend(str(tmp_path))

    with pytest.raises(ValueError, match="Invalid method definition"):
        discover_calculation_methods({"standard_enthalpy_of_formation": "badpkg"})


def test_discovery_rejects_duplicate_method_ids(tmp_path, monkeypatch) -> None:
    package_dir = tmp_path / "duppkg"
    package_dir.mkdir()
    (package_dir / "__init__.py").write_text("", encoding="utf-8")
    yaml_text = """
method_id: duplicate
version: 1.0
display_name: Duplicate
property_id: standard_enthalpy_of_formation
property_name: Standard enthalpy of formation
conformer_selection:
  strategy: lowest
  settings:
    max_selected: 1
model_requirement:
  model_required: false
  model_binding:
compatibility:
  property: standard_enthalpy_of_formation
  baseline_program: mopac
  baseline_hamiltonian:
  feature_schema:
xtb:
  optional: true
  enabled_by_default: false
output:
  units:
    - kcal/mol
"""
    (package_dir / "a.yaml").write_text(yaml_text, encoding="utf-8")
    (package_dir / "b.yaml").write_text(yaml_text, encoding="utf-8")
    monkeypatch.syspath_prepend(str(tmp_path))

    with pytest.raises(ValueError, match="Duplicate calculation method_id"):
        discover_calculation_methods({"standard_enthalpy_of_formation": "duppkg"})
