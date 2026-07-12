"""Static registry for declarative calculation method definitions."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from importlib import resources
from pathlib import PurePosixPath, PureWindowsPath
from typing import Any

from grimperium.calculation.methods.feature_schema import get_feature_schema

STANDARD_ENTHALPY_PACKAGE = (
    "grimperium.calculation.methods.standard_enthalpy_of_formation"
)
METHOD_PROPERTY_PACKAGES = {
    "standard_enthalpy_of_formation": STANDARD_ENTHALPY_PACKAGE,
}


@dataclass(frozen=True)
class CalculationPropertyDefinition:
    """Property catalogue entry visible in the CLI."""

    property_id: str
    display_name: str
    package: str


@dataclass(frozen=True)
class ConformerSelectionDefinition:
    """Conformer selection strategy declared by a calculation method."""

    strategy: str
    settings: dict[str, Any]


@dataclass(frozen=True)
class ModelRequirementDefinition:
    """Model binding requirements for a calculation method."""

    model_required: bool
    model_binding: str | None


@dataclass(frozen=True)
class CompatibilityDefinition:
    """Compatibility metadata used before applying a model-backed method."""

    property: str
    baseline_program: str | None
    baseline_hamiltonian: str | None
    feature_schema: str | None


@dataclass(frozen=True)
class XtbDefinition:
    """xTB behavior declared by a calculation method."""

    optional: bool
    enabled_by_default: bool


@dataclass(frozen=True)
class CalculationMethodDefinition:
    """Validated calculation method definition loaded from package YAML."""

    method_id: str
    version: str
    display_name: str
    property_id: str
    property_name: str
    conformer_selection: ConformerSelectionDefinition
    model_requirement: ModelRequirementDefinition
    compatibility: CompatibilityDefinition
    xtb: XtbDefinition
    output: dict[str, Any]


PROPERTY_CATALOG = (
    CalculationPropertyDefinition(
        property_id="standard_enthalpy_of_formation",
        display_name="Standard enthalpy of formation",
        package=STANDARD_ENTHALPY_PACKAGE,
    ),
)


def _require_mapping(value: Any, field_name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Method definition field {field_name!r} must be a mapping")
    return value


def _require_str(payload: Mapping[str, Any], field_name: str) -> str:
    value = payload.get(field_name)
    if not isinstance(value, str) or not value:
        raise ValueError(f"Method definition field {field_name!r} must be a string")
    return value


def _optional_str(payload: Mapping[str, Any], field_name: str) -> str | None:
    value = payload.get(field_name)
    if value is None:
        return None
    if not isinstance(value, str) or not value:
        raise ValueError(f"Method definition field {field_name!r} must be a string")
    return value


def _require_bool(payload: Mapping[str, Any], field_name: str) -> bool:
    value = payload.get(field_name)
    if not isinstance(value, bool):
        raise ValueError(f"Method definition field {field_name!r} must be a boolean")
    return value


def _is_absolute_path(value: str) -> bool:
    return PurePosixPath(value).is_absolute() or PureWindowsPath(value).is_absolute()


def _reject_absolute_paths(value: Any, path: str = "root") -> None:
    if isinstance(value, str):
        if _is_absolute_path(value):
            raise ValueError(f"Method definition contains absolute path at {path}")
    elif isinstance(value, Mapping):
        for key, nested_value in value.items():
            _reject_absolute_paths(nested_value, f"{path}.{key}")
    elif isinstance(value, list):
        for index, nested_value in enumerate(value):
            _reject_absolute_paths(nested_value, f"{path}[{index}]")


def parse_method_definition(payload: Mapping[str, Any]) -> CalculationMethodDefinition:
    """Validate and convert raw YAML payload into a method definition."""
    _reject_absolute_paths(payload)

    conformer_payload = _require_mapping(
        payload.get("conformer_selection"), "conformer_selection"
    )
    model_payload = _require_mapping(
        payload.get("model_requirement"), "model_requirement"
    )
    compatibility_payload = _require_mapping(
        payload.get("compatibility"), "compatibility"
    )
    xtb_payload = _require_mapping(payload.get("xtb"), "xtb")
    output_payload = _require_mapping(payload.get("output"), "output")

    feature_schema = _optional_str(compatibility_payload, "feature_schema")
    if feature_schema is not None:
        get_feature_schema(feature_schema)

    return CalculationMethodDefinition(
        method_id=_require_str(payload, "method_id"),
        version=_require_str(payload, "version"),
        display_name=_require_str(payload, "display_name"),
        property_id=_require_str(payload, "property_id"),
        property_name=_require_str(payload, "property_name"),
        conformer_selection=ConformerSelectionDefinition(
            strategy=_require_str(conformer_payload, "strategy"),
            settings=dict(
                _require_mapping(conformer_payload.get("settings"), "settings")
            ),
        ),
        model_requirement=ModelRequirementDefinition(
            model_required=_require_bool(model_payload, "model_required"),
            model_binding=_optional_str(model_payload, "model_binding"),
        ),
        compatibility=CompatibilityDefinition(
            property=_require_str(compatibility_payload, "property"),
            baseline_program=_optional_str(compatibility_payload, "baseline_program"),
            baseline_hamiltonian=_optional_str(
                compatibility_payload,
                "baseline_hamiltonian",
            ),
            feature_schema=feature_schema,
        ),
        xtb=XtbDefinition(
            optional=_require_bool(xtb_payload, "optional"),
            enabled_by_default=_require_bool(xtb_payload, "enabled_by_default"),
        ),
        output=dict(output_payload),
    )


def _parse_scalar(value: str) -> str | bool | int | None:
    if value == "":
        return None
    if value == "true":
        return True
    if value == "false":
        return False
    if value.isdecimal():
        return int(value)
    return value


def _prepared_yaml_lines(text: str) -> list[tuple[int, str]]:
    lines: list[tuple[int, str]] = []
    for raw_line in text.splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        lines.append((indent, raw_line.strip()))
    return lines


def _parse_yaml_list(
    lines: list[tuple[int, str]],
    index: int,
    indent: int,
) -> tuple[list[Any], int]:
    values: list[Any] = []
    while index < len(lines):
        line_indent, content = lines[index]
        if line_indent < indent:
            break
        if line_indent != indent or not content.startswith("- "):
            raise ValueError("Unsupported YAML list structure")
        values.append(_parse_scalar(content[2:].strip()))
        index += 1
    return values, index


def _parse_yaml_mapping(
    lines: list[tuple[int, str]],
    index: int,
    indent: int,
) -> tuple[dict[str, Any], int]:
    values: dict[str, Any] = {}
    while index < len(lines):
        line_indent, content = lines[index]
        if line_indent < indent:
            break
        if line_indent != indent or content.startswith("- ") or ":" not in content:
            raise ValueError("Unsupported YAML mapping structure")

        key, raw_value = content.split(":", 1)
        key = key.strip()
        value = raw_value.strip()
        index += 1

        if value:
            values[key] = _parse_scalar(value)
            continue

        if index >= len(lines) or lines[index][0] <= line_indent:
            values[key] = None
        elif lines[index][1].startswith("- "):
            values[key], index = _parse_yaml_list(lines, index, line_indent + 2)
        else:
            values[key], index = _parse_yaml_mapping(lines, index, line_indent + 2)

    return values, index


def _load_yaml_mapping(text: str) -> Mapping[str, Any]:
    lines = _prepared_yaml_lines(text)
    payload, index = _parse_yaml_mapping(lines, 0, 0)
    if index != len(lines):
        raise ValueError("Unsupported YAML trailing content")
    return payload


def load_method_definition(
    package: str,
    resource_name: str,
) -> CalculationMethodDefinition:
    """Load a method definition from a package resource YAML file."""
    resource = resources.files(package).joinpath(resource_name)
    loaded = _load_yaml_mapping(resource.read_text(encoding="utf-8"))
    if not loaded:
        raise ValueError(f"Method resource {resource_name} must contain a mapping")
    return parse_method_definition(loaded)


def list_calculation_properties() -> list[CalculationPropertyDefinition]:
    """List supported calculation properties."""
    return sorted(PROPERTY_CATALOG, key=lambda item: item.property_id)


def get_calculation_property(property_id: str) -> CalculationPropertyDefinition:
    """Return a single supported calculation property by ID."""
    for definition in PROPERTY_CATALOG:
        if definition.property_id == property_id:
            return definition
    raise ValueError(f"Unknown calculation property: {property_id}")


def _yaml_resource_names(package: str) -> list[str]:
    package_files = resources.files(package)
    return sorted(
        resource.name
        for resource in package_files.iterdir()
        if resource.name.endswith(".yaml")
    )


def discover_calculation_methods(
    property_packages: Mapping[str, str] | None = None,
) -> list[CalculationMethodDefinition]:
    """Discover and validate YAML method definitions from property packages."""
    packages = property_packages or METHOD_PROPERTY_PACKAGES
    methods: list[CalculationMethodDefinition] = []
    seen: dict[str, str] = {}
    for package_property_id, package in sorted(packages.items()):
        for resource_name in _yaml_resource_names(package):
            try:
                method = load_method_definition(package, resource_name)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid method definition {package}:{resource_name}: {exc}"
                ) from exc
            if method.property_id != package_property_id:
                raise ValueError(
                    f"Method {method.method_id} declares property "
                    f"{method.property_id!r}, expected {package_property_id!r}"
                )
            if method.method_id in seen:
                raise ValueError(
                    f"Duplicate calculation method_id {method.method_id!r} in "
                    f"{package}:{resource_name} and {seen[method.method_id]}"
                )
            seen[method.method_id] = f"{package}:{resource_name}"
            methods.append(method)
    return methods


def list_calculation_methods(
    property_id: str | None = None,
) -> list[CalculationMethodDefinition]:
    """List known calculation methods, optionally filtered by property ID."""
    methods = discover_calculation_methods()
    if property_id is not None:
        methods = [method for method in methods if method.property_id == property_id]
    return methods


def get_calculation_method(
    method_id: str,
    *,
    property_id: str | None = None,
) -> CalculationMethodDefinition:
    """Return a single calculation method by ID."""
    for method in list_calculation_methods(property_id):
        if method.method_id == method_id:
            return method
    raise ValueError(f"Unknown calculation method: {method_id}")
