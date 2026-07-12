"""Calculation method definitions and registries."""

from grimperium.calculation.methods.feature_schema import (
    FEATURE_SCHEMA_HASH,
    FEATURE_SCHEMA_ID,
    FeatureSchemaDefinition,
    calculate_feature_schema_hash,
    get_feature_schema,
    validate_feature_schema,
)
from grimperium.calculation.methods.registry import (
    CalculationMethodDefinition,
    CalculationPropertyDefinition,
    CompatibilityDefinition,
    ConformerSelectionDefinition,
    ModelRequirementDefinition,
    XtbDefinition,
    discover_calculation_methods,
    get_calculation_method,
    get_calculation_property,
    list_calculation_methods,
    list_calculation_properties,
    load_method_definition,
    parse_method_definition,
)

__all__ = [
    "FEATURE_SCHEMA_HASH",
    "FEATURE_SCHEMA_ID",
    "CalculationMethodDefinition",
    "CalculationPropertyDefinition",
    "CompatibilityDefinition",
    "ConformerSelectionDefinition",
    "FeatureSchemaDefinition",
    "ModelRequirementDefinition",
    "XtbDefinition",
    "calculate_feature_schema_hash",
    "discover_calculation_methods",
    "get_calculation_method",
    "get_calculation_property",
    "get_feature_schema",
    "list_calculation_properties",
    "list_calculation_methods",
    "load_method_definition",
    "parse_method_definition",
    "validate_feature_schema",
]
