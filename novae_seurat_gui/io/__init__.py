"""I/O utilities for loading and validating AnnData objects."""

from .loader import load_h5ad, detect_mappings
from .validator import validate_schema, check_required_fields
from .converter import ensure_spatial_coords, normalize_metadata
from .convert import (
    check_r_available,
    check_r_packages,
    convert_rds_to_h5ad,
    ensure_anndata_schema,
    convert_rds_to_h5ad_with_validation,
)

__all__ = [
    "load_h5ad",
    "detect_mappings",
    "validate_schema",
    "check_required_fields",
    "ensure_spatial_coords",
    "normalize_metadata",
    "check_r_available",
    "check_r_packages",
    "convert_rds_to_h5ad",
    "ensure_anndata_schema",
    "convert_rds_to_h5ad_with_validation",
]
