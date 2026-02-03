"""I/O utilities for loading and validating AnnData objects."""

from .loader import load_h5ad, detect_mappings
from .validator import validate_schema, check_required_fields
from .converter import ensure_spatial_coords, normalize_metadata

__all__ = [
    "load_h5ad",
    "detect_mappings",
    "validate_schema",
    "check_required_fields",
    "ensure_spatial_coords",
    "normalize_metadata",
]
