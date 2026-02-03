"""Validator for H5AD schema and required fields."""

import logging
from typing import Dict, List, Optional, Tuple

import anndata
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Custom exception for validation errors."""

    pass


def validate_schema(adata: anndata.AnnData, strict: bool = False) -> Tuple[bool, List[str]]:
    """
    Validate that the AnnData object has the required schema.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    strict : bool
        If True, enforce stricter requirements.

    Returns
    -------
    tuple of (bool, list)
        (is_valid, list of warning/error messages)
    """
    messages = []
    is_valid = True

    # Check basic structure
    if adata.n_obs == 0:
        messages.append("ERROR: No cells (observations) in the dataset.")
        is_valid = False

    if adata.n_vars == 0:
        messages.append("ERROR: No features (variables) in the dataset.")
        is_valid = False

    # Check for data matrix
    if adata.X is None and not adata.layers:
        messages.append("ERROR: No data matrix found (neither adata.X nor adata.layers).")
        is_valid = False

    # Check obs
    if adata.obs.empty:
        messages.append("WARNING: adata.obs is empty (no metadata).")

    # Check for cell_id or index
    if not adata.obs.index.is_unique:
        messages.append("ERROR: Cell IDs (obs.index) are not unique.")
        is_valid = False

    # Check for spatial coordinates
    has_spatial = False
    if "spatial" in adata.obsm or "X_spatial" in adata.obsm:
        has_spatial = True
        spatial_key = "spatial" if "spatial" in adata.obsm else "X_spatial"
        spatial_coords = adata.obsm[spatial_key]

        if spatial_coords.shape[1] != 2:
            messages.append(
                f"ERROR: Spatial coordinates in obsm['{spatial_key}'] should have 2 columns (x, y), "
                f"found {spatial_coords.shape[1]}."
            )
            is_valid = False

        if np.any(np.isnan(spatial_coords)):
            messages.append(
                f"WARNING: Spatial coordinates contain NaN values in obsm['{spatial_key}']."
            )

    if not has_spatial and strict:
        messages.append("ERROR: No spatial coordinates found in adata.obsm.")
        is_valid = False
    elif not has_spatial:
        messages.append(
            "WARNING: No spatial coordinates found. Expected 'spatial' or 'X_spatial' in adata.obsm."
        )

    logger.info(f"Validation completed: {'PASSED' if is_valid else 'FAILED'}")
    for msg in messages:
        if msg.startswith("ERROR"):
            logger.error(msg)
        else:
            logger.warning(msg)

    return is_valid, messages


def check_required_fields(
    adata: anndata.AnnData,
    required_obs_cols: Optional[List[str]] = None,
    required_obsm_keys: Optional[List[str]] = None,
) -> Tuple[bool, Dict[str, List[str]]]:
    """
    Check for required fields in adata.obs and adata.obsm.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    required_obs_cols : list of str, optional
        List of required column names in adata.obs.
    required_obsm_keys : list of str, optional
        List of required keys in adata.obsm.

    Returns
    -------
    tuple of (bool, dict)
        (all_present, dict with 'missing_obs' and 'missing_obsm' lists)
    """
    missing = {"missing_obs": [], "missing_obsm": []}

    if required_obs_cols:
        for col in required_obs_cols:
            if col not in adata.obs.columns:
                missing["missing_obs"].append(col)

    if required_obsm_keys:
        for key in required_obsm_keys:
            if key not in adata.obsm:
                missing["missing_obsm"].append(key)

    all_present = not missing["missing_obs"] and not missing["missing_obsm"]

    if not all_present:
        logger.warning(f"Missing required fields: {missing}")

    return all_present, missing


def check_data_types(adata: anndata.AnnData) -> Dict[str, str]:
    """
    Check and report data types for key fields.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.

    Returns
    -------
    dict
        Dictionary mapping field names to data types.
    """
    dtypes = {}

    if adata.X is not None:
        dtypes["X"] = str(adata.X.dtype)

    for layer_name, layer_data in adata.layers.items():
        dtypes[f"layers/{layer_name}"] = str(layer_data.dtype)

    for obsm_key, obsm_data in adata.obsm.items():
        dtypes[f"obsm/{obsm_key}"] = str(obsm_data.dtype)

    return dtypes


def check_counts_data(adata: anndata.AnnData, layer: Optional[str] = None) -> Dict:
    """
    Check if the data looks like raw counts.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    layer : str, optional
        Layer to check. If None, checks adata.X.

    Returns
    -------
    dict
        Dictionary with statistics about the data:
        - 'is_integer': whether all values are integers
        - 'has_negative': whether there are negative values
        - 'max_value': maximum value
        - 'mean_value': mean value
        - 'sparsity': fraction of zero values
    """
    if layer is not None:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' not found in adata.layers")
        data = adata.layers[layer]
    else:
        data = adata.X

    if data is None:
        raise ValueError("No data matrix found")

    # Convert to dense if sparse
    import scipy.sparse as sp

    if sp.issparse(data):
        # Sample for efficiency
        sample_size = min(10000, data.shape[0] * data.shape[1])
        if data.shape[0] * data.shape[1] > sample_size:
            row_idx = np.random.choice(data.shape[0], size=100, replace=False)
            sample_data = data[row_idx, :].toarray().flatten()
        else:
            sample_data = data.toarray().flatten()
    else:
        sample_data = data.flatten()

    stats = {
        "is_integer": np.allclose(sample_data, np.round(sample_data)),
        "has_negative": np.any(sample_data < 0),
        "max_value": float(np.max(sample_data)),
        "mean_value": float(np.mean(sample_data)),
        "sparsity": float(np.mean(sample_data == 0)),
    }

    return stats
