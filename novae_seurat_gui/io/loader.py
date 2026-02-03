"""Loader for H5AD files with automatic mapping detection."""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata
import pandas as pd

logger = logging.getLogger(__name__)


def load_h5ad(file_path: str) -> anndata.AnnData:
    """
    Load an H5AD file.

    Parameters
    ----------
    file_path : str
        Path to H5AD file.

    Returns
    -------
    anndata.AnnData
        Loaded AnnData object.
    """
    logger.info(f"Loading H5AD file: {file_path}")
    adata = anndata.read_h5ad(file_path)
    logger.info(
        f"Loaded {adata.n_obs} cells × {adata.n_vars} features from {file_path}"
    )
    return adata


def detect_mappings(adata: anndata.AnnData) -> Dict[str, Optional[str]]:
    """
    Auto-detect column mappings for coordinates, sample_id, and cell_type.

    Looks for common column names in adata.obs and spatial coordinates in adata.obsm.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.

    Returns
    -------
    dict
        Dictionary with detected mappings:
        - 'x_col': x coordinate column
        - 'y_col': y coordinate column
        - 'sample_id_col': sample identifier column
        - 'cell_type_col': cell type column (optional)
        - 'confidence_col': confidence/probability column (optional)
        - 'spatial_key': key in obsm for spatial coords
        - 'units': coordinate units ('mm' or 'px')
    """
    mappings = {
        "x_col": None,
        "y_col": None,
        "sample_id_col": None,
        "cell_type_col": None,
        "confidence_col": None,
        "spatial_key": None,
        "units": None,
    }

    obs_cols = adata.obs.columns.tolist()

    # Detect spatial coordinates in obsm
    if "spatial" in adata.obsm:
        mappings["spatial_key"] = "spatial"
        logger.info("Detected spatial coordinates in adata.obsm['spatial']")
    elif "X_spatial" in adata.obsm:
        mappings["spatial_key"] = "X_spatial"
        logger.info("Detected spatial coordinates in adata.obsm['X_spatial']")

    # Detect x,y coordinate columns in obs
    x_candidates = [
        "x_slide_mm",
        "x_location",
        "x_centroid",
        "x_FOV_px",
        "x",
        "X",
    ]
    y_candidates = [
        "y_slide_mm",
        "y_location",
        "y_centroid",
        "y_FOV_px",
        "y",
        "Y",
    ]

    for x_col in x_candidates:
        if x_col in obs_cols:
            mappings["x_col"] = x_col
            # Infer units from column name
            if "_mm" in x_col:
                mappings["units"] = "mm"
            elif "_px" in x_col:
                mappings["units"] = "px"
            else:
                mappings["units"] = "unknown"
            logger.info(f"Detected x coordinate column: {x_col}")
            break

    for y_col in y_candidates:
        if y_col in obs_cols:
            mappings["y_col"] = y_col
            logger.info(f"Detected y coordinate column: {y_col}")
            break

    # Detect sample_id column
    sample_id_candidates = ["SampleID", "sample_id", "Tissue", "tissue", "sample"]
    for col in sample_id_candidates:
        if col in obs_cols:
            mappings["sample_id_col"] = col
            logger.info(f"Detected sample_id column: {col}")
            break

    # Detect cell_type column
    cell_type_candidates = [
        "cell_type_1",
        "cell_type",
        "celltype",
        "CellType",
        "cell_label",
    ]
    for col in cell_type_candidates:
        if col in obs_cols:
            mappings["cell_type_col"] = col
            logger.info(f"Detected cell_type column: {col}")
            break

    # Detect confidence/probability column
    confidence_candidates = [
        "posterior_probability",
        "confidence",
        "probability",
        "cell_type_confidence",
    ]
    for col in confidence_candidates:
        if col in obs_cols:
            mappings["confidence_col"] = col
            logger.info(f"Detected confidence column: {col}")
            break

    return mappings


def get_available_columns(adata: anndata.AnnData) -> Dict[str, List[str]]:
    """
    Get lists of available columns for user selection.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.

    Returns
    -------
    dict
        Dictionary with lists of available columns:
        - 'obs_columns': all observation columns
        - 'numeric_columns': numeric observation columns
        - 'categorical_columns': categorical observation columns
        - 'obsm_keys': keys in obsm
        - 'layers': available layers
    """
    obs_cols = adata.obs.columns.tolist()
    numeric_cols = adata.obs.select_dtypes(include=["number"]).columns.tolist()
    categorical_cols = adata.obs.select_dtypes(
        include=["object", "category"]
    ).columns.tolist()

    return {
        "obs_columns": obs_cols,
        "numeric_columns": numeric_cols,
        "categorical_columns": categorical_cols,
        "obsm_keys": list(adata.obsm.keys()),
        "layers": list(adata.layers.keys()) if adata.layers else [],
    }


def summarize_adata(adata: anndata.AnnData) -> Dict:
    """
    Generate a summary of the AnnData object.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.

    Returns
    -------
    dict
        Summary statistics and metadata.
    """
    summary = {
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "obs_columns": adata.obs.columns.tolist(),
        "obsm_keys": list(adata.obsm.keys()),
        "layers": list(adata.layers.keys()) if adata.layers else [],
        "uns_keys": list(adata.uns.keys()) if adata.uns else [],
    }

    # Sample statistics
    if "sample_id" in adata.obs.columns:
        summary["n_samples"] = adata.obs["sample_id"].nunique()
        summary["cells_per_sample"] = adata.obs["sample_id"].value_counts().to_dict()

    return summary
