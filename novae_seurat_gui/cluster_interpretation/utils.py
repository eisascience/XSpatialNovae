"""Utility functions for cluster interpretation."""

import logging
from typing import List, Optional, Tuple

import anndata
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def get_candidate_label_columns(adata: anndata.AnnData) -> List[str]:
    """
    Get candidate label/grouping columns from adata.obs.
    
    Prioritizes columns containing keywords: domain, cluster, leiden, louvain, cell_type.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    
    Returns
    -------
    list of str
        List of candidate column names, with prioritized ones first.
    """
    keywords = ["domain", "cluster", "leiden", "louvain", "cell_type"]
    
    # Find columns matching keywords
    priority_cols = []
    for col in adata.obs.columns:
        col_lower = col.lower()
        if any(kw in col_lower for kw in keywords):
            priority_cols.append(col)
    
    # Remove duplicates while preserving order
    priority_cols = list(dict.fromkeys(priority_cols))
    
    # Get all other columns (excluding priority ones)
    all_cols = [col for col in adata.obs.columns if col not in priority_cols]
    
    # Return priority columns first, then all others
    return priority_cols + all_cols


def prepare_expression_data(
    adata: anndata.AnnData,
    use_layer: Optional[str] = None,
    normalize: bool = True,
) -> np.ndarray:
    """
    Prepare expression data for marker gene computation.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    use_layer : str, optional
        Layer to use for expression data. If None, uses adata.X.
    normalize : bool
        If True, applies log1p normalization (if data is not already normalized).
    
    Returns
    -------
    np.ndarray
        Expression matrix (cells × genes).
    """
    # Get expression data
    if use_layer is not None and use_layer in adata.layers:
        expr = adata.layers[use_layer]
        logger.info(f"Using layer '{use_layer}' for expression data")
    else:
        expr = adata.X
        logger.info("Using adata.X for expression data")
    
    # Convert to dense if sparse
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    
    # Normalize if requested and data looks like raw counts
    if normalize:
        # Simple heuristic: if max value is large, assume raw counts
        if expr.max() > 100:
            logger.info("Applying log1p normalization")
            # Normalize to counts per 10k, then log1p
            expr = expr.copy()
            # Normalize each cell
            cell_sums = expr.sum(axis=1, keepdims=True)
            cell_sums[cell_sums == 0] = 1  # Avoid division by zero
            expr = expr / cell_sums * 1e4
            expr = np.log1p(expr)
        else:
            logger.info("Data appears normalized, skipping normalization")
    
    return expr


def check_spatial_coords(adata: anndata.AnnData) -> Tuple[bool, Optional[str]]:
    """
    Check if spatial coordinates are available.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    
    Returns
    -------
    has_spatial : bool
        True if spatial coordinates are available.
    message : str or None
        Message describing the status or issue.
    """
    if "spatial" in adata.obsm:
        return True, "Spatial coordinates found in adata.obsm['spatial']"
    
    # Check for x/y columns in obs
    x_candidates = [col for col in adata.obs.columns if col.lower() in ['x', 'x_coord', 'x_coordinate']]
    y_candidates = [col for col in adata.obs.columns if col.lower() in ['y', 'y_coord', 'y_coordinate']]
    
    if x_candidates and y_candidates:
        return True, f"Spatial coordinates found in obs: {x_candidates[0]}, {y_candidates[0]}"
    
    return False, "No spatial coordinates found. Please ensure adata.obsm['spatial'] exists or x/y columns are in adata.obs."


def get_sample_column(adata: anndata.AnnData) -> Optional[str]:
    """
    Detect sample ID column in adata.obs.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    
    Returns
    -------
    str or None
        Name of sample column, or None if not found.
    """
    sample_candidates = ["sample_id", "SampleID", "sample", "Sample", "fov", "FOV"]
    
    for candidate in sample_candidates:
        if candidate in adata.obs.columns:
            return candidate
    
    return None


def get_celltype_column(adata: anndata.AnnData) -> Optional[str]:
    """
    Detect cell type column in adata.obs.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    
    Returns
    -------
    str or None
        Name of cell type column, or None if not found.
    """
    for col in adata.obs.columns:
        if "cell_type" in col.lower():
            return col
    
    return None


def get_qc_columns(adata: anndata.AnnData) -> List[str]:
    """
    Get QC metric columns from adata.obs.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    
    Returns
    -------
    list of str
        List of QC column names.
    """
    qc_keywords = [
        "ncount", "nfeature", "n_genes", "n_counts",
        "area", "area.um2", "posterior_probability",
        "pct_", "percent_", "mito", "ribo"
    ]
    
    qc_cols = []
    for col in adata.obs.columns:
        col_lower = col.lower()
        if any(kw in col_lower for kw in qc_keywords):
            # Only include numeric columns
            if pd.api.types.is_numeric_dtype(adata.obs[col]):
                qc_cols.append(col)
    
    return qc_cols
