"""Compute neighborhood composition features for niche analysis."""

import logging
from typing import Optional

import anndata
import numpy as np
import pandas as pd
from scipy import sparse

logger = logging.getLogger(__name__)


def compute_neighborhood_composition(
    adata: anndata.AnnData,
    cell_type_col: str = "cell_type",
    confidence_col: Optional[str] = None,
    connectivities_key: str = "spatial_connectivities",
    normalize: bool = True,
) -> pd.DataFrame:
    """
    Compute neighborhood composition features for each cell.

    For each cell, computes the proportion of each cell type among its spatial neighbors.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with spatial neighbor graph.
    cell_type_col : str
        Column in adata.obs containing cell type labels.
    confidence_col : str, optional
        Column in adata.obs containing confidence weights (e.g., posterior_probability).
        If provided, uses weighted composition.
    connectivities_key : str
        Key in adata.obsp containing connectivity matrix.
    normalize : bool
        If True, normalize composition to sum to 1 for each cell.

    Returns
    -------
    pd.DataFrame
        DataFrame with rows=cells and columns=cell types, containing composition fractions.
    """
    if cell_type_col not in adata.obs.columns:
        raise ValueError(f"Cell type column '{cell_type_col}' not found in adata.obs")

    if connectivities_key not in adata.obsp:
        raise ValueError(f"Connectivity key '{connectivities_key}' not found in adata.obsp")

    logger.info(f"Computing neighborhood composition from '{cell_type_col}'")

    connectivities = adata.obsp[connectivities_key]
    cell_types = adata.obs[cell_type_col].values
    unique_types = np.unique(cell_types)

    # Create cell type indicator matrix
    n_cells = len(cell_types)
    n_types = len(unique_types)

    type_to_idx = {ct: i for i, ct in enumerate(unique_types)}
    cell_type_matrix = np.zeros((n_cells, n_types))

    for i, ct in enumerate(cell_types):
        cell_type_matrix[i, type_to_idx[ct]] = 1.0

    # Apply confidence weights if provided
    if confidence_col and confidence_col in adata.obs.columns:
        logger.info(f"Using confidence weights from '{confidence_col}'")
        confidence = adata.obs[confidence_col].values
        cell_type_matrix *= confidence[:, np.newaxis]

    # Compute neighborhood composition
    # For each cell, sum the cell types of its neighbors
    if sparse.issparse(connectivities):
        neighbor_composition = connectivities.dot(cell_type_matrix)
    else:
        neighbor_composition = connectivities @ cell_type_matrix

    # Normalize
    if normalize:
        row_sums = neighbor_composition.sum(axis=1)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        neighbor_composition = neighbor_composition / row_sums[:, np.newaxis]

    # Convert to DataFrame
    composition_df = pd.DataFrame(
        neighbor_composition,
        index=adata.obs.index,
        columns=[f"niche_comp_{ct}" for ct in unique_types],
    )

    logger.info(f"Computed neighborhood composition for {n_types} cell types")

    return composition_df


def add_composition_to_adata(
    adata: anndata.AnnData,
    composition_df: pd.DataFrame,
    obsm_key: str = "neighborhood_composition",
) -> anndata.AnnData:
    """
    Add neighborhood composition to adata.obsm.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    composition_df : pd.DataFrame
        Composition dataframe from compute_neighborhood_composition.
    obsm_key : str
        Key to store composition in adata.obsm.

    Returns
    -------
    anndata.AnnData
        Modified AnnData object.
    """
    adata.obsm[obsm_key] = composition_df.values
    adata.uns[f"{obsm_key}_columns"] = composition_df.columns.tolist()

    logger.info(f"Added neighborhood composition to adata.obsm['{obsm_key}']")

    return adata


def compute_diversity_metrics(composition_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute diversity metrics from neighborhood composition.

    Parameters
    ----------
    composition_df : pd.DataFrame
        Composition dataframe.

    Returns
    -------
    pd.DataFrame
        DataFrame with diversity metrics per cell.
    """
    # Shannon entropy
    def shannon_entropy(row):
        p = row[row > 0]
        return -np.sum(p * np.log(p))

    # Simpson diversity
    def simpson_diversity(row):
        return 1 - np.sum(row ** 2)

    # Richness (number of cell types present)
    def richness(row):
        return np.sum(row > 0)

    diversity_df = pd.DataFrame(index=composition_df.index)
    diversity_df["shannon_entropy"] = composition_df.apply(shannon_entropy, axis=1)
    diversity_df["simpson_diversity"] = composition_df.apply(simpson_diversity, axis=1)
    diversity_df["richness"] = composition_df.apply(richness, axis=1)

    return diversity_df
