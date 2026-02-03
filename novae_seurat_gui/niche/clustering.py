"""Clustering neighborhoods into niche states."""

import logging
from typing import Literal, Optional

import anndata
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


def cluster_niches(
    composition_df: pd.DataFrame,
    n_niches: int = 10,
    method: Literal["kmeans", "hierarchical"] = "kmeans",
    standardize: bool = True,
    random_state: int = 42,
) -> np.ndarray:
    """
    Cluster neighborhood composition vectors into niche states.

    Parameters
    ----------
    composition_df : pd.DataFrame
        Neighborhood composition dataframe.
    n_niches : int
        Number of niche clusters.
    method : {'kmeans', 'hierarchical'}
        Clustering method.
    standardize : bool
        If True, standardize features before clustering.
    random_state : int
        Random seed.

    Returns
    -------
    np.ndarray
        Array of niche labels.
    """
    logger.info(f"Clustering neighborhoods into {n_niches} niches using {method}")

    X = composition_df.values

    # Standardize
    if standardize:
        scaler = StandardScaler()
        X = scaler.fit_transform(X)

    # Cluster
    if method == "kmeans":
        clusterer = KMeans(n_clusters=n_niches, random_state=random_state, n_init=10)
        labels = clusterer.fit_predict(X)
    elif method == "hierarchical":
        clusterer = AgglomerativeClustering(n_clusters=n_niches)
        labels = clusterer.fit_predict(X)
    else:
        raise ValueError(f"Unknown clustering method: {method}")

    logger.info(f"Identified {len(np.unique(labels))} niche clusters")

    return labels


def assign_niche_labels(
    adata: anndata.AnnData,
    composition_df: pd.DataFrame,
    n_niches: int = 10,
    method: Literal["kmeans", "hierarchical"] = "kmeans",
    niche_col_name: str = "niche",
    random_state: int = 42,
) -> anndata.AnnData:
    """
    Cluster niches and assign labels to adata.obs.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    composition_df : pd.DataFrame
        Neighborhood composition dataframe.
    n_niches : int
        Number of niche clusters.
    method : str
        Clustering method.
    niche_col_name : str
        Column name for niche labels in adata.obs.
    random_state : int
        Random seed.

    Returns
    -------
    anndata.AnnData
        Modified AnnData object with niche labels.
    """
    labels = cluster_niches(
        composition_df, n_niches=n_niches, method=method, random_state=random_state
    )

    adata.obs[niche_col_name] = labels.astype(str)
    adata.obs[niche_col_name] = adata.obs[niche_col_name].astype("category")

    logger.info(f"Assigned niche labels to adata.obs['{niche_col_name}']")

    return adata


def summarize_niche_composition(
    adata: anndata.AnnData,
    composition_df: pd.DataFrame,
    niche_col: str = "niche",
    cell_type_col: str = "cell_type",
) -> pd.DataFrame:
    """
    Summarize cell type composition for each niche.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with niche labels.
    composition_df : pd.DataFrame
        Neighborhood composition dataframe.
    niche_col : str
        Column in adata.obs containing niche labels.
    cell_type_col : str
        Column in adata.obs containing cell type labels.

    Returns
    -------
    pd.DataFrame
        DataFrame with niche × cell type composition.
    """
    if niche_col not in adata.obs.columns:
        raise ValueError(f"Niche column '{niche_col}' not found in adata.obs")

    logger.info(f"Summarizing composition for niches in '{niche_col}'")

    # Merge niche labels with composition
    niche_labels = adata.obs[niche_col].values
    composition_with_niche = composition_df.copy()
    composition_with_niche["niche"] = niche_labels

    # Average composition per niche
    niche_summary = composition_with_niche.groupby("niche").mean()

    # Sort columns for readability
    niche_summary = niche_summary.sort_index(axis=1)

    return niche_summary


def compute_niche_enrichment(
    adata: anndata.AnnData,
    niche_col: str = "niche",
    cell_type_col: str = "cell_type",
) -> pd.DataFrame:
    """
    Compute enrichment of cell types within each niche.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    niche_col : str
        Column containing niche labels.
    cell_type_col : str
        Column containing cell type labels.

    Returns
    -------
    pd.DataFrame
        DataFrame with enrichment scores (log odds ratio).
    """
    if niche_col not in adata.obs.columns:
        raise ValueError(f"Niche column '{niche_col}' not found in adata.obs")
    if cell_type_col not in adata.obs.columns:
        raise ValueError(f"Cell type column '{cell_type_col}' not found in adata.obs")

    logger.info("Computing niche enrichment for cell types")

    # Crosstab
    contingency = pd.crosstab(adata.obs[niche_col], adata.obs[cell_type_col])

    # Compute proportions
    niche_totals = contingency.sum(axis=1)
    type_totals = contingency.sum(axis=0)
    total = contingency.sum().sum()

    # Enrichment (log odds ratio)
    enrichment = pd.DataFrame(index=contingency.index, columns=contingency.columns)

    for niche in contingency.index:
        for cell_type in contingency.columns:
            observed = contingency.loc[niche, cell_type]
            expected = (niche_totals[niche] * type_totals[cell_type]) / total

            if expected > 0:
                enrichment.loc[niche, cell_type] = np.log2((observed + 1) / (expected + 1))
            else:
                enrichment.loc[niche, cell_type] = 0

    return enrichment.astype(float)


def export_niche_results(
    adata: anndata.AnnData,
    output_file: str,
    niche_col: str = "niche",
    cell_id_col: str = "cell_id",
) -> None:
    """
    Export niche assignments to CSV.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with niche labels.
    output_file : str
        Output CSV file path.
    niche_col : str
        Column containing niche labels.
    cell_id_col : str
        Column containing cell IDs.
    """
    logger.info(f"Exporting niche assignments to {output_file}")

    if niche_col not in adata.obs.columns:
        raise ValueError(f"Niche column '{niche_col}' not found in adata.obs")

    # Prepare export
    export_cols = [cell_id_col, niche_col]
    if "sample_id" in adata.obs.columns:
        export_cols.insert(1, "sample_id")

    export_df = adata.obs[export_cols].copy()
    export_df.to_csv(output_file, index=False)

    logger.info(f"Exported {len(export_df)} cells with niche labels")
