"""Cluster summary computation."""

import logging
from typing import Optional

import anndata
import pandas as pd

logger = logging.getLogger(__name__)


def compute_cluster_summary(
    adata: anndata.AnnData,
    label_col: str,
    sample_col: Optional[str] = None,
    exclude_na: bool = True,
) -> pd.DataFrame:
    """
    Compute cluster summary statistics.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    label_col : str
        Column in adata.obs containing group labels.
    sample_col : str, optional
        Column in adata.obs containing sample IDs.
        If provided, per-sample counts will be included.
    exclude_na : bool
        If True, exclude cells with NA/missing values in label_col.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with cluster summary statistics.
        Columns: group_id, n_cells, percent_of_total, [sample-specific counts]
    """
    if label_col not in adata.obs.columns:
        raise ValueError(f"Label column '{label_col}' not found in adata.obs")
    
    # Filter out NA values if requested
    if exclude_na:
        mask = ~adata.obs[label_col].isna()
        obs_data = adata.obs[mask]
    else:
        obs_data = adata.obs
    
    total_cells = len(obs_data)
    
    # Compute basic statistics
    group_counts = obs_data[label_col].value_counts()
    
    summary = pd.DataFrame({
        "group_id": group_counts.index,
        "n_cells": group_counts.values,
        "percent_of_total": (group_counts.values / total_cells * 100).round(2),
    })
    
    # Add per-sample counts if sample column is provided
    if sample_col and sample_col in obs_data.columns:
        # Create crosstab
        crosstab = pd.crosstab(obs_data[label_col], obs_data[sample_col])
        
        # Merge with summary
        # Convert crosstab to have group_id as a column
        crosstab_reset = crosstab.reset_index()
        crosstab_reset.rename(columns={label_col: "group_id"}, inplace=True)
        
        # Merge
        summary = summary.merge(crosstab_reset, on="group_id", how="left")
        
        logger.info(f"Added per-sample counts for {len(crosstab.columns)} samples")
    
    # Sort by group_id
    summary = summary.sort_values("group_id").reset_index(drop=True)
    
    logger.info(f"Computed summary for {len(summary)} groups in '{label_col}'")
    
    return summary


def compute_group_composition(
    adata: anndata.AnnData,
    label_col: str,
    group_id: str,
    celltype_col: str,
) -> pd.DataFrame:
    """
    Compute cell type composition for a specific group.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    label_col : str
        Column in adata.obs containing group labels.
    group_id : str
        Specific group ID to analyze.
    celltype_col : str
        Column in adata.obs containing cell type labels.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with cell type composition.
        Columns: cell_type, n_cells, percent
    """
    if label_col not in adata.obs.columns:
        raise ValueError(f"Label column '{label_col}' not found in adata.obs")
    if celltype_col not in adata.obs.columns:
        raise ValueError(f"Cell type column '{celltype_col}' not found in adata.obs")
    
    # Filter to selected group
    group_mask = adata.obs[label_col] == group_id
    group_obs = adata.obs[group_mask]
    
    if len(group_obs) == 0:
        logger.warning(f"Group '{group_id}' has no cells")
        return pd.DataFrame(columns=["cell_type", "n_cells", "percent"])
    
    # Compute composition
    celltype_counts = group_obs[celltype_col].value_counts()
    
    composition = pd.DataFrame({
        "cell_type": celltype_counts.index,
        "n_cells": celltype_counts.values,
        "percent": (celltype_counts.values / len(group_obs) * 100).round(2),
    })
    
    composition = composition.sort_values("n_cells", ascending=False).reset_index(drop=True)
    
    return composition


def compare_qc_metrics(
    adata: anndata.AnnData,
    label_col: str,
    group_id: str,
    qc_columns: list,
) -> pd.DataFrame:
    """
    Compare QC metrics between selected group and all other cells.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    label_col : str
        Column in adata.obs containing group labels.
    group_id : str
        Specific group ID to analyze.
    qc_columns : list of str
        List of QC metric columns to compare.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with QC metric comparisons.
        Columns: metric, mean_in_group, median_in_group, mean_other, median_other
    """
    if label_col not in adata.obs.columns:
        raise ValueError(f"Label column '{label_col}' not found in adata.obs")
    
    # Create masks
    group_mask = adata.obs[label_col] == group_id
    other_mask = (adata.obs[label_col] != group_id) & (~adata.obs[label_col].isna())
    
    group_obs = adata.obs[group_mask]
    other_obs = adata.obs[other_mask]
    
    # Compute metrics
    results = []
    for col in qc_columns:
        if col not in adata.obs.columns:
            logger.warning(f"QC column '{col}' not found, skipping")
            continue
        
        if not pd.api.types.is_numeric_dtype(adata.obs[col]):
            logger.warning(f"QC column '{col}' is not numeric, skipping")
            continue
        
        group_values = group_obs[col].dropna()
        other_values = other_obs[col].dropna()
        
        if len(group_values) == 0 or len(other_values) == 0:
            continue
        
        results.append({
            "metric": col,
            "mean_in_group": group_values.mean(),
            "median_in_group": group_values.median(),
            "mean_other": other_values.mean(),
            "median_other": other_values.median(),
        })
    
    if not results:
        return pd.DataFrame(columns=["metric", "mean_in_group", "median_in_group", "mean_other", "median_other"])
    
    comparison = pd.DataFrame(results)
    
    return comparison
