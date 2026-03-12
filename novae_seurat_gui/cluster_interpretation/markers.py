"""Marker gene computation for cluster interpretation."""

import logging
from typing import Optional

import anndata
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import issparse

logger = logging.getLogger(__name__)

# Constants for numerical stability
EPSILON_FOR_LOG = 1e-9  # Small constant to avoid log(0)


def compute_marker_genes(
    adata: anndata.AnnData,
    label_col: str,
    group_id: str,
    n_genes: int = 25,
    use_layer: Optional[str] = None,
    normalize: bool = True,
) -> pd.DataFrame:
    """
    Compute marker genes for a specific group vs all other cells.
    
    Uses Wilcoxon rank-sum test with Benjamini-Hochberg correction.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    label_col : str
        Column in adata.obs containing group labels.
    group_id : str
        Specific group ID to compute markers for.
    n_genes : int
        Number of top marker genes to return (default: 25).
    use_layer : str, optional
        Layer to use for expression data. If None, uses adata.X.
    normalize : bool
        If True, applies log1p normalization if needed (default: True).
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: gene, logFC, mean_diff, p_value, adj_p_value,
        pct_in_group, pct_out_group. Sorted by adjusted p-value and logFC.
    """
    if label_col not in adata.obs.columns:
        raise ValueError(f"Label column '{label_col}' not found in adata.obs")
    
    # Get expression data
    from .utils import prepare_expression_data
    expr = prepare_expression_data(adata, use_layer=use_layer, normalize=normalize)
    
    # Create group mask (exclude NA/missing values)
    group_mask = (adata.obs[label_col] == group_id).values
    other_mask = (adata.obs[label_col] != group_id) & (~adata.obs[label_col].isna())
    
    n_in = group_mask.sum()
    n_out = other_mask.sum()
    
    if n_in == 0:
        raise ValueError(f"Group '{group_id}' has no cells")
    if n_out == 0:
        raise ValueError("No cells in 'other' group for comparison")
    
    logger.info(f"Computing markers for {group_id}: {n_in} cells vs {n_out} other cells")
    
    # Get expression for each group
    expr_in = expr[group_mask, :]
    expr_out = expr[other_mask, :]
    
    # Compute statistics for each gene
    n_genes_total = expr.shape[1]
    p_values = np.zeros(n_genes_total)
    log_fcs = np.zeros(n_genes_total)
    mean_diffs = np.zeros(n_genes_total)
    pct_in = np.zeros(n_genes_total)
    pct_out = np.zeros(n_genes_total)
    
    for i in range(n_genes_total):
        gene_in = expr_in[:, i]
        gene_out = expr_out[:, i]
        
        # Wilcoxon rank-sum test (Mann-Whitney U)
        try:
            stat, p_val = stats.ranksums(gene_in, gene_out)
            p_values[i] = p_val
        except Exception as e:
            logger.warning(f"Error computing p-value for gene {i}: {e}")
            p_values[i] = 1.0
        
        # Compute mean difference and log fold change
        mean_in = gene_in.mean()
        mean_out = gene_out.mean()
        mean_diffs[i] = mean_in - mean_out
        
        # Log fold change (add small constant to avoid log(0))
        log_fcs[i] = np.log2((mean_in + EPSILON_FOR_LOG) / (mean_out + EPSILON_FOR_LOG))
        
        # Compute percentage of cells expressing (>0)
        pct_in[i] = (gene_in > 0).sum() / n_in * 100
        pct_out[i] = (gene_out > 0).sum() / n_out * 100
    
    # Benjamini-Hochberg correction
    adj_p_values = _benjamini_hochberg_correction(p_values)
    
    # Create results dataframe
    gene_names = adata.var_names.tolist()
    
    results = pd.DataFrame({
        "gene": gene_names,
        "logFC": log_fcs,
        "mean_diff": mean_diffs,
        "p_value": p_values,
        "adj_p_value": adj_p_values,
        "pct_in_group": pct_in,
        "pct_out_group": pct_out,
    })
    
    # Filter for upregulated genes and sort
    results = results[results["logFC"] > 0]  # Only upregulated
    results = results.sort_values(["adj_p_value", "logFC"], ascending=[True, False])
    
    # Return top N
    results = results.head(n_genes).reset_index(drop=True)
    
    logger.info(f"Computed {len(results)} marker genes for group {group_id}")
    
    return results


def _benjamini_hochberg_correction(p_values: np.ndarray) -> np.ndarray:
    """
    Apply Benjamini-Hochberg FDR correction.
    
    Parameters
    ----------
    p_values : np.ndarray
        Array of p-values.
    
    Returns
    -------
    np.ndarray
        Array of adjusted p-values.
    """
    n = len(p_values)
    
    # Sort p-values and keep track of original indices
    sorted_indices = np.argsort(p_values)
    sorted_p_values = p_values[sorted_indices]
    
    # Compute adjusted p-values
    adj_p_values = np.zeros(n)
    
    # BH correction: p_adj[i] = p[i] * n / (i+1)
    for i in range(n):
        adj_p_values[sorted_indices[i]] = min(
            sorted_p_values[i] * n / (i + 1),
            1.0
        )
    
    # Enforce monotonicity (adjusted p-values should be non-decreasing)
    for i in range(n - 2, -1, -1):
        if adj_p_values[sorted_indices[i]] > adj_p_values[sorted_indices[i + 1]]:
            adj_p_values[sorted_indices[i]] = adj_p_values[sorted_indices[i + 1]]
    
    return adj_p_values


def compute_fold_change(
    adata: anndata.AnnData,
    label_col: str,
    group_id: str,
    use_layer: Optional[str] = None,
) -> pd.DataFrame:
    """
    Compute simple fold change for all genes.
    
    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    label_col : str
        Column in adata.obs containing group labels.
    group_id : str
        Specific group ID to compute fold change for.
    use_layer : str, optional
        Layer to use for expression data.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: gene, logFC, mean_in, mean_out.
    """
    from .utils import prepare_expression_data
    expr = prepare_expression_data(adata, use_layer=use_layer, normalize=True)
    
    group_mask = (adata.obs[label_col] == group_id).values
    other_mask = (adata.obs[label_col] != group_id) & (~adata.obs[label_col].isna())
    
    expr_in = expr[group_mask, :]
    expr_out = expr[other_mask, :]
    
    mean_in = expr_in.mean(axis=0)
    mean_out = expr_out.mean(axis=0)
    
    # Convert to 1D arrays if needed
    if hasattr(mean_in, "A1"):
        mean_in = mean_in.A1
    if hasattr(mean_out, "A1"):
        mean_out = mean_out.A1
    
    log_fcs = np.log2((mean_in + EPSILON_FOR_LOG) / (mean_out + EPSILON_FOR_LOG))
    
    results = pd.DataFrame({
        "gene": adata.var_names.tolist(),
        "logFC": log_fcs,
        "mean_in": mean_in,
        "mean_out": mean_out,
    })
    
    return results.sort_values("logFC", ascending=False)


def compute_domain_weights(
    adata: anndata.AnnData,
    label_col: str,
    group_id: str,
    n_top: int = 20,
    use_layer: Optional[str] = None,
    store_in_uns: bool = True,
) -> pd.DataFrame:
    """
    Compute proxy gene weights that drive a domain/cluster assignment.

    When the real Novae package is not available this uses log2 fold change
    (group vs. rest) as a proxy weight.  Both positive and negative
    contributors are returned so the caller can display top positive, top
    negative and top absolute drivers.

    Results are cached in ``adata.uns['xspatialnovae_domain_weights']`` for
    reproducibility and fast repeated access.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    label_col : str
        Column in ``adata.obs`` containing domain/cluster labels.
    group_id : str
        Specific domain/cluster to explain.
    n_top : int
        Number of top genes to return per direction.  The returned table
        contains at most ``2 * n_top`` rows (top positive + top negative by
        ``logFC``).
    use_layer : str, optional
        Layer to use for expression data.  Defaults to ``adata.X``.
    store_in_uns : bool
        If True, save the full ranked list in
        ``adata.uns['xspatialnovae_domain_weights'][label_col][str(group_id)]``.

    Returns
    -------
    pd.DataFrame
        DataFrame sorted by absolute ``logFC`` descending, with columns:
        ``gene``, ``logFC``, ``mean_in``, ``mean_out``, ``direction``
        (``"positive"`` | ``"negative"``).
    """
    if label_col not in adata.obs.columns:
        raise ValueError(f"Label column '{label_col}' not found in adata.obs")

    from .utils import prepare_expression_data

    expr = prepare_expression_data(adata, use_layer=use_layer, normalize=True)

    group_mask = (adata.obs[label_col] == group_id).values
    other_mask = (adata.obs[label_col] != group_id) & (~adata.obs[label_col].isna())

    n_in = int(group_mask.sum())
    n_out = int(other_mask.sum())

    if n_in == 0:
        raise ValueError(f"Group '{group_id}' has no cells")
    if n_out == 0:
        raise ValueError("No cells in 'other' group for comparison")

    logger.info(
        "Computing domain weights for %s: %d cells vs %d other cells",
        group_id,
        n_in,
        n_out,
    )

    expr_in = expr[group_mask, :]
    expr_out = expr[other_mask, :]

    mean_in = np.asarray(expr_in.mean(axis=0)).ravel()
    mean_out = np.asarray(expr_out.mean(axis=0)).ravel()

    log_fcs = np.log2((mean_in + EPSILON_FOR_LOG) / (mean_out + EPSILON_FOR_LOG))

    all_results = pd.DataFrame(
        {
            "gene": adata.var_names.tolist(),
            "logFC": log_fcs,
            "mean_in": mean_in,
            "mean_out": mean_out,
        }
    )
    all_results["direction"] = np.where(all_results["logFC"] >= 0, "positive", "negative")

    # Persist full list before subsetting
    if store_in_uns:
        if "xspatialnovae_domain_weights" not in adata.uns:
            adata.uns["xspatialnovae_domain_weights"] = {}
        if label_col not in adata.uns["xspatialnovae_domain_weights"]:
            adata.uns["xspatialnovae_domain_weights"][label_col] = {}
        adata.uns["xspatialnovae_domain_weights"][label_col][str(group_id)] = (
            all_results.to_dict(orient="records")
        )

    # Return top positive + top negative sorted by absolute logFC
    top_pos = (
        all_results[all_results["logFC"] >= 0]
        .sort_values("logFC", ascending=False)
        .head(n_top)
    )
    top_neg = (
        all_results[all_results["logFC"] < 0]
        .sort_values("logFC", ascending=True)
        .head(n_top)
    )

    combined = pd.concat([top_pos, top_neg], ignore_index=True)
    combined = combined.sort_values("logFC", key=np.abs, ascending=False).reset_index(drop=True)

    logger.info(
        "Domain weights for %s: %d positive, %d negative drivers returned",
        group_id,
        len(top_pos),
        len(top_neg),
    )

    return combined

