"""QC summary and statistics functions."""

import logging
from typing import Dict, List, Optional

import anndata
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def compute_cell_statistics(adata: anndata.AnnData, layer: Optional[str] = None) -> pd.DataFrame:
    """
    Compute per-cell statistics from expression data.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    layer : str, optional
        Layer to use. If None, uses adata.X.

    Returns
    -------
    pd.DataFrame
        DataFrame with per-cell statistics.
    """
    if layer is not None:
        data = adata.layers[layer]
    else:
        data = adata.X

    import scipy.sparse as sp

    if sp.issparse(data):
        n_counts = np.array(data.sum(axis=1)).flatten()
        n_features = np.array((data > 0).sum(axis=1)).flatten()
    else:
        n_counts = data.sum(axis=1)
        n_features = (data > 0).sum(axis=1)

    stats_df = pd.DataFrame(
        {
            "n_counts": n_counts,
            "n_features": n_features,
        },
        index=adata.obs.index,
    )

    return stats_df


def compute_qc_summary(
    adata: anndata.AnnData,
    qc_columns: Optional[List[str]] = None,
    group_by: Optional[str] = None,
) -> Dict:
    """
    Compute summary statistics for QC metrics.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    qc_columns : list of str, optional
        List of QC metric columns to summarize. If None, auto-detects numeric columns.
    group_by : str, optional
        Column name to group by (e.g., 'sample_id').

    Returns
    -------
    dict
        Dictionary containing summary statistics.
    """
    summary = {}

    # Auto-detect QC columns if not specified
    if qc_columns is None:
        numeric_cols = adata.obs.select_dtypes(include=[np.number]).columns.tolist()
        qc_columns = [
            col
            for col in numeric_cols
            if any(
                keyword in col.lower()
                for keyword in ["count", "feature", "area", "qc", "probability"]
            )
        ]

    summary["qc_columns"] = qc_columns
    summary["n_cells"] = adata.n_obs

    # Overall statistics
    summary["overall"] = {}
    for col in qc_columns:
        if col in adata.obs.columns:
            col_data = adata.obs[col]
            summary["overall"][col] = {
                "mean": float(col_data.mean()),
                "median": float(col_data.median()),
                "std": float(col_data.std()),
                "min": float(col_data.min()),
                "max": float(col_data.max()),
                "q25": float(col_data.quantile(0.25)),
                "q75": float(col_data.quantile(0.75)),
            }

    # Group-wise statistics
    if group_by and group_by in adata.obs.columns:
        summary["by_group"] = {}
        for group_name, group_data in adata.obs.groupby(group_by):
            summary["by_group"][str(group_name)] = {}
            summary["by_group"][str(group_name)]["n_cells"] = len(group_data)

            for col in qc_columns:
                if col in group_data.columns:
                    col_data = group_data[col]
                    summary["by_group"][str(group_name)][col] = {
                        "mean": float(col_data.mean()),
                        "median": float(col_data.median()),
                    }

    return summary


def compute_filter_stats(
    adata: anndata.AnnData, mask: np.ndarray, group_by: Optional[str] = None
) -> Dict:
    """
    Compute statistics about filtering results.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    mask : np.ndarray
        Boolean mask indicating cells that passed filters.
    group_by : str, optional
        Column name to group by for per-group statistics.

    Returns
    -------
    dict
        Dictionary with filtering statistics.
    """
    n_kept = np.sum(mask)
    n_filtered = np.sum(~mask)
    total = len(mask)

    stats = {
        "n_total": total,
        "n_kept": int(n_kept),
        "n_filtered": int(n_filtered),
        "percent_kept": float(100 * n_kept / total),
        "percent_filtered": float(100 * n_filtered / total),
    }

    if group_by and group_by in adata.obs.columns:
        stats["by_group"] = {}
        for group_name in adata.obs[group_by].unique():
            group_mask = adata.obs[group_by] == group_name
            group_total = np.sum(group_mask)
            group_kept = np.sum(mask & group_mask)
            group_filtered = group_total - group_kept

            stats["by_group"][str(group_name)] = {
                "n_total": int(group_total),
                "n_kept": int(group_kept),
                "n_filtered": int(group_filtered),
                "percent_kept": float(100 * group_kept / group_total) if group_total > 0 else 0.0,
            }

    return stats


def identify_problematic_cells(
    adata: anndata.AnnData, thresholds: Optional[Dict[str, tuple]] = None
) -> pd.DataFrame:
    """
    Identify cells that fail QC thresholds and report reasons.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    thresholds : dict, optional
        Dictionary mapping column names to (min, max) thresholds.

    Returns
    -------
    pd.DataFrame
        DataFrame with cell_id and failure reasons.
    """
    if thresholds is None:
        # Default thresholds
        thresholds = {
            "nCount_RNA": (10, None),
            "nFeature_RNA": (5, None),
        }

    problematic = []

    for idx, row in adata.obs.iterrows():
        reasons = []
        for col, (min_val, max_val) in thresholds.items():
            if col not in adata.obs.columns:
                continue

            value = row[col]

            if min_val is not None and value < min_val:
                reasons.append(f"{col} < {min_val} (value: {value})")

            if max_val is not None and value > max_val:
                reasons.append(f"{col} > {max_val} (value: {value})")

        if reasons:
            problematic.append(
                {
                    "cell_id": idx,
                    "n_failures": len(reasons),
                    "reasons": "; ".join(reasons),
                }
            )

    if problematic:
        return pd.DataFrame(problematic)
    else:
        return pd.DataFrame(columns=["cell_id", "n_failures", "reasons"])


def compare_pre_post_filtering(
    adata_pre: anndata.AnnData,
    adata_post: anndata.AnnData,
    metrics: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Compare statistics before and after filtering.

    Parameters
    ----------
    adata_pre : anndata.AnnData
        AnnData object before filtering.
    adata_post : anndata.AnnData
        AnnData object after filtering.
    metrics : list of str, optional
        Metrics to compare. If None, uses common numeric columns.

    Returns
    -------
    pd.DataFrame
        DataFrame comparing pre and post filtering statistics.
    """
    if metrics is None:
        # Use intersection of numeric columns
        numeric_cols_pre = set(adata_pre.obs.select_dtypes(include=[np.number]).columns)
        numeric_cols_post = set(adata_post.obs.select_dtypes(include=[np.number]).columns)
        metrics = list(numeric_cols_pre & numeric_cols_post)

    comparison = []

    for metric in metrics:
        pre_mean = adata_pre.obs[metric].mean()
        post_mean = adata_post.obs[metric].mean()

        comparison.append(
            {
                "metric": metric,
                "pre_mean": pre_mean,
                "post_mean": post_mean,
                "change": post_mean - pre_mean,
                "percent_change": 100 * (post_mean - pre_mean) / pre_mean if pre_mean != 0 else 0,
            }
        )

    return pd.DataFrame(comparison)
