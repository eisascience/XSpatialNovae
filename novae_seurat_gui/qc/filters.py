"""QC filtering functions."""

import logging
from typing import Dict, List, Optional

import anndata
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def create_filter_mask(
    adata: anndata.AnnData, filter_criteria: Dict[str, tuple]
) -> np.ndarray:
    """
    Create a boolean mask for filtering cells based on metadata criteria.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    filter_criteria : dict
        Dictionary mapping column names to (min_value, max_value) tuples.
        Use None for unbounded. Example:
        {'nCount_RNA': (10, None), 'nFeature_RNA': (5, 1000)}

    Returns
    -------
    np.ndarray
        Boolean mask where True indicates cells that pass all filters.
    """
    mask = np.ones(adata.n_obs, dtype=bool)

    for col_name, (min_val, max_val) in filter_criteria.items():
        if col_name not in adata.obs.columns:
            logger.warning(f"Column '{col_name}' not found in adata.obs. Skipping.")
            continue

        col_data = adata.obs[col_name]

        if min_val is not None:
            mask &= col_data >= min_val

        if max_val is not None:
            mask &= col_data <= max_val

        n_filtered = np.sum(~mask)
        logger.info(
            f"Filter '{col_name}' [{min_val}, {max_val}]: "
            f"{n_filtered} cells filtered, {np.sum(mask)} remaining"
        )

    return mask


def apply_qc_filters(
    adata: anndata.AnnData,
    filter_criteria: Dict[str, tuple],
    inplace: bool = False,
    add_qc_column: bool = True,
) -> anndata.AnnData:
    """
    Apply QC filters to AnnData object.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    filter_criteria : dict
        Dictionary mapping column names to (min_value, max_value) tuples.
    inplace : bool
        If True, modify adata in place and return it.
        If False, return a filtered copy.
    add_qc_column : bool
        If True, add a 'qc_pass' column to adata.obs before filtering.

    Returns
    -------
    anndata.AnnData
        Filtered AnnData object.
    """
    mask = create_filter_mask(adata, filter_criteria)

    if add_qc_column:
        adata.obs["qc_pass"] = mask

    n_filtered = np.sum(~mask)
    n_kept = np.sum(mask)

    logger.info(
        f"QC filtering complete: {n_kept} cells kept, {n_filtered} cells filtered "
        f"({100 * n_filtered / adata.n_obs:.1f}%)"
    )

    if inplace:
        # Subset in place
        adata._inplace_subset_obs(mask)
        return adata
    else:
        # Return filtered copy
        return adata[mask, :].copy()


def filter_by_boolean_flags(
    adata: anndata.AnnData,
    flag_columns: List[str],
    require_all: bool = True,
    invert: bool = False,
) -> np.ndarray:
    """
    Filter cells based on boolean QC flag columns.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    flag_columns : list of str
        List of boolean column names in adata.obs.
    require_all : bool
        If True, require all flags to be True (AND logic).
        If False, require at least one flag to be True (OR logic).
    invert : bool
        If True, invert the logic (keep cells where flags are False).

    Returns
    -------
    np.ndarray
        Boolean mask indicating cells that pass the filter.
    """
    masks = []
    for col in flag_columns:
        if col not in adata.obs.columns:
            logger.warning(f"Flag column '{col}' not found. Skipping.")
            continue

        col_mask = adata.obs[col].astype(bool)
        masks.append(col_mask)

    if not masks:
        logger.warning("No valid flag columns found. Returning all True mask.")
        return np.ones(adata.n_obs, dtype=bool)

    if require_all:
        combined_mask = np.all(masks, axis=0)
    else:
        combined_mask = np.any(masks, axis=0)

    if invert:
        combined_mask = ~combined_mask

    return combined_mask


def filter_by_sample(
    adata: anndata.AnnData, sample_ids: List[str], sample_col: str = "sample_id"
) -> np.ndarray:
    """
    Create a mask to filter cells by sample ID.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    sample_ids : list of str
        List of sample IDs to keep.
    sample_col : str
        Column name containing sample IDs.

    Returns
    -------
    np.ndarray
        Boolean mask indicating cells that belong to the specified samples.
    """
    if sample_col not in adata.obs.columns:
        raise ValueError(f"Sample column '{sample_col}' not found in adata.obs")

    mask = adata.obs[sample_col].isin(sample_ids)
    n_kept = np.sum(mask)
    logger.info(
        f"Sample filter: keeping {n_kept} cells from samples {sample_ids}"
    )

    return mask.values


def filter_outliers_mad(
    adata: anndata.AnnData,
    column: str,
    n_mads: float = 5.0,
    only_upper: bool = False,
) -> np.ndarray:
    """
    Filter outliers based on median absolute deviation (MAD).

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    column : str
        Column name in adata.obs to filter.
    n_mads : float
        Number of MADs from the median to use as threshold.
    only_upper : bool
        If True, only filter upper outliers (values above median + n_mads*MAD).

    Returns
    -------
    np.ndarray
        Boolean mask indicating cells that are not outliers.
    """
    if column not in adata.obs.columns:
        raise ValueError(f"Column '{column}' not found in adata.obs")

    values = adata.obs[column]
    median = values.median()
    mad = (values - median).abs().median()

    if only_upper:
        upper_threshold = median + n_mads * mad
        mask = values <= upper_threshold
        logger.info(
            f"MAD filter ({column}): removed {np.sum(~mask)} upper outliers "
            f"(threshold: {upper_threshold:.2f})"
        )
    else:
        lower_threshold = median - n_mads * mad
        upper_threshold = median + n_mads * mad
        mask = (values >= lower_threshold) & (values <= upper_threshold)
        logger.info(
            f"MAD filter ({column}): removed {np.sum(~mask)} outliers "
            f"(thresholds: [{lower_threshold:.2f}, {upper_threshold:.2f}])"
        )

    return mask.values
