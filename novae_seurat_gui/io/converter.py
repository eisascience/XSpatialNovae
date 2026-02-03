"""Converter utilities for spatial coordinates and metadata normalization."""

import logging
from typing import Optional, Tuple

import anndata
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def ensure_spatial_coords(
    adata: anndata.AnnData,
    x_col: Optional[str] = None,
    y_col: Optional[str] = None,
    spatial_key: str = "spatial",
    overwrite: bool = False,
) -> anndata.AnnData:
    """
    Ensure spatial coordinates are in adata.obsm[spatial_key].

    If x_col and y_col are provided, extract from adata.obs and store in obsm.
    Otherwise, check if spatial_key already exists in obsm.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    x_col : str, optional
        Column name for x coordinates in adata.obs.
    y_col : str, optional
        Column name for y coordinates in adata.obs.
    spatial_key : str
        Key to use in adata.obsm for spatial coordinates.
    overwrite : bool
        If True, overwrite existing spatial coordinates.

    Returns
    -------
    anndata.AnnData
        Modified AnnData object with spatial coordinates in obsm.
    """
    if spatial_key in adata.obsm and not overwrite:
        logger.info(f"Spatial coordinates already exist in adata.obsm['{spatial_key}']")
        return adata

    if x_col and y_col:
        if x_col not in adata.obs.columns:
            raise ValueError(f"Column '{x_col}' not found in adata.obs")
        if y_col not in adata.obs.columns:
            raise ValueError(f"Column '{y_col}' not found in adata.obs")

        spatial_coords = adata.obs[[x_col, y_col]].values.astype(float)
        adata.obsm[spatial_key] = spatial_coords
        logger.info(
            f"Created spatial coordinates in adata.obsm['{spatial_key}'] from {x_col}, {y_col}"
        )
    elif spatial_key not in adata.obsm:
        raise ValueError(
            f"No spatial coordinates found. Please provide x_col and y_col, "
            f"or ensure adata.obsm['{spatial_key}'] exists."
        )

    return adata


def normalize_metadata(
    adata: anndata.AnnData,
    cell_id_col: Optional[str] = None,
    sample_id_col: Optional[str] = None,
) -> anndata.AnnData:
    """
    Normalize metadata columns to standard names.

    Creates 'cell_id' and 'sample_id' columns in adata.obs if they don't exist.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    cell_id_col : str, optional
        Column name to use as cell_id. If None, uses the index.
    sample_id_col : str, optional
        Column name to use as sample_id.

    Returns
    -------
    anndata.AnnData
        Modified AnnData object with normalized metadata.
    """
    # Ensure cell_id
    if "cell_id" not in adata.obs.columns:
        if cell_id_col and cell_id_col in adata.obs.columns:
            adata.obs["cell_id"] = adata.obs[cell_id_col]
            logger.info(f"Created 'cell_id' from column '{cell_id_col}'")
        else:
            adata.obs["cell_id"] = adata.obs.index.astype(str)
            logger.info("Created 'cell_id' from adata.obs.index")

    # Ensure sample_id
    if "sample_id" not in adata.obs.columns:
        if sample_id_col and sample_id_col in adata.obs.columns:
            adata.obs["sample_id"] = adata.obs[sample_id_col].astype(str)
            logger.info(f"Created 'sample_id' from column '{sample_id_col}'")
        else:
            adata.obs["sample_id"] = "sample_0"
            logger.warning("No sample_id column found. Creating default 'sample_0'")

    return adata


def convert_units(
    coords: np.ndarray, from_units: str, to_units: str, conversion_factor: Optional[float] = None
) -> np.ndarray:
    """
    Convert spatial coordinates between units.

    Parameters
    ----------
    coords : np.ndarray
        Nx2 array of spatial coordinates.
    from_units : str
        Original units ('px' or 'mm').
    to_units : str
        Target units ('px' or 'mm').
    conversion_factor : float, optional
        Conversion factor. If None, uses default (1 mm = 1000 px).

    Returns
    -------
    np.ndarray
        Converted coordinates.
    """
    if from_units == to_units:
        return coords

    if conversion_factor is None:
        # Default: 1 mm = 1000 px
        if from_units == "px" and to_units == "mm":
            conversion_factor = 1.0 / 1000.0
        elif from_units == "mm" and to_units == "px":
            conversion_factor = 1000.0
        else:
            raise ValueError(f"Unsupported unit conversion: {from_units} to {to_units}")

    converted_coords = coords * conversion_factor
    logger.info(
        f"Converted coordinates from {from_units} to {to_units} "
        f"using factor {conversion_factor}"
    )

    return converted_coords


def add_cell_id_column(adata: anndata.AnnData, prefix: str = "cell_") -> anndata.AnnData:
    """
    Add a unique cell_id column to adata.obs if it doesn't exist.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    prefix : str
        Prefix for generated cell IDs.

    Returns
    -------
    anndata.AnnData
        Modified AnnData object.
    """
    if "cell_id" not in adata.obs.columns:
        adata.obs["cell_id"] = [f"{prefix}{i}" for i in range(adata.n_obs)]
        logger.info(f"Created unique cell_id column with prefix '{prefix}'")

    return adata


def standardize_column_names(adata: anndata.AnnData, mapping: dict) -> anndata.AnnData:
    """
    Rename columns in adata.obs according to a mapping dictionary.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    mapping : dict
        Dictionary mapping old column names to new column names.

    Returns
    -------
    anndata.AnnData
        Modified AnnData object.
    """
    rename_dict = {}
    for old_name, new_name in mapping.items():
        if old_name in adata.obs.columns and new_name not in adata.obs.columns:
            rename_dict[old_name] = new_name

    if rename_dict:
        adata.obs = adata.obs.rename(columns=rename_dict)
        logger.info(f"Renamed columns: {rename_dict}")

    return adata
