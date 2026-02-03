"""Manifest creation for documenting run parameters and metadata."""

import hashlib
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import anndata

logger = logging.getLogger(__name__)


def compute_file_hash(file_path: str, algorithm: str = "sha256") -> str:
    """
    Compute hash of a file.

    Parameters
    ----------
    file_path : str
        Path to file.
    algorithm : str
        Hash algorithm ('md5', 'sha256').

    Returns
    -------
    str
        Hex digest of file hash.
    """
    hash_func = hashlib.new(algorithm)

    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_func.update(chunk)

    return hash_func.hexdigest()


def create_manifest(
    adata: anndata.AnnData,
    input_files: list,
    parameters: Optional[Dict[str, Any]] = None,
    qc_filters: Optional[Dict[str, Any]] = None,
    n_cells_pre_qc: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Create a manifest documenting the analysis run.

    Parameters
    ----------
    adata : anndata.AnnData
        Processed AnnData object.
    input_files : list
        List of input file paths.
    parameters : dict, optional
        Analysis parameters (Novae, preprocessing, etc.).
    qc_filters : dict, optional
        QC filter criteria applied.
    n_cells_pre_qc : int, optional
        Number of cells before QC filtering.

    Returns
    -------
    dict
        Manifest dictionary.
    """
    manifest = {
        "timestamp": datetime.now().isoformat(),
        "version": "0.1.0",
        "input": {
            "files": [],
            "n_cells_total": n_cells_pre_qc or adata.n_obs,
            "n_features": adata.n_vars,
        },
        "processing": {
            "qc_filters": qc_filters or {},
            "n_cells_kept": adata.n_obs,
            "n_cells_filtered": (n_cells_pre_qc - adata.n_obs) if n_cells_pre_qc else 0,
            "percent_kept": (
                100 * adata.n_obs / n_cells_pre_qc if n_cells_pre_qc else 100.0
            ),
        },
        "parameters": parameters or {},
        "output": {
            "n_domains": int(adata.obs["domain"].nunique()) if "domain" in adata.obs else None,
            "embedding_dimensions": (
                adata.obsm["X_novae"].shape[1] if "X_novae" in adata.obsm else None
            ),
        },
        "data_schema": {
            "cell_id_column": "cell_id",
            "sample_id_column": "sample_id",
            "spatial_key": "spatial" if "spatial" in adata.obsm else None,
            "coordinate_units": parameters.get("coordinate_units", "unknown") if parameters else "unknown",
        },
    }

    # Add input file information
    for file_path in input_files:
        if Path(file_path).exists():
            file_info = {
                "path": str(file_path),
                "name": Path(file_path).name,
                "size_bytes": Path(file_path).stat().st_size,
                "sha256": compute_file_hash(file_path, "sha256"),
            }
            manifest["input"]["files"].append(file_info)

    # Add sample information
    if "sample_id" in adata.obs:
        manifest["samples"] = {
            "n_samples": int(adata.obs["sample_id"].nunique()),
            "sample_ids": adata.obs["sample_id"].unique().tolist(),
            "cells_per_sample": adata.obs["sample_id"].value_counts().to_dict(),
        }

    # Add software versions
    import sys

    try:
        import scanpy as sc
        import anndata as ad

        manifest["software"] = {
            "python_version": sys.version,
            "novae_seurat_gui_version": "0.1.0",
            "scanpy_version": sc.__version__,
            "anndata_version": ad.__version__,
        }
    except Exception as e:
        logger.warning(f"Could not retrieve software versions: {e}")

    return manifest


def save_manifest(manifest: Dict[str, Any], output_file: str) -> None:
    """
    Save manifest to JSON file.

    Parameters
    ----------
    manifest : dict
        Manifest dictionary.
    output_file : str
        Output JSON file path.
    """
    logger.info(f"Saving manifest to {output_file}")

    with open(output_file, "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    logger.info("Manifest saved")


def add_qc_summary_to_manifest(
    manifest: Dict[str, Any], qc_summary: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Add QC summary statistics to manifest.

    Parameters
    ----------
    manifest : dict
        Existing manifest dictionary.
    qc_summary : dict
        QC summary from qc.summaries.compute_qc_summary.

    Returns
    -------
    dict
        Updated manifest.
    """
    manifest["qc_summary"] = qc_summary
    return manifest


def add_spatial_summary_to_manifest(
    manifest: Dict[str, Any], spatial_diagnostics: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Add spatial graph diagnostics to manifest.

    Parameters
    ----------
    manifest : dict
        Existing manifest dictionary.
    spatial_diagnostics : dict
        Spatial diagnostics from spatial.diagnostics.graph_diagnostics.

    Returns
    -------
    dict
        Updated manifest.
    """
    manifest["spatial_graph"] = spatial_diagnostics
    return manifest


def validate_manifest(manifest: Dict[str, Any]) -> tuple[bool, list]:
    """
    Validate manifest structure.

    Parameters
    ----------
    manifest : dict
        Manifest dictionary to validate.

    Returns
    -------
    tuple of (bool, list)
        (is_valid, list of error messages)
    """
    errors = []

    required_keys = ["timestamp", "version", "input", "processing", "parameters"]
    for key in required_keys:
        if key not in manifest:
            errors.append(f"Missing required key: {key}")

    # Validate input section
    if "input" in manifest:
        if "files" not in manifest["input"]:
            errors.append("Missing 'files' in input section")

    # Validate processing section
    if "processing" in manifest:
        if "n_cells_kept" not in manifest["processing"]:
            errors.append("Missing 'n_cells_kept' in processing section")

    is_valid = len(errors) == 0

    return is_valid, errors
