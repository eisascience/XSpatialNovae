"""Conversion utilities for .rds Seurat objects to H5AD format."""

import hashlib
import json
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple

import anndata
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def check_r_available() -> Tuple[bool, str]:
    """
    Check if R and Rscript are available in PATH.

    Returns
    -------
    tuple[bool, str]
        (is_available, message)
    """
    try:
        result = subprocess.run(
            ["Rscript", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0:
            version_info = result.stderr.strip() or result.stdout.strip()
            return True, f"R is available: {version_info}"
        else:
            return False, "Rscript command failed"
    except FileNotFoundError:
        return False, (
            "Rscript not found in PATH. Please install R (>=4.2) and ensure "
            "it is available in your system PATH."
        )
    except subprocess.TimeoutExpired:
        return False, "Rscript command timed out"
    except Exception as e:
        return False, f"Error checking R availability: {str(e)}"


def check_r_packages() -> Tuple[bool, str, list]:
    """
    Check if required R packages are installed.

    Returns
    -------
    tuple[bool, str, list]
        (all_available, message, missing_packages)
    """
    required_packages = ["Seurat", "SeuratDisk", "hdf5r", "optparse"]
    
    # Create R script to check packages
    check_script = """
    packages <- c({})
    missing <- c()
    for (pkg in packages) {{
        if (!requireNamespace(pkg, quietly = TRUE)) {{
            missing <- c(missing, pkg)
        }}
    }}
    if (length(missing) > 0) {{
        cat("MISSING:", paste(missing, collapse = ","))
    }} else {{
        cat("OK")
    }}
    """.format(", ".join([f'"{pkg}"' for pkg in required_packages]))
    
    try:
        result = subprocess.run(
            ["Rscript", "-e", check_script],
            capture_output=True,
            text=True,
            timeout=30,
        )
        
        output = result.stdout.strip()
        if output.startswith("MISSING:"):
            missing = output.replace("MISSING:", "").strip().split(",")
            msg = (
                f"Missing R packages: {', '.join(missing)}\n"
                f"Install with:\n"
                f"  R -e \"install.packages(c({', '.join([repr(p) for p in missing if p not in ['SeuratDisk']])}))\"\\n"
            )
            if "SeuratDisk" in missing:
                msg += "  R -e \"remotes::install_github('mojaveazure/seurat-disk')\"\n"
            return False, msg, missing
        elif output == "OK":
            return True, "All required R packages are installed", []
        else:
            return False, f"Unexpected output from R package check: {output}", []
            
    except Exception as e:
        return False, f"Error checking R packages: {str(e)}", []


def generate_cache_key(input_rds: str, params: Dict) -> str:
    """
    Generate a stable cache key from input file and parameters.

    Parameters
    ----------
    input_rds : str
        Path to input .rds file.
    params : dict
        Conversion parameters.

    Returns
    -------
    str
        Hash-based cache key.
    """
    # Get file size and modification time
    path = Path(input_rds)
    file_stats = {
        "size": path.stat().st_size,
        "mtime": path.stat().st_mtime,
    }
    
    # Combine with parameters
    cache_data = {
        "file_stats": file_stats,
        "params": params,
    }
    
    # Create hash
    cache_str = json.dumps(cache_data, sort_keys=True)
    cache_hash = hashlib.sha256(cache_str.encode()).hexdigest()[:16]
    
    return cache_hash


def convert_rds_to_h5ad(
    input_rds: str,
    output_dir: str,
    assay: str = "RNA",
    counts_slot: str = "counts",
    x_col: str = "auto",
    y_col: str = "auto",
    sample_id_col: str = "auto",
    cell_id_col: Optional[str] = None,
    keep_meta_regex: Optional[str] = None,
    overwrite: bool = False,
    use_cache: bool = True,
) -> Tuple[str, Dict]:
    """
    Convert a Seurat .rds object to H5AD format using R script.

    Parameters
    ----------
    input_rds : str
        Path to input .rds file containing Seurat object.
    output_dir : str
        Directory for output files.
    assay : str
        Assay name to extract (default: "RNA").
    counts_slot : str
        Slot to use for counts: "counts" or "data" (default: "counts").
    x_col : str
        Column name for x coordinates (default: "auto").
    y_col : str
        Column name for y coordinates (default: "auto").
    sample_id_col : str
        Column name for sample ID (default: "auto").
    cell_id_col : str, optional
        Column name for cell ID (default: use rownames).
    keep_meta_regex : str, optional
        Regex pattern to filter metadata columns.
    overwrite : bool
        Whether to overwrite existing files.
    use_cache : bool
        Whether to use cached conversion if available.

    Returns
    -------
    tuple[str, dict]
        (output_h5ad_path, conversion_info)

    Raises
    ------
    RuntimeError
        If R is not available or conversion fails.
    """
    # Check R availability
    r_available, r_msg = check_r_available()
    if not r_available:
        raise RuntimeError(f"R is not available: {r_msg}")
    
    logger.info(f"R check: {r_msg}")
    
    # Check R packages
    packages_ok, packages_msg, missing = check_r_packages()
    if not packages_ok:
        raise RuntimeError(f"Required R packages missing:\n{packages_msg}")
    
    logger.info("All required R packages are installed")
    
    # Prepare output paths
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Generate cache key
    params = {
        "assay": assay,
        "counts_slot": counts_slot,
        "x_col": x_col,
        "y_col": y_col,
        "sample_id_col": sample_id_col,
        "cell_id_col": cell_id_col,
        "keep_meta_regex": keep_meta_regex,
    }
    cache_key = generate_cache_key(input_rds, params)
    
    input_basename = Path(input_rds).stem
    output_h5ad = output_dir_path / f"{input_basename}_{cache_key}.h5ad"
    
    # Check cache
    if use_cache and output_h5ad.exists() and not overwrite:
        logger.info(f"Using cached conversion: {output_h5ad}")
        
        # Load and return info
        try:
            adata = anndata.read_h5ad(output_h5ad)
            conversion_info = {
                "input_rds": str(input_rds),
                "output_h5ad": str(output_h5ad),
                "cached": True,
                "n_cells": adata.n_obs,
                "n_features": adata.n_vars,
                "params": params,
            }
            return str(output_h5ad), conversion_info
        except Exception as e:
            logger.warning(f"Failed to load cached file: {e}. Re-running conversion.")
    
    # Build R script command
    script_path = Path(__file__).parent.parent.parent / "scripts" / "convert_seurat_rds_to_h5ad.R"
    
    if not script_path.exists():
        raise RuntimeError(f"Conversion script not found: {script_path}")
    
    cmd = [
        "Rscript",
        str(script_path),
        "--input_rds", str(input_rds),
        "--output_h5ad", str(output_h5ad),
        "--assay", assay,
        "--counts_slot", counts_slot,
        "--x_col", x_col,
        "--y_col", y_col,
        "--sample_id_col", sample_id_col,
    ]
    
    if cell_id_col:
        cmd.extend(["--cell_id_col", cell_id_col])
    
    if keep_meta_regex:
        cmd.extend(["--keep_meta_regex", keep_meta_regex])
    
    if overwrite:
        cmd.append("--overwrite=TRUE")
    
    # Run conversion
    logger.info(f"Running R conversion: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 minutes timeout
        )
        
        # Log output
        if result.stdout:
            logger.info(f"R script output:\n{result.stdout}")
        
        if result.returncode != 0:
            error_msg = result.stderr or "Unknown error"
            raise RuntimeError(
                f"R conversion failed (exit code {result.returncode}):\n{error_msg}\n\n"
                f"stdout:\n{result.stdout}"
            )
        
        # Check if output file was created
        if not output_h5ad.exists():
            raise RuntimeError(
                f"R conversion completed but output file not found: {output_h5ad}\n"
                f"stdout:\n{result.stdout}\n"
                f"stderr:\n{result.stderr}"
            )
        
        # Load to get info
        adata = anndata.read_h5ad(output_h5ad)
        
        conversion_info = {
            "input_rds": str(input_rds),
            "output_h5ad": str(output_h5ad),
            "cached": False,
            "n_cells": adata.n_obs,
            "n_features": adata.n_vars,
            "params": params,
            "r_stdout": result.stdout,
        }
        
        logger.info(f"Conversion successful: {output_h5ad}")
        logger.info(f"Cells: {adata.n_obs}, Features: {adata.n_vars}")
        
        return str(output_h5ad), conversion_info
        
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            "R conversion timed out after 10 minutes. "
            "The input file may be too large or there may be an issue with the Seurat object."
        )
    except Exception as e:
        raise RuntimeError(f"Error during R conversion: {str(e)}")


def ensure_anndata_schema(
    h5ad_path: str,
    x_col: Optional[str] = None,
    y_col: Optional[str] = None,
    sample_id_col: Optional[str] = None,
    overwrite: bool = True,
) -> anndata.AnnData:
    """
    Ensure AnnData has required schema after conversion from Seurat.

    This function:
    - Creates obsm["spatial"] from obs columns if missing
    - Ensures obs["cell_id"] exists
    - Ensures obs["sample_id"] exists
    - Validates counts data

    Parameters
    ----------
    h5ad_path : str
        Path to H5AD file.
    x_col : str, optional
        Column name for x coordinates in obs.
    y_col : str, optional
        Column name for y coordinates in obs.
    sample_id_col : str, optional
        Column name for sample ID in obs.
    overwrite : bool
        Whether to save changes back to the file.

    Returns
    -------
    anndata.AnnData
        Modified AnnData object.

    Raises
    ------
    ValueError
        If required data is missing and cannot be created.
    """
    logger.info(f"Ensuring AnnData schema for {h5ad_path}")
    
    adata = anndata.read_h5ad(h5ad_path)
    modified = False
    
    # 1. Ensure obsm["spatial"] exists
    if "spatial" not in adata.obsm:
        logger.info("obsm['spatial'] not found, attempting to create from obs columns")
        
        # Auto-detect x/y columns if not provided
        if x_col is None or y_col is None:
            x_candidates = ["x_slide_mm", "x_location", "x_centroid", "x_FOV_px", "x", "X"]
            y_candidates = ["y_slide_mm", "y_location", "y_centroid", "y_FOV_px", "y", "Y"]
            
            for x_cand in x_candidates:
                if x_cand in adata.obs.columns:
                    x_col = x_cand
                    logger.info(f"Auto-detected x column: {x_col}")
                    break
            
            for y_cand in y_candidates:
                if y_cand in adata.obs.columns:
                    y_col = y_cand
                    logger.info(f"Auto-detected y column: {y_col}")
                    break
        
        if x_col and y_col and x_col in adata.obs.columns and y_col in adata.obs.columns:
            spatial_coords = adata.obs[[x_col, y_col]].values.astype(float)
            adata.obsm["spatial"] = spatial_coords
            logger.info(f"Created obsm['spatial'] from {x_col}, {y_col}")
            modified = True
        else:
            available_cols = list(adata.obs.columns)
            raise ValueError(
                f"Cannot create obsm['spatial']: x_col='{x_col}', y_col='{y_col}' not found.\n"
                f"Available columns: {available_cols}"
            )
    else:
        logger.info("obsm['spatial'] already exists")
    
    # 2. Ensure obs["cell_id"] exists
    if "cell_id" not in adata.obs.columns:
        adata.obs["cell_id"] = adata.obs.index.astype(str)
        logger.info("Created obs['cell_id'] from index")
        modified = True
    else:
        logger.info("obs['cell_id'] already exists")
    
    # 3. Ensure obs["sample_id"] exists
    if "sample_id" not in adata.obs.columns:
        if sample_id_col and sample_id_col in adata.obs.columns:
            adata.obs["sample_id"] = adata.obs[sample_id_col].astype(str)
            logger.info(f"Created obs['sample_id'] from {sample_id_col}")
            modified = True
        else:
            # Check common column names
            sample_candidates = ["SampleID", "Tissue", "tissue", "sample"]
            found = False
            for candidate in sample_candidates:
                if candidate in adata.obs.columns:
                    adata.obs["sample_id"] = adata.obs[candidate].astype(str)
                    logger.info(f"Created obs['sample_id'] from {candidate}")
                    modified = True
                    found = True
                    break
            
            if not found:
                # Default to single sample
                adata.obs["sample_id"] = "sample_0"
                logger.warning("No sample_id column found, created default 'sample_0'")
                modified = True
    else:
        logger.info("obs['sample_id'] already exists")
    
    # 4. Check counts availability
    has_counts_layer = "counts" in adata.layers if adata.layers else False
    has_x = adata.X is not None
    
    if not has_counts_layer and not has_x:
        raise ValueError("No counts data found in adata.X or adata.layers['counts']")
    
    if has_counts_layer:
        logger.info("Counts available in adata.layers['counts']")
    elif has_x:
        logger.info("Counts available in adata.X")
    
    # Save if modified
    if modified and overwrite:
        logger.info(f"Saving modified AnnData to {h5ad_path}")
        adata.write_h5ad(h5ad_path)
    
    return adata


def convert_rds_to_h5ad_with_validation(
    input_rds: str,
    output_dir: str,
    assay: str = "RNA",
    counts_slot: str = "counts",
    x_col: str = "auto",
    y_col: str = "auto",
    sample_id_col: str = "auto",
    cell_id_col: Optional[str] = None,
    keep_meta_regex: Optional[str] = None,
    overwrite: bool = False,
    use_cache: bool = True,
) -> Tuple[str, anndata.AnnData, Dict]:
    """
    Convert .rds to H5AD with full validation and schema enforcement.

    This is the main entry point that:
    1. Runs R conversion script
    2. Ensures AnnData schema is correct
    3. Returns validated AnnData object

    Parameters
    ----------
    See convert_rds_to_h5ad for parameter descriptions.

    Returns
    -------
    tuple[str, anndata.AnnData, dict]
        (output_path, adata, conversion_info)
    """
    # Run R conversion
    output_path, conversion_info = convert_rds_to_h5ad(
        input_rds=input_rds,
        output_dir=output_dir,
        assay=assay,
        counts_slot=counts_slot,
        x_col=x_col,
        y_col=y_col,
        sample_id_col=sample_id_col,
        cell_id_col=cell_id_col,
        keep_meta_regex=keep_meta_regex,
        overwrite=overwrite,
        use_cache=use_cache,
    )
    
    # Post-process to ensure schema
    # Extract actual column names from params (handle "auto" case)
    x_col_actual = x_col if x_col != "auto" else None
    y_col_actual = y_col if y_col != "auto" else None
    sample_id_col_actual = sample_id_col if sample_id_col != "auto" else None
    
    adata = ensure_anndata_schema(
        output_path,
        x_col=x_col_actual,
        y_col=y_col_actual,
        sample_id_col=sample_id_col_actual,
        overwrite=True,
    )
    
    # Add schema validation info to conversion_info
    conversion_info["schema_validated"] = True
    conversion_info["has_spatial"] = "spatial" in adata.obsm
    conversion_info["has_cell_id"] = "cell_id" in adata.obs.columns
    conversion_info["has_sample_id"] = "sample_id" in adata.obs.columns
    
    return output_path, adata, conversion_info
