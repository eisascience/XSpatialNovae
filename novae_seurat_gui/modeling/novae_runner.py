"""Novae model runner and preprocessing."""

import logging
from typing import Optional

import anndata
import numpy as np
import scanpy as sc

logger = logging.getLogger(__name__)


def preprocess_for_novae(
    adata: anndata.AnnData,
    target_sum: Optional[float] = 1e4,
    n_top_genes: int = 2000,
    n_comps: int = 50,
    use_layer: Optional[str] = None,
    skip_normalize: bool = False,
    skip_log: bool = False,
    random_state: int = 42,
) -> anndata.AnnData:
    """
    Preprocess AnnData for Novae: normalize, log-transform, feature selection, PCA.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    target_sum : float, optional
        Target sum for normalization. If None, skips normalization.
    n_top_genes : int
        Number of highly variable genes to select.
    n_comps : int
        Number of principal components to compute.
    use_layer : str, optional
        Layer to use as input. If None, uses adata.X.
    skip_normalize : bool
        If True, skip normalization (for proteomics with quantile scaling).
    skip_log : bool
        If True, skip log transformation.
    random_state : int
        Random seed for reproducibility.

    Returns
    -------
    anndata.AnnData
        Preprocessed AnnData object with PCA in obsm['X_pca'].
    """
    logger.info("Starting preprocessing for Novae")

    # Make a copy to avoid modifying original
    adata = adata.copy()

    # Use specified layer
    if use_layer is not None:
        if use_layer not in adata.layers:
            raise ValueError(f"Layer '{use_layer}' not found in adata.layers")
        adata.X = adata.layers[use_layer].copy()
        logger.info(f"Using layer '{use_layer}' as input")

    # Store raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # Normalize
    if not skip_normalize and target_sum is not None:
        logger.info(f"Normalizing to target sum {target_sum}")
        sc.pp.normalize_total(adata, target_sum=target_sum)

    # Log transform
    if not skip_log:
        logger.info("Log-transforming (log1p)")
        sc.pp.log1p(adata)

    # Store normalized data
    adata.layers["normalized"] = adata.X.copy()

    # Highly variable genes
    logger.info(f"Selecting {n_top_genes} highly variable genes")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor="seurat_v3", layer="counts"
    )

    # Scale
    logger.info("Scaling data")
    sc.pp.scale(adata, max_value=10)

    # PCA
    logger.info(f"Computing {n_comps} principal components")
    sc.tl.pca(adata, n_comps=n_comps, random_state=random_state)

    logger.info("Preprocessing complete")

    return adata


def run_novae_zeroshot(
    adata: anndata.AnnData,
    model_name: str = "MICS-Lab/novae-human-0",
    n_domains: int = 10,
    use_pca: bool = True,
    n_comps: int = 50,
    random_state: int = 42,
) -> anndata.AnnData:
    """
    Run Novae in zero-shot mode using a pretrained model.

    NOTE: This is a placeholder implementation. The actual Novae model integration
    would require the novae package to be installed and proper model loading.

    Parameters
    ----------
    adata : anndata.AnnData
        Preprocessed AnnData object.
    model_name : str
        Pretrained model identifier.
    n_domains : int
        Number of spatial domains to assign.
    use_pca : bool
        If True, use PCA representation as input.
    n_comps : int
        Number of PCs to use if use_pca=True.
    random_state : int
        Random seed.

    Returns
    -------
    anndata.AnnData
        AnnData with Novae embeddings in obsm['X_novae'] and domains in obs['domain'].
    """
    logger.warning(
        "Novae integration is a placeholder. Install 'novae' package for full functionality."
    )

    # Placeholder: Use UMAP + Leiden clustering as a proxy
    logger.info(f"Running Novae zero-shot with model {model_name}")

    # Ensure PCA is computed
    if "X_pca" not in adata.obsm:
        logger.info("Computing PCA as input representation")
        sc.tl.pca(adata, n_comps=n_comps, random_state=random_state)

    # Compute neighborhood graph on PCA
    logger.info("Computing neighborhood graph")
    sc.pp.neighbors(adata, n_pcs=n_comps, random_state=random_state)

    # Compute UMAP as embedding (proxy for Novae embeddings)
    logger.info("Computing embeddings (UMAP as proxy)")
    sc.tl.umap(adata, random_state=random_state)

    # Store as "novae" embeddings
    adata.obsm["X_novae"] = adata.obsm["X_umap"].copy()

    # Leiden clustering for domain assignment
    logger.info(f"Assigning {n_domains} spatial domains (Leiden clustering)")
    sc.tl.leiden(adata, resolution=1.0, n_iterations=2, random_state=random_state)

    # Rename to domain
    adata.obs["domain"] = adata.obs["leiden"].astype(str)

    # Add hierarchical domain levels if needed
    adata.obs["domain_level_0"] = adata.obs["domain"]

    # Add metadata
    adata.uns["novae_params"] = {
        "model_name": model_name,
        "n_domains": n_domains,
        "use_pca": use_pca,
        "n_comps": n_comps,
        "random_state": random_state,
    }

    logger.info(
        f"Novae zero-shot complete. Assigned {adata.obs['domain'].nunique()} domains"
    )

    return adata


def run_novae_training(
    adata: anndata.AnnData,
    n_epochs: int = 100,
    batch_size: int = 256,
    learning_rate: float = 1e-3,
    n_domains: int = 10,
    random_state: int = 42,
) -> anndata.AnnData:
    """
    Train Novae model from scratch (for proteomics/PhenoCycler data).

    NOTE: This is a placeholder implementation.

    Parameters
    ----------
    adata : anndata.AnnData
        Preprocessed AnnData object.
    n_epochs : int
        Number of training epochs.
    batch_size : int
        Batch size for training.
    learning_rate : float
        Learning rate.
    n_domains : int
        Number of spatial domains.
    random_state : int
        Random seed.

    Returns
    -------
    anndata.AnnData
        AnnData with trained Novae model outputs.
    """
    logger.warning(
        "Novae training from scratch is a placeholder. "
        "Install 'novae' package for full functionality."
    )

    logger.info(
        f"Training Novae from scratch: {n_epochs} epochs, "
        f"batch_size={batch_size}, lr={learning_rate}"
    )

    # Placeholder: Use similar approach as zero-shot
    adata = run_novae_zeroshot(
        adata,
        model_name="trained_from_scratch",
        n_domains=n_domains,
        random_state=random_state,
    )

    adata.uns["novae_params"]["training"] = {
        "n_epochs": n_epochs,
        "batch_size": batch_size,
        "learning_rate": learning_rate,
    }

    logger.info("Novae training complete")

    return adata


def apply_quantile_scaling(adata: anndata.AnnData, axis: int = 0) -> anndata.AnnData:
    """
    Apply quantile (rank) scaling for proteomics data.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    axis : int
        Axis to scale (0=features, 1=cells).

    Returns
    -------
    anndata.AnnData
        AnnData with quantile-scaled data.
    """
    from scipy.stats import rankdata

    logger.info("Applying quantile scaling")

    if axis == 0:
        # Scale features (proteins)
        for i in range(adata.n_vars):
            adata.X[:, i] = rankdata(adata.X[:, i]) / adata.n_obs
    else:
        # Scale cells
        for i in range(adata.n_obs):
            adata.X[i, :] = rankdata(adata.X[i, :]) / adata.n_vars

    logger.info("Quantile scaling complete")

    return adata
