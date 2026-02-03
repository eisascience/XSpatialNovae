"""Spatial graph diagnostics and statistics."""

import logging
from typing import Dict, Optional

import anndata
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def graph_diagnostics(adata: anndata.AnnData, connectivities_key: str = "spatial_connectivities") -> Dict:
    """
    Compute diagnostic statistics for spatial neighbor graph.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with neighbor graph.
    connectivities_key : str
        Key in adata.obsp containing connectivity matrix.

    Returns
    -------
    dict
        Dictionary with diagnostic statistics.
    """
    if connectivities_key not in adata.obsp:
        raise ValueError(f"Connectivity key '{connectivities_key}' not found in adata.obsp")

    connectivities = adata.obsp[connectivities_key]

    # Degree (number of neighbors per cell)
    degree = np.array(connectivities.sum(axis=1)).flatten()

    # Connected components
    from scipy.sparse.csgraph import connected_components

    n_components, labels = connected_components(
        connectivities, directed=False, return_labels=True
    )

    # Component sizes
    component_sizes = pd.Series(labels).value_counts().sort_values(ascending=False)

    diagnostics = {
        "n_cells": adata.n_obs,
        "n_edges": int(connectivities.nnz // 2),  # Divide by 2 for undirected
        "degree": {
            "mean": float(degree.mean()),
            "median": float(np.median(degree)),
            "min": int(degree.min()),
            "max": int(degree.max()),
            "std": float(degree.std()),
        },
        "connected_components": {
            "n_components": int(n_components),
            "largest_component_size": int(component_sizes.iloc[0]) if len(component_sizes) > 0 else 0,
            "isolated_cells": int(np.sum(degree == 0)),
        },
        "sparsity": float(1 - (connectivities.nnz / (adata.n_obs ** 2))),
    }

    logger.info(
        f"Graph diagnostics: {diagnostics['n_edges']} edges, "
        f"{diagnostics['degree']['mean']:.1f} avg degree, "
        f"{diagnostics['connected_components']['n_components']} components"
    )

    return diagnostics


def compute_degree_distribution(
    adata: anndata.AnnData, connectivities_key: str = "spatial_connectivities"
) -> pd.Series:
    """
    Compute degree distribution for spatial neighbor graph.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with neighbor graph.
    connectivities_key : str
        Key in adata.obsp containing connectivity matrix.

    Returns
    -------
    pd.Series
        Series with degree as index and count as values.
    """
    if connectivities_key not in adata.obsp:
        raise ValueError(f"Connectivity key '{connectivities_key}' not found in adata.obsp")

    connectivities = adata.obsp[connectivities_key]
    degree = np.array(connectivities.sum(axis=1)).flatten()

    degree_dist = pd.Series(degree).value_counts().sort_index()

    return degree_dist


def plot_degree_distribution(
    adata: anndata.AnnData,
    connectivities_key: str = "spatial_connectivities",
    figsize: tuple = (8, 5),
):
    """
    Plot degree distribution histogram.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with neighbor graph.
    connectivities_key : str
        Key in adata.obsp containing connectivity matrix.
    figsize : tuple
        Figure size.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object.
    """
    import matplotlib.pyplot as plt

    degree_dist = compute_degree_distribution(adata, connectivities_key)

    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(degree_dist.index, degree_dist.values)
    ax.set_xlabel("Number of Neighbors")
    ax.set_ylabel("Number of Cells")
    ax.set_title("Spatial Neighbor Graph Degree Distribution")
    ax.grid(True, alpha=0.3)

    return fig


def identify_isolated_cells(
    adata: anndata.AnnData, connectivities_key: str = "spatial_connectivities"
) -> np.ndarray:
    """
    Identify cells with no spatial neighbors.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with neighbor graph.
    connectivities_key : str
        Key in adata.obsp containing connectivity matrix.

    Returns
    -------
    np.ndarray
        Array of cell indices (positions in adata.obs) that are isolated.
    """
    if connectivities_key not in adata.obsp:
        raise ValueError(f"Connectivity key '{connectivities_key}' not found in adata.obsp")

    connectivities = adata.obsp[connectivities_key]
    degree = np.array(connectivities.sum(axis=1)).flatten()

    isolated_indices = np.where(degree == 0)[0]

    logger.info(f"Found {len(isolated_indices)} isolated cells (no neighbors)")

    return isolated_indices


def compute_local_density(
    adata: anndata.AnnData,
    spatial_key: str = "spatial",
    radius: Optional[float] = None,
) -> np.ndarray:
    """
    Compute local cell density for each cell.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    spatial_key : str
        Key in adata.obsm containing spatial coordinates.
    radius : float, optional
        Radius for density computation. If None, uses median nearest neighbor distance * 2.

    Returns
    -------
    np.ndarray
        Array of local density values (cells per unit area).
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")

    coords = adata.obsm[spatial_key]

    from scipy.spatial import cKDTree

    tree = cKDTree(coords)

    if radius is None:
        # Estimate radius from median NN distance
        distances, _ = tree.query(coords, k=2)  # k=2 to get first NN (not self)
        median_dist = np.median(distances[:, 1])
        radius = median_dist * 2
        logger.info(f"Auto-selected density radius: {radius:.2f}")

    # Count neighbors within radius
    neighbor_counts = tree.query_ball_point(coords, radius, return_length=True)

    # Compute density (cells per unit area)
    area = np.pi * (radius ** 2)
    density = (neighbor_counts - 1) / area  # Subtract self

    return density


def compute_edge_lengths(
    adata: anndata.AnnData,
    spatial_key: str = "spatial",
    connectivities_key: str = "spatial_connectivities",
    sample_size: Optional[int] = None,
) -> np.ndarray:
    """
    Compute edge lengths in spatial neighbor graph.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    spatial_key : str
        Key in adata.obsm containing spatial coordinates.
    connectivities_key : str
        Key in adata.obsp containing connectivity matrix.
    sample_size : int, optional
        Sample a subset of edges for efficiency. If None, computes all edges.

    Returns
    -------
    np.ndarray
        Array of edge lengths.
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")
    if connectivities_key not in adata.obsp:
        raise ValueError(f"Connectivity key '{connectivities_key}' not found in adata.obsp")

    coords = adata.obsm[spatial_key]
    connectivities = adata.obsp[connectivities_key]

    # Get edges
    edges = np.array(connectivities.nonzero()).T

    # Remove duplicate edges (keep only i < j for undirected graph)
    edges = edges[edges[:, 0] < edges[:, 1]]

    if sample_size is not None and len(edges) > sample_size:
        sample_indices = np.random.choice(len(edges), size=sample_size, replace=False)
        edges = edges[sample_indices]

    # Compute edge lengths
    edge_lengths = np.linalg.norm(
        coords[edges[:, 0]] - coords[edges[:, 1]], axis=1
    )

    return edge_lengths


def get_neighbor_edges_for_visualization(
    adata: anndata.AnnData,
    spatial_key: str = "spatial",
    connectivities_key: str = "spatial_connectivities",
    max_edges: int = 10000,
) -> tuple:
    """
    Get edge coordinates for visualization, downsampled if necessary.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    spatial_key : str
        Key in adata.obsm containing spatial coordinates.
    connectivities_key : str
        Key in adata.obsp containing connectivity matrix.
    max_edges : int
        Maximum number of edges to return.

    Returns
    -------
    tuple
        (edge_coords_x, edge_coords_y) where each is a list of (start, end, None) for plotting.
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")
    if connectivities_key not in adata.obsp:
        raise ValueError(f"Connectivity key '{connectivities_key}' not found in adata.obsp")

    coords = adata.obsm[spatial_key]
    connectivities = adata.obsp[connectivities_key]

    # Get edges
    edges = np.array(connectivities.nonzero()).T
    edges = edges[edges[:, 0] < edges[:, 1]]  # Remove duplicates

    # Downsample if needed
    if len(edges) > max_edges:
        sample_indices = np.random.choice(len(edges), size=max_edges, replace=False)
        edges = edges[sample_indices]
        logger.info(f"Downsampled edges from {connectivities.nnz // 2} to {max_edges} for visualization")

    # Create line segments
    edge_x = []
    edge_y = []

    for i, j in edges:
        edge_x.extend([coords[i, 0], coords[j, 0], None])
        edge_y.extend([coords[i, 1], coords[j, 1], None])

    return edge_x, edge_y
