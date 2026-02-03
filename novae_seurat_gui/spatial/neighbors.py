"""Spatial neighbor graph construction."""

import logging
from typing import Literal, Optional

import anndata
import numpy as np
from scipy.spatial import cKDTree
from sklearn.neighbors import NearestNeighbors, radius_neighbors_graph

logger = logging.getLogger(__name__)


def compute_neighbors(
    adata: anndata.AnnData,
    spatial_key: str = "spatial",
    method: Literal["radius", "knn"] = "radius",
    radius: Optional[float] = None,
    n_neighbors: Optional[int] = None,
    coord_type: str = "generic",
) -> anndata.AnnData:
    """
    Compute spatial neighbors and add to adata.obsp and adata.uns.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object with spatial coordinates in obsm[spatial_key].
    spatial_key : str
        Key in adata.obsm containing spatial coordinates.
    method : {'radius', 'knn'}
        Method for computing neighbors.
    radius : float, optional
        Radius for radius-based neighbors (required if method='radius').
    n_neighbors : int, optional
        Number of neighbors for KNN (required if method='knn').
    coord_type : str
        Type of coordinates ('generic', 'um', 'mm', 'px').

    Returns
    -------
    anndata.AnnData
        Modified AnnData object with neighbor graph in obsp['spatial_connectivities']
        and 'spatial_distances', and metadata in uns['spatial_neighbors'].
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")

    coords = adata.obsm[spatial_key]

    if method == "radius":
        if radius is None:
            raise ValueError("Must specify 'radius' when method='radius'")

        logger.info(f"Computing radius-based neighbors with radius={radius}")
        connectivities = radius_neighbors_graph(
            coords, radius=radius, mode="connectivity", include_self=False
        )
        distances = radius_neighbors_graph(
            coords, radius=radius, mode="distance", include_self=False
        )

    elif method == "knn":
        if n_neighbors is None:
            raise ValueError("Must specify 'n_neighbors' when method='knn'")

        logger.info(f"Computing KNN with n_neighbors={n_neighbors}")
        nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm="auto")
        nbrs.fit(coords)
        distances, indices = nbrs.kneighbors(coords)

        # Convert to sparse matrix
        from scipy.sparse import csr_matrix

        n_cells = coords.shape[0]
        connectivities = csr_matrix((n_cells, n_cells))
        distances_sparse = csr_matrix((n_cells, n_cells))

        row_indices = np.repeat(np.arange(n_cells), n_neighbors)
        col_indices = indices.flatten()
        data_ones = np.ones(len(row_indices))
        data_dists = distances.flatten()

        connectivities = csr_matrix(
            (data_ones, (row_indices, col_indices)), shape=(n_cells, n_cells)
        )
        distances_sparse = csr_matrix(
            (data_dists, (row_indices, col_indices)), shape=(n_cells, n_cells)
        )

        connectivities = connectivities.maximum(connectivities.T)  # Make symmetric
        distances_sparse = distances_sparse.maximum(distances_sparse.T)

    else:
        raise ValueError(f"Unknown method: {method}")

    # Store in AnnData
    adata.obsp["spatial_connectivities"] = connectivities
    adata.obsp["spatial_distances"] = distances_sparse

    # Store metadata
    adata.uns["spatial_neighbors"] = {
        "connectivities_key": "spatial_connectivities",
        "distances_key": "spatial_distances",
        "params": {
            "method": method,
            "radius": radius,
            "n_neighbors": n_neighbors,
            "coord_type": coord_type,
            "spatial_key": spatial_key,
        },
    }

    n_edges = connectivities.nnz // 2  # Divide by 2 for undirected graph
    logger.info(
        f"Spatial neighbor graph constructed: {n_edges} edges, "
        f"avg {connectivities.nnz / adata.n_obs:.1f} neighbors per cell"
    )

    return adata


def build_spatial_graph(
    adata: anndata.AnnData,
    spatial_key: str = "spatial",
    radius: Optional[float] = None,
    n_neighbors: Optional[int] = None,
    coord_type: str = "generic",
) -> anndata.AnnData:
    """
    Build spatial neighbor graph using automatic method selection.

    If radius is provided, uses radius-based neighbors.
    If n_neighbors is provided, uses KNN.
    If both are provided, uses radius.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    spatial_key : str
        Key in adata.obsm containing spatial coordinates.
    radius : float, optional
        Radius for radius-based neighbors.
    n_neighbors : int, optional
        Number of neighbors for KNN.
    coord_type : str
        Type of coordinates.

    Returns
    -------
    anndata.AnnData
        Modified AnnData object with spatial neighbor graph.
    """
    if radius is not None:
        return compute_neighbors(
            adata,
            spatial_key=spatial_key,
            method="radius",
            radius=radius,
            coord_type=coord_type,
        )
    elif n_neighbors is not None:
        return compute_neighbors(
            adata,
            spatial_key=spatial_key,
            method="knn",
            n_neighbors=n_neighbors,
            coord_type=coord_type,
        )
    else:
        raise ValueError("Must provide either 'radius' or 'n_neighbors'")


def compute_spatial_distance_matrix(
    coords: np.ndarray, max_distance: Optional[float] = None
) -> np.ndarray:
    """
    Compute pairwise Euclidean distance matrix.

    Parameters
    ----------
    coords : np.ndarray
        Nx2 array of spatial coordinates.
    max_distance : float, optional
        Maximum distance to compute. Distances beyond this are set to infinity.

    Returns
    -------
    np.ndarray
        NxN distance matrix.
    """
    from scipy.spatial.distance import pdist, squareform

    distances = squareform(pdist(coords, metric="euclidean"))

    if max_distance is not None:
        distances[distances > max_distance] = np.inf

    return distances


def find_nearest_neighbors(
    adata: anndata.AnnData,
    cell_indices: np.ndarray,
    n_neighbors: int = 10,
    spatial_key: str = "spatial",
) -> np.ndarray:
    """
    Find nearest neighbors for specific cells.

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object.
    cell_indices : np.ndarray
        Indices of cells to find neighbors for.
    n_neighbors : int
        Number of neighbors to find.
    spatial_key : str
        Key in adata.obsm containing spatial coordinates.

    Returns
    -------
    np.ndarray
        Array of shape (len(cell_indices), n_neighbors) with neighbor indices.
    """
    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")

    coords = adata.obsm[spatial_key]
    query_coords = coords[cell_indices]

    tree = cKDTree(coords)
    distances, indices = tree.query(query_coords, k=n_neighbors + 1)

    # Remove self (first neighbor)
    indices = indices[:, 1:]

    return indices
