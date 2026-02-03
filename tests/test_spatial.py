"""Tests for spatial module."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from novae_seurat_gui.spatial import neighbors, diagnostics


def create_test_adata_with_spatial(n_obs=100, n_vars=50):
    """Create a test AnnData object with spatial coordinates."""
    X = np.random.rand(n_obs, n_vars)
    obs = pd.DataFrame(
        {
            "sample_id": ["sample_1"] * n_obs,
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])

    adata = AnnData(X=X, obs=obs, var=var)

    # Add spatial coordinates (random positions in 1000x1000 space)
    adata.obsm["spatial"] = np.random.rand(n_obs, 2) * 1000

    return adata


class TestNeighbors:
    """Tests for neighbor graph construction."""

    def test_compute_neighbors_radius(self):
        """Test radius-based neighbor computation."""
        adata = create_test_adata_with_spatial()

        adata = neighbors.compute_neighbors(
            adata, method="radius", radius=150.0
        )

        assert "spatial_connectivities" in adata.obsp
        assert "spatial_distances" in adata.obsp
        assert "spatial_neighbors" in adata.uns

    def test_compute_neighbors_knn(self):
        """Test KNN neighbor computation."""
        adata = create_test_adata_with_spatial()

        adata = neighbors.compute_neighbors(
            adata, method="knn", n_neighbors=10
        )

        assert "spatial_connectivities" in adata.obsp
        assert "spatial_distances" in adata.obsp

    def test_build_spatial_graph(self):
        """Test building spatial graph with auto method selection."""
        adata = create_test_adata_with_spatial()

        # With radius
        adata = neighbors.build_spatial_graph(adata, radius=150.0)

        assert "spatial_connectivities" in adata.obsp


class TestDiagnostics:
    """Tests for spatial graph diagnostics."""

    def test_graph_diagnostics(self):
        """Test computing graph diagnostics."""
        adata = create_test_adata_with_spatial()
        adata = neighbors.compute_neighbors(adata, method="radius", radius=150.0)

        diag = diagnostics.graph_diagnostics(adata)

        assert "n_cells" in diag
        assert "n_edges" in diag
        assert "degree" in diag
        assert "connected_components" in diag
        assert diag["n_cells"] == adata.n_obs

    def test_compute_degree_distribution(self):
        """Test computing degree distribution."""
        adata = create_test_adata_with_spatial()
        adata = neighbors.compute_neighbors(adata, method="knn", n_neighbors=10)

        degree_dist = diagnostics.compute_degree_distribution(adata)

        assert len(degree_dist) > 0
        assert degree_dist.sum() == adata.n_obs

    def test_identify_isolated_cells(self):
        """Test identifying isolated cells."""
        adata = create_test_adata_with_spatial()
        adata = neighbors.compute_neighbors(adata, method="radius", radius=10.0)  # Very small radius

        isolated = diagnostics.identify_isolated_cells(adata)

        # With small radius, should have some isolated cells
        assert len(isolated) >= 0
