"""Tests for cluster interpretation functionality."""

import numpy as np
import pandas as pd
import pytest
import anndata as ad

from novae_seurat_gui import cluster_interpretation


@pytest.fixture
def sample_adata():
    """Create a simple test AnnData object."""
    np.random.seed(42)
    n_cells = 100
    n_genes = 50
    
    X = np.random.lognormal(0, 1, size=(n_cells, n_genes))
    spatial = np.random.rand(n_cells, 2) * 1000
    
    obs = pd.DataFrame({
        'domain': np.random.choice(['A', 'B', 'C'], n_cells),
        'cluster_leiden': np.random.choice([0, 1, 2], n_cells),
        'cell_type': np.random.choice(['T cell', 'B cell', 'Macrophage'], n_cells),
        'sample_id': np.random.choice(['sample1', 'sample2'], n_cells),
        'nCount_RNA': np.random.randint(500, 5000, n_cells),
        'nFeature_RNA': np.random.randint(100, 1000, n_cells),
    })
    
    var = pd.DataFrame({'gene_name': [f'gene_{i}' for i in range(n_genes)]})
    var.index = var['gene_name']
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm['spatial'] = spatial
    
    return adata


class TestUtils:
    """Test utility functions."""
    
    def test_get_candidate_label_columns(self, sample_adata):
        """Test getting candidate label columns."""
        candidates = cluster_interpretation.get_candidate_label_columns(sample_adata)
        
        # Should prioritize domain and cluster columns
        assert 'domain' in candidates
        assert 'cluster_leiden' in candidates
        
        # Priority columns should come first
        priority_cols = [c for c in candidates[:2] if 'domain' in c.lower() or 'cluster' in c.lower()]
        assert len(priority_cols) > 0
    
    def test_prepare_expression_data(self, sample_adata):
        """Test expression data preparation."""
        expr = cluster_interpretation.prepare_expression_data(sample_adata, normalize=True)
        
        assert expr.shape == sample_adata.X.shape
        assert isinstance(expr, np.ndarray)
        
    def test_check_spatial_coords(self, sample_adata):
        """Test spatial coordinate checking."""
        has_spatial, msg = cluster_interpretation.utils.check_spatial_coords(sample_adata)
        
        assert has_spatial is True
        assert 'spatial' in msg.lower()
    
    def test_get_sample_column(self, sample_adata):
        """Test sample column detection."""
        sample_col = cluster_interpretation.utils.get_sample_column(sample_adata)
        
        assert sample_col == 'sample_id'
    
    def test_get_celltype_column(self, sample_adata):
        """Test cell type column detection."""
        celltype_col = cluster_interpretation.utils.get_celltype_column(sample_adata)
        
        assert celltype_col == 'cell_type'
    
    def test_get_qc_columns(self, sample_adata):
        """Test QC column detection."""
        qc_cols = cluster_interpretation.utils.get_qc_columns(sample_adata)
        
        assert 'nCount_RNA' in qc_cols
        assert 'nFeature_RNA' in qc_cols


class TestSummaries:
    """Test cluster summary functions."""
    
    def test_compute_cluster_summary(self, sample_adata):
        """Test cluster summary computation."""
        summary = cluster_interpretation.compute_cluster_summary(
            sample_adata, label_col='domain', sample_col='sample_id'
        )
        
        assert 'group_id' in summary.columns
        assert 'n_cells' in summary.columns
        assert 'percent_of_total' in summary.columns
        assert len(summary) > 0
        
        # Check that percentages sum to ~100
        assert abs(summary['percent_of_total'].sum() - 100.0) < 1.0
    
    def test_compute_group_composition(self, sample_adata):
        """Test cell type composition."""
        composition = cluster_interpretation.summaries.compute_group_composition(
            sample_adata, label_col='domain', group_id='A', celltype_col='cell_type'
        )
        
        assert 'cell_type' in composition.columns
        assert 'n_cells' in composition.columns
        assert 'percent' in composition.columns
    
    def test_compare_qc_metrics(self, sample_adata):
        """Test QC metrics comparison."""
        qc_cols = ['nCount_RNA', 'nFeature_RNA']
        comparison = cluster_interpretation.summaries.compare_qc_metrics(
            sample_adata, label_col='domain', group_id='A', qc_columns=qc_cols
        )
        
        assert 'metric' in comparison.columns
        assert 'mean_in_group' in comparison.columns
        assert 'mean_other' in comparison.columns
        assert len(comparison) == 2


class TestMarkers:
    """Test marker gene computation."""
    
    def test_compute_marker_genes(self, sample_adata):
        """Test marker gene computation."""
        markers = cluster_interpretation.compute_marker_genes(
            sample_adata, label_col='domain', group_id='A', n_genes=10, normalize=True
        )
        
        assert 'gene' in markers.columns
        assert 'logFC' in markers.columns
        assert 'p_value' in markers.columns
        assert 'adj_p_value' in markers.columns
        assert 'pct_in_group' in markers.columns
        assert 'pct_out_group' in markers.columns
        
        # Should return at most n_genes
        assert len(markers) <= 10
        
        # logFC should be positive (upregulated)
        assert (markers['logFC'] > 0).all()
    
    def test_compute_fold_change(self, sample_adata):
        """Test fold change computation."""
        fc = cluster_interpretation.markers.compute_fold_change(
            sample_adata, label_col='domain', group_id='A'
        )
        
        assert 'gene' in fc.columns
        assert 'logFC' in fc.columns
        assert len(fc) == sample_adata.n_vars


class TestVisualization:
    """Test visualization functions."""
    
    def test_plot_spatial_highlight(self, sample_adata):
        """Test spatial highlight plot."""
        fig = cluster_interpretation.plot_spatial_highlight(
            sample_adata, label_col='domain', group_id='A'
        )
        
        # Check that figure was created
        assert fig is not None
        assert hasattr(fig, 'data')
        
        # Should have 2 traces (background and highlighted)
        assert len(fig.data) == 2
    
    def test_plot_celltype_composition(self, sample_adata):
        """Test cell type composition plot."""
        composition = cluster_interpretation.summaries.compute_group_composition(
            sample_adata, label_col='domain', group_id='A', celltype_col='cell_type'
        )
        
        fig = cluster_interpretation.plot_celltype_composition(composition)
        
        assert fig is not None
        assert hasattr(fig, 'data')
    
    def test_create_qc_comparison_table(self, sample_adata):
        """Test QC comparison table formatting."""
        qc_cols = ['nCount_RNA', 'nFeature_RNA']
        comparison = cluster_interpretation.summaries.compare_qc_metrics(
            sample_adata, label_col='domain', group_id='A', qc_columns=qc_cols
        )
        
        display_df = cluster_interpretation.create_qc_comparison_table(comparison)
        
        assert 'Metric' in display_df.columns
        assert 'Mean (Group)' in display_df.columns
