"""Tests for QC module."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from novae_seurat_gui.qc import filters, summaries


def create_test_adata(n_obs=100, n_vars=50):
    """Create a minimal test AnnData object."""
    X = np.random.rand(n_obs, n_vars)
    obs = pd.DataFrame(
        {
            "nCount_RNA": np.random.randint(5, 1000, n_obs),
            "nFeature_RNA": np.random.randint(1, 100, n_obs),
            "Area_um2": np.random.uniform(20, 500, n_obs),
            "sample_id": ["sample_1"] * (n_obs // 2) + ["sample_2"] * (n_obs // 2),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])

    adata = AnnData(X=X, obs=obs, var=var)
    return adata


class TestFilters:
    """Tests for filter functions."""

    def test_create_filter_mask(self):
        """Test creating filter mask."""
        adata = create_test_adata()

        filter_criteria = {
            "nCount_RNA": (10, 500),
            "nFeature_RNA": (5, None),
        }

        mask = filters.create_filter_mask(adata, filter_criteria)

        assert len(mask) == adata.n_obs
        assert mask.dtype == bool
        assert np.sum(mask) <= adata.n_obs  # Some cells should be filtered

    def test_apply_qc_filters(self):
        """Test applying QC filters."""
        adata = create_test_adata(n_obs=100)

        filter_criteria = {
            "nCount_RNA": (20, 800),
        }

        adata_filtered = filters.apply_qc_filters(
            adata, filter_criteria, inplace=False
        )

        assert adata_filtered.n_obs < adata.n_obs
        assert adata_filtered.n_vars == adata.n_vars

    def test_filter_by_boolean_flags(self):
        """Test filtering by boolean flags."""
        adata = create_test_adata()
        adata.obs["qc_pass"] = np.random.choice([True, False], size=adata.n_obs)

        mask = filters.filter_by_boolean_flags(adata, ["qc_pass"], require_all=True)

        assert len(mask) == adata.n_obs
        assert np.sum(mask) == np.sum(adata.obs["qc_pass"])

    def test_filter_by_sample(self):
        """Test filtering by sample."""
        adata = create_test_adata()

        mask = filters.filter_by_sample(
            adata, sample_ids=["sample_1"], sample_col="sample_id"
        )

        assert np.sum(mask) == adata.n_obs // 2

    def test_filter_outliers_mad(self):
        """Test MAD-based outlier filtering."""
        adata = create_test_adata()

        mask = filters.filter_outliers_mad(
            adata, column="nCount_RNA", n_mads=3.0, only_upper=True
        )

        assert len(mask) == adata.n_obs
        # Most cells should pass
        assert np.sum(mask) > adata.n_obs * 0.9


class TestSummaries:
    """Tests for summary functions."""

    def test_compute_cell_statistics(self):
        """Test computing per-cell statistics."""
        adata = create_test_adata()

        stats = summaries.compute_cell_statistics(adata)

        assert len(stats) == adata.n_obs
        assert "n_counts" in stats.columns
        assert "n_features" in stats.columns
        assert all(stats["n_counts"] >= 0)

    def test_compute_qc_summary(self):
        """Test computing QC summary."""
        adata = create_test_adata()

        summary = summaries.compute_qc_summary(adata)

        assert "qc_columns" in summary
        assert "n_cells" in summary
        assert "overall" in summary
        assert summary["n_cells"] == adata.n_obs

    def test_compute_filter_stats(self):
        """Test computing filter statistics."""
        adata = create_test_adata(n_obs=100)
        mask = np.random.choice([True, False], size=100, p=[0.7, 0.3])

        stats = summaries.compute_filter_stats(adata, mask)

        assert stats["n_total"] == 100
        assert stats["n_kept"] + stats["n_filtered"] == 100
        assert 0 <= stats["percent_kept"] <= 100

    def test_identify_problematic_cells(self):
        """Test identifying problematic cells."""
        adata = create_test_adata()

        thresholds = {
            "nCount_RNA": (100, 900),
            "nFeature_RNA": (10, 90),
        }

        problematic = summaries.identify_problematic_cells(adata, thresholds)

        assert "cell_id" in problematic.columns
        assert "reasons" in problematic.columns

    def test_compare_pre_post_filtering(self):
        """Test comparing pre/post filtering statistics."""
        adata = create_test_adata(n_obs=100)

        # Create filtered version
        mask = adata.obs["nCount_RNA"] > 100
        adata_filtered = adata[mask, :].copy()

        comparison = summaries.compare_pre_post_filtering(adata, adata_filtered)

        assert "metric" in comparison.columns
        assert "pre_mean" in comparison.columns
        assert "post_mean" in comparison.columns
