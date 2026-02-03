"""Tests for RDS to H5AD conversion module."""

import tempfile
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from novae_seurat_gui.io import convert


class TestConversionUtils:
    """Tests for conversion utility functions."""

    def test_check_r_available(self):
        """Test R availability check."""
        is_available, message = convert.check_r_available()
        
        # This test will pass or fail depending on environment
        # We just check the function returns expected types
        assert isinstance(is_available, bool)
        assert isinstance(message, str)
        assert len(message) > 0
    
    def test_generate_cache_key(self, tmp_path):
        """Test cache key generation."""
        # Create a dummy file
        test_file = tmp_path / "test.rds"
        test_file.write_text("dummy content")
        
        params1 = {"assay": "RNA", "x_col": "x"}
        params2 = {"assay": "RNA", "x_col": "y"}
        
        key1 = convert.generate_cache_key(str(test_file), params1)
        key2 = convert.generate_cache_key(str(test_file), params1)
        key3 = convert.generate_cache_key(str(test_file), params2)
        
        # Same file + params should give same key
        assert key1 == key2
        
        # Different params should give different key
        assert key1 != key3
        
        # Keys should be hex strings of reasonable length
        assert len(key1) == 16
        assert all(c in "0123456789abcdef" for c in key1)


class TestEnsureAnnDataSchema:
    """Tests for AnnData schema enforcement after conversion."""
    
    def create_test_adata(self, has_spatial=False, has_cell_id=False, 
                         has_sample_id=False, n_obs=100, n_vars=50):
        """Create a minimal test AnnData object."""
        X = np.random.rand(n_obs, n_vars)
        obs = pd.DataFrame(
            {
                "x_slide_mm": np.random.rand(n_obs) * 1000,
                "y_slide_mm": np.random.rand(n_obs) * 1000,
                "SampleID": ["sample_1"] * (n_obs // 2) + ["sample_2"] * (n_obs // 2),
            },
            index=[f"cell_{i}" for i in range(n_obs)],
        )
        
        if has_cell_id:
            obs["cell_id"] = obs.index
        
        if has_sample_id:
            obs["sample_id"] = obs["SampleID"]
        
        var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])
        
        adata = AnnData(X=X, obs=obs, var=var)
        
        if has_spatial:
            adata.obsm["spatial"] = obs[["x_slide_mm", "y_slide_mm"]].values
        
        return adata
    
    def test_ensure_spatial_coords_from_obs(self, tmp_path):
        """Test creation of obsm['spatial'] from obs columns."""
        adata = self.create_test_adata(has_spatial=False)
        
        # Save to file
        h5ad_path = tmp_path / "test.h5ad"
        adata.write_h5ad(h5ad_path)
        
        # Ensure schema
        adata_processed = convert.ensure_anndata_schema(
            str(h5ad_path),
            x_col="x_slide_mm",
            y_col="y_slide_mm",
            overwrite=True,
        )
        
        assert "spatial" in adata_processed.obsm
        assert adata_processed.obsm["spatial"].shape == (adata.n_obs, 2)
        
        # Check values are correct
        expected_coords = adata.obs[["x_slide_mm", "y_slide_mm"]].values
        np.testing.assert_array_equal(adata_processed.obsm["spatial"], expected_coords)
    
    def test_ensure_spatial_coords_auto_detect(self, tmp_path):
        """Test auto-detection of coordinate columns."""
        adata = self.create_test_adata(has_spatial=False)
        
        h5ad_path = tmp_path / "test.h5ad"
        adata.write_h5ad(h5ad_path)
        
        # Don't specify columns, should auto-detect
        adata_processed = convert.ensure_anndata_schema(
            str(h5ad_path),
            overwrite=True,
        )
        
        assert "spatial" in adata_processed.obsm
    
    def test_ensure_spatial_coords_already_exists(self, tmp_path):
        """Test that existing spatial coords are preserved."""
        adata = self.create_test_adata(has_spatial=True)
        
        h5ad_path = tmp_path / "test.h5ad"
        adata.write_h5ad(h5ad_path)
        
        original_spatial = adata.obsm["spatial"].copy()
        
        adata_processed = convert.ensure_anndata_schema(
            str(h5ad_path),
            overwrite=True,
        )
        
        # Should preserve existing spatial coords
        np.testing.assert_array_equal(adata_processed.obsm["spatial"], original_spatial)
    
    def test_ensure_cell_id(self, tmp_path):
        """Test creation of cell_id from index."""
        adata = self.create_test_adata(has_cell_id=False)
        
        h5ad_path = tmp_path / "test.h5ad"
        adata.write_h5ad(h5ad_path)
        
        adata_processed = convert.ensure_anndata_schema(
            str(h5ad_path),
            x_col="x_slide_mm",
            y_col="y_slide_mm",
            overwrite=True,
        )
        
        assert "cell_id" in adata_processed.obs.columns
        # Should be string version of index
        assert list(adata_processed.obs["cell_id"]) == list(adata.obs.index.astype(str))
    
    def test_ensure_sample_id_from_column(self, tmp_path):
        """Test creation of sample_id from specified column."""
        adata = self.create_test_adata(has_sample_id=False)
        
        h5ad_path = tmp_path / "test.h5ad"
        adata.write_h5ad(h5ad_path)
        
        adata_processed = convert.ensure_anndata_schema(
            str(h5ad_path),
            x_col="x_slide_mm",
            y_col="y_slide_mm",
            sample_id_col="SampleID",
            overwrite=True,
        )
        
        assert "sample_id" in adata_processed.obs.columns
        # Should match the SampleID column
        assert list(adata_processed.obs["sample_id"]) == list(adata.obs["SampleID"].astype(str))
    
    def test_ensure_sample_id_auto_detect(self, tmp_path):
        """Test auto-detection of sample_id column."""
        adata = self.create_test_adata(has_sample_id=False)
        
        h5ad_path = tmp_path / "test.h5ad"
        adata.write_h5ad(h5ad_path)
        
        adata_processed = convert.ensure_anndata_schema(
            str(h5ad_path),
            x_col="x_slide_mm",
            y_col="y_slide_mm",  # Fixed typo
            overwrite=True,
        )
        
        assert "sample_id" in adata_processed.obs.columns
    
    def test_ensure_sample_id_default(self, tmp_path):
        """Test default sample_id creation when no column found."""
        adata = self.create_test_adata(has_sample_id=False)
        
        # Remove SampleID column
        adata.obs = adata.obs.drop(columns=["SampleID"])
        
        h5ad_path = tmp_path / "test.h5ad"
        adata.write_h5ad(h5ad_path)
        
        adata_processed = convert.ensure_anndata_schema(
            str(h5ad_path),
            x_col="x_slide_mm",
            y_col="y_slide_mm",
            overwrite=True,
        )
        
        assert "sample_id" in adata_processed.obs.columns
        # Should all be "sample_0"
        assert (adata_processed.obs["sample_id"] == "sample_0").all()
    
    def test_missing_coords_raises_error(self, tmp_path):
        """Test that missing coords without detection raises error."""
        adata = self.create_test_adata(has_spatial=False)
        
        # Remove coord columns
        adata.obs = adata.obs.drop(columns=["x_slide_mm", "y_slide_mm"])
        
        h5ad_path = tmp_path / "test.h5ad"
        adata.write_h5ad(h5ad_path)
        
        with pytest.raises(ValueError, match="Cannot create obsm"):
            convert.ensure_anndata_schema(
                str(h5ad_path),
                x_col="nonexistent_x",
                y_col="nonexistent_y",
                overwrite=True,
            )
    
    def test_full_workflow(self, tmp_path):
        """Test full workflow: missing everything gets fixed."""
        adata = self.create_test_adata(
            has_spatial=False,
            has_cell_id=False,
            has_sample_id=False,
        )
        
        h5ad_path = tmp_path / "test.h5ad"
        adata.write_h5ad(h5ad_path)
        
        adata_processed = convert.ensure_anndata_schema(
            str(h5ad_path),
            x_col="x_slide_mm",
            y_col="y_slide_mm",
            sample_id_col="SampleID",
            overwrite=True,
        )
        
        # All required fields should now exist
        assert "spatial" in adata_processed.obsm
        assert "cell_id" in adata_processed.obs.columns
        assert "sample_id" in adata_processed.obs.columns
        
        # Verify file was saved with modifications
        adata_reloaded = anndata.read_h5ad(h5ad_path)
        assert "spatial" in adata_reloaded.obsm
        assert "cell_id" in adata_reloaded.obs.columns
        assert "sample_id" in adata_reloaded.obs.columns


class TestConvertRdsToH5ad:
    """Tests for RDS to H5AD conversion (integration tests)."""
    
    def test_r_not_available_raises_error(self, tmp_path, monkeypatch):
        """Test that missing R raises appropriate error."""
        # Mock R check to return False
        def mock_check_r():
            return False, "R not found"
        
        monkeypatch.setattr(convert, "check_r_available", mock_check_r)
        
        dummy_rds = tmp_path / "test.rds"
        dummy_rds.write_text("dummy")
        
        with pytest.raises(RuntimeError, match="R is not available"):
            convert.convert_rds_to_h5ad(
                input_rds=str(dummy_rds),
                output_dir=str(tmp_path),
            )
    
    def test_missing_packages_raises_error(self, tmp_path, monkeypatch):
        """Test that missing R packages raises appropriate error."""
        # Mock R available but packages missing
        def mock_check_r():
            return True, "R version 4.2"
        
        def mock_check_packages():
            return False, "Missing: Seurat, SeuratDisk", ["Seurat", "SeuratDisk"]
        
        monkeypatch.setattr(convert, "check_r_available", mock_check_r)
        monkeypatch.setattr(convert, "check_r_packages", mock_check_packages)
        
        dummy_rds = tmp_path / "test.rds"
        dummy_rds.write_text("dummy")
        
        with pytest.raises(RuntimeError, match="Required R packages missing"):
            convert.convert_rds_to_h5ad(
                input_rds=str(dummy_rds),
                output_dir=str(tmp_path),
            )
