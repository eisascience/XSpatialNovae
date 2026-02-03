"""Tests for I/O module."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from novae_seurat_gui.io import loader, validator, converter


def create_test_adata(n_obs=100, n_vars=50):
    """Create a minimal test AnnData object."""
    X = np.random.rand(n_obs, n_vars)
    obs = pd.DataFrame(
        {
            "x_slide_mm": np.random.rand(n_obs) * 1000,
            "y_slide_mm": np.random.rand(n_obs) * 1000,
            "SampleID": ["sample_1"] * (n_obs // 2) + ["sample_2"] * (n_obs // 2),
            "nCount_RNA": np.random.randint(10, 1000, n_obs),
            "nFeature_RNA": np.random.randint(5, 100, n_obs),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])

    adata = AnnData(X=X, obs=obs, var=var)
    return adata


class TestLoader:
    """Tests for loader functions."""

    def test_detect_mappings(self):
        """Test automatic mapping detection."""
        adata = create_test_adata()

        mappings = loader.detect_mappings(adata)

        assert mappings["x_col"] == "x_slide_mm"
        assert mappings["y_col"] == "y_slide_mm"
        assert mappings["sample_id_col"] == "SampleID"
        assert mappings["units"] == "mm"

    def test_get_available_columns(self):
        """Test getting available columns."""
        adata = create_test_adata()

        cols = loader.get_available_columns(adata)

        assert "obs_columns" in cols
        assert "numeric_columns" in cols
        assert "categorical_columns" in cols
        assert len(cols["numeric_columns"]) >= 3  # x, y, nCount_RNA

    def test_summarize_adata(self):
        """Test AnnData summarization."""
        adata = create_test_adata(n_obs=100, n_vars=50)

        summary = loader.summarize_adata(adata)

        assert summary["n_obs"] == 100
        assert summary["n_vars"] == 50
        assert "obs_columns" in summary


class TestValidator:
    """Tests for validator functions."""

    def test_validate_schema_valid(self):
        """Test validation of valid AnnData."""
        adata = create_test_adata()
        adata.obsm["spatial"] = adata.obs[["x_slide_mm", "y_slide_mm"]].values

        is_valid, messages = validator.validate_schema(adata, strict=False)

        assert is_valid
        assert len(messages) == 0

    def test_validate_schema_no_spatial(self):
        """Test validation without spatial coordinates."""
        adata = create_test_adata()

        is_valid, messages = validator.validate_schema(adata, strict=False)

        assert is_valid  # Should pass in non-strict mode
        assert any("spatial" in msg.lower() for msg in messages)

    def test_check_required_fields(self):
        """Test checking required fields."""
        adata = create_test_adata()
        adata.obs["cell_id"] = adata.obs.index
        adata.obs["sample_id"] = adata.obs["SampleID"]

        all_present, missing = validator.check_required_fields(
            adata,
            required_obs_cols=["cell_id", "sample_id"],
        )

        assert all_present
        assert len(missing["missing_obs"]) == 0

    def test_check_counts_data(self):
        """Test checking if data looks like counts."""
        adata = create_test_adata()
        adata.X = np.random.poisson(5, size=adata.shape)  # Integer counts

        stats = validator.check_counts_data(adata)

        assert stats["is_integer"]
        assert not stats["has_negative"]
        assert stats["max_value"] >= 0


class TestConverter:
    """Tests for converter functions."""

    def test_ensure_spatial_coords(self):
        """Test creation of spatial coordinates in obsm."""
        adata = create_test_adata()

        adata = converter.ensure_spatial_coords(
            adata, x_col="x_slide_mm", y_col="y_slide_mm"
        )

        assert "spatial" in adata.obsm
        assert adata.obsm["spatial"].shape == (adata.n_obs, 2)

    def test_normalize_metadata(self):
        """Test metadata normalization."""
        adata = create_test_adata()

        adata = converter.normalize_metadata(
            adata, sample_id_col="SampleID"
        )

        assert "cell_id" in adata.obs.columns
        assert "sample_id" in adata.obs.columns

    def test_convert_units(self):
        """Test unit conversion."""
        coords = np.array([[1.0, 2.0], [3.0, 4.0]])

        # mm to px (default: 1mm = 1000px)
        coords_px = converter.convert_units(coords, from_units="mm", to_units="px")

        assert coords_px[0, 0] == 1000.0
        assert coords_px[0, 1] == 2000.0

    def test_add_cell_id_column(self):
        """Test adding cell_id column."""
        adata = create_test_adata()

        # Remove index names
        adata.obs.index = range(adata.n_obs)

        adata = converter.add_cell_id_column(adata, prefix="cell_")

        assert "cell_id" in adata.obs.columns
        assert adata.obs["cell_id"][0] == "cell_0"
