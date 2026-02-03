"""Tests for export module."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from novae_seurat_gui.export import writers, manifest


def create_test_adata(n_obs=100, n_vars=50):
    """Create a test AnnData object with Novae results."""
    X = np.random.rand(n_obs, n_vars)
    obs = pd.DataFrame(
        {
            "cell_id": [f"cell_{i}" for i in range(n_obs)],
            "sample_id": ["sample_1"] * (n_obs // 2) + ["sample_2"] * (n_obs // 2),
            "domain": np.random.choice(["0", "1", "2", "3"], size=n_obs),
            "domain_level_0": np.random.choice(["0", "1", "2", "3"], size=n_obs),
        },
        index=[f"cell_{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])

    adata = AnnData(X=X, obs=obs, var=var)

    # Add Novae embeddings
    adata.obsm["X_novae"] = np.random.rand(n_obs, 10)
    adata.obsm["spatial"] = np.random.rand(n_obs, 2) * 1000

    return adata


class TestWriters:
    """Tests for export writers."""

    def test_export_domains(self):
        """Test exporting domains to CSV."""
        adata = create_test_adata()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "domains.csv"

            writers.export_domains(adata, str(output_file))

            assert output_file.exists()

            # Read and check
            df = pd.read_csv(output_file)
            assert "cell_id" in df.columns
            assert "domain" in df.columns or "domain_level_0" in df.columns
            assert len(df) == adata.n_obs

    def test_export_embeddings_parquet(self):
        """Test exporting embeddings to Parquet."""
        adata = create_test_adata()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "embeddings.parquet"

            writers.export_embeddings(
                adata, str(output_file), embedding_key="X_novae", format="parquet"
            )

            assert output_file.exists()

            # Read and check
            df = pd.read_parquet(output_file)
            assert "cell_id" in df.columns
            assert len(df) == adata.n_obs
            assert df.shape[1] == adata.obsm["X_novae"].shape[1] + 1  # +1 for cell_id

    def test_export_embeddings_csv(self):
        """Test exporting embeddings to CSV."""
        adata = create_test_adata()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "embeddings.csv"

            writers.export_embeddings(
                adata, str(output_file), embedding_key="X_novae", format="csv"
            )

            assert output_file.exists()

            # Read and check
            df = pd.read_csv(output_file)
            assert "cell_id" in df.columns
            assert len(df) == adata.n_obs

    def test_export_filtered_cells(self):
        """Test exporting filtered cell IDs."""
        adata = create_test_adata()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "filtered_cells.txt"

            writers.export_filtered_cells(adata, str(output_file))

            assert output_file.exists()

            # Read and check
            with open(output_file) as f:
                cell_ids = [line.strip() for line in f]

            assert len(cell_ids) == adata.n_obs

    def test_write_r_import_script(self):
        """Test writing R import script."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "import.R"

            writers.write_r_import_script(str(output_file))

            assert output_file.exists()

            # Check content
            with open(output_file) as f:
                content = f.read()

            assert "library(Seurat)" in content
            assert "read.csv" in content
            assert "read_parquet" in content

    def test_export_all(self):
        """Test exporting all files at once."""
        adata = create_test_adata()

        with tempfile.TemporaryDirectory() as tmpdir:
            exported = writers.export_all(adata, tmpdir)

            assert "domains" in exported
            assert "embeddings" in exported
            assert "filtered_cells" in exported
            assert "r_script" in exported

            # Check all files exist
            for file_path in exported.values():
                assert Path(file_path).exists()


class TestManifest:
    """Tests for manifest creation."""

    def test_create_manifest(self):
        """Test creating a run manifest."""
        adata = create_test_adata()

        manifest_data = manifest.create_manifest(
            adata,
            input_files=[],
            parameters={"n_domains": 10},
            qc_filters={"nCount_RNA": (10, None)},
            n_cells_pre_qc=150,
        )

        assert "timestamp" in manifest_data
        assert "version" in manifest_data
        assert "input" in manifest_data
        assert "processing" in manifest_data
        assert "parameters" in manifest_data

        # Check processing stats
        assert manifest_data["processing"]["n_cells_kept"] == adata.n_obs
        assert manifest_data["processing"]["n_cells_filtered"] == 50

    def test_save_manifest(self):
        """Test saving manifest to JSON."""
        adata = create_test_adata()

        manifest_data = manifest.create_manifest(adata, input_files=[])

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "manifest.json"

            manifest.save_manifest(manifest_data, str(output_file))

            assert output_file.exists()

            # Read and check
            with open(output_file) as f:
                loaded = json.load(f)

            assert loaded["version"] == manifest_data["version"]
            assert loaded["input"]["n_features"] == manifest_data["input"]["n_features"]

    def test_validate_manifest(self):
        """Test manifest validation."""
        adata = create_test_adata()

        manifest_data = manifest.create_manifest(adata, input_files=[])

        is_valid, errors = manifest.validate_manifest(manifest_data)

        assert is_valid
        assert len(errors) == 0

    def test_validate_manifest_invalid(self):
        """Test validation of invalid manifest."""
        invalid_manifest = {"timestamp": "2024-01-01"}  # Missing required keys

        is_valid, errors = manifest.validate_manifest(invalid_manifest)

        assert not is_valid
        assert len(errors) > 0
