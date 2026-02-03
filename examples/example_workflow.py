"""
Example workflow demonstrating the Novae-Seurat-GUI pipeline.

This script shows how to:
1. Load an H5AD file from Seurat
2. Validate and detect mappings
3. Apply QC filters
4. Build spatial neighbor graph
5. Preprocess data
6. Run Novae model
7. Export results for R
"""

import logging
from pathlib import Path

from novae_seurat_gui import io, qc, spatial, modeling, export

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main():
    """Run the example workflow."""

    # ==================== 1. Load Data ====================
    logger.info("Step 1: Loading H5AD file")

    input_file = "data/seurat_export.h5ad"  # Replace with your file
    adata = io.load_h5ad(input_file)

    logger.info(f"Loaded {adata.n_obs} cells × {adata.n_vars} features")

    # ==================== 2. Validate & Detect Mappings ====================
    logger.info("Step 2: Validating schema and detecting mappings")

    # Validate
    is_valid, messages = io.validate_schema(adata, strict=False)
    if not is_valid:
        logger.error("Validation failed!")
        for msg in messages:
            logger.error(msg)
        return

    # Detect mappings
    mappings = io.detect_mappings(adata)
    logger.info(f"Detected mappings: {mappings}")

    # Apply mappings
    if mappings["x_col"] and mappings["y_col"]:
        adata = io.ensure_spatial_coords(
            adata, x_col=mappings["x_col"], y_col=mappings["y_col"]
        )

    adata = io.normalize_metadata(adata, sample_id_col=mappings["sample_id_col"])

    # ==================== 3. QC Filtering ====================
    logger.info("Step 3: Applying QC filters")

    # Define filter criteria
    filter_criteria = {
        "nCount_RNA": (10, None),
        "nFeature_RNA": (5, 1000),
    }

    # Apply filters
    mask = qc.create_filter_mask(adata, filter_criteria)
    stats = qc.compute_filter_stats(adata, mask)

    logger.info(
        f"QC filtering: {stats['n_kept']} cells kept ({stats['percent_kept']:.1f}%)"
    )

    # Create filtered AnnData
    adata_filtered = adata[mask, :].copy()

    # ==================== 4. Build Spatial Neighbor Graph ====================
    logger.info("Step 4: Building spatial neighbor graph")

    adata_filtered = spatial.compute_neighbors(
        adata_filtered, method="radius", radius=150.0
    )

    # Compute diagnostics
    diag = spatial.graph_diagnostics(adata_filtered)
    logger.info(
        f"Graph: {diag['n_edges']} edges, "
        f"avg degree {diag['degree']['mean']:.1f}"
    )

    # ==================== 5. Preprocess Data ====================
    logger.info("Step 5: Preprocessing (normalize, PCA)")

    adata_filtered = modeling.preprocess_for_novae(
        adata_filtered,
        target_sum=1e4,
        n_top_genes=2000,
        n_comps=50,
        random_state=42,
    )

    logger.info("Preprocessing complete")

    # ==================== 6. Run Novae Model ====================
    logger.info("Step 6: Running Novae model")

    adata_filtered = modeling.run_novae_zeroshot(
        adata_filtered,
        model_name="MICS-Lab/novae-human-0",
        n_domains=10,
        random_state=42,
    )

    logger.info(
        f"Novae complete: {adata_filtered.obs['domain'].nunique()} domains assigned"
    )

    # ==================== 7. Export Results ====================
    logger.info("Step 7: Exporting results for R")

    output_dir = "results"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Export all
    exported = export.export_all(
        adata_filtered,
        output_dir=output_dir,
        embedding_format="parquet",
    )

    # Create manifest
    manifest = export.create_manifest(
        adata_filtered,
        input_files=[input_file],
        parameters=adata_filtered.uns.get("novae_params", {}),
        qc_filters=filter_criteria,
        n_cells_pre_qc=adata.n_obs,
    )

    # Add spatial diagnostics
    manifest = export.add_spatial_summary_to_manifest(manifest, diag)

    # Save manifest
    manifest_file = Path(output_dir) / "run_manifest.json"
    export.save_manifest(manifest, str(manifest_file))

    logger.info("Workflow complete!")
    logger.info(f"Results exported to {output_dir}/")
    logger.info(f"  - domains.csv: Domain assignments")
    logger.info(f"  - embeddings.parquet: Novae embeddings")
    logger.info(f"  - filtered_cell_ids.txt: QC-passed cells")
    logger.info(f"  - run_manifest.json: Run metadata")
    logger.info(f"  - import_to_seurat.R: R import script")

    # ==================== Optional: Save processed H5AD ====================
    output_h5ad = Path(output_dir) / "processed_data.h5ad"
    adata_filtered.write_h5ad(output_h5ad)
    logger.info(f"Saved processed H5AD to {output_h5ad}")


if __name__ == "__main__":
    main()
