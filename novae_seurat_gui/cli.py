"""Command-line interface for novae-seurat-gui."""

import logging
import sys
from pathlib import Path
from typing import Optional

import click

from . import io, qc, spatial, modeling, export


# Configure logging
def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


@click.group()
@click.version_option(version="0.1.0")
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
def main(verbose):
    """Novae-Seurat-GUI: Python workflow for Novae spatial model on Seurat objects."""
    setup_logging(verbose)


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--strict", is_flag=True, help="Enable strict validation")
def validate(input_file, strict):
    """
    Validate H5AD file schema and required fields.

    INPUT_FILE: Path to H5AD file
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Validating {input_file}")

    # Load
    adata = io.load_h5ad(input_file)

    # Validate schema
    is_valid, messages = io.validate_schema(adata, strict=strict)

    # Check mappings
    mappings = io.detect_mappings(adata)

    # Print results
    click.echo("\n=== Validation Results ===")
    click.echo(f"Status: {'PASSED' if is_valid else 'FAILED'}")
    click.echo(f"\nMessages:")
    for msg in messages:
        click.echo(f"  {msg}")

    click.echo(f"\n=== Detected Mappings ===")
    for key, value in mappings.items():
        click.echo(f"  {key}: {value}")

    sys.exit(0 if is_valid else 1)


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--output", "-o", type=click.Path(), help="Output H5AD file")
@click.option("--neighbors-radius", type=float, default=150.0, help="Neighbor radius")
@click.option("--neighbors-k", type=int, help="Number of neighbors (KNN)")
@click.option("--n-top-genes", type=int, default=2000, help="Number of HVGs")
@click.option("--n-comps", type=int, default=50, help="Number of PCA components")
@click.option("--random-state", type=int, default=42, help="Random seed")
def preprocess(
    input_file, output, neighbors_radius, neighbors_k, n_top_genes, n_comps, random_state
):
    """
    Preprocess H5AD: normalize, PCA, spatial neighbors.

    INPUT_FILE: Path to H5AD file
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Preprocessing {input_file}")

    # Load
    adata = io.load_h5ad(input_file)

    # Detect mappings
    mappings = io.detect_mappings(adata)

    # Ensure spatial coords
    if mappings["x_col"] and mappings["y_col"]:
        adata = io.ensure_spatial_coords(
            adata, x_col=mappings["x_col"], y_col=mappings["y_col"]
        )

    # Normalize metadata
    adata = io.normalize_metadata(
        adata, sample_id_col=mappings["sample_id_col"]
    )

    # Preprocess
    adata = modeling.preprocess_for_novae(
        adata,
        n_top_genes=n_top_genes,
        n_comps=n_comps,
        random_state=random_state,
    )

    # Build spatial neighbors
    if neighbors_k:
        adata = spatial.compute_neighbors(
            adata, method="knn", n_neighbors=neighbors_k
        )
    else:
        adata = spatial.compute_neighbors(
            adata, method="radius", radius=neighbors_radius
        )

    # Save
    output_file = output or input_file.replace(".h5ad", "_preprocessed.h5ad")
    logger.info(f"Saving to {output_file}")
    adata.write_h5ad(output_file)

    click.echo(f"Preprocessing complete: {output_file}")


@main.command("run-novae")
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--output", "-o", type=click.Path(), help="Output H5AD file")
@click.option("--model", default="MICS-Lab/novae-human-0", help="Pretrained model name")
@click.option("--n-domains", type=int, default=10, help="Number of domains")
@click.option("--random-state", type=int, default=42, help="Random seed")
def run_novae(input_file, output, model, n_domains, random_state):
    """
    Run Novae model (zero-shot or train).

    INPUT_FILE: Path to preprocessed H5AD file
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Running Novae on {input_file}")

    # Load
    adata = io.load_h5ad(input_file)

    # Run Novae
    adata = modeling.run_novae_zeroshot(
        adata,
        model_name=model,
        n_domains=n_domains,
        random_state=random_state,
    )

    # Save
    output_file = output or input_file.replace(".h5ad", "_novae.h5ad")
    logger.info(f"Saving to {output_file}")
    adata.write_h5ad(output_file)

    click.echo(f"Novae complete: {output_file}")
    click.echo(f"Assigned {adata.obs['domain'].nunique()} domains")


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--output", "-o", type=click.Path(), help="Output JSON file")
def summarize(input_file, output):
    """
    Generate QC and domain summaries.

    INPUT_FILE: Path to H5AD file
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Summarizing {input_file}")

    # Load
    adata = io.load_h5ad(input_file)

    # QC summary
    qc_summary = qc.compute_qc_summary(adata)

    # Spatial diagnostics (if graph exists)
    if "spatial_connectivities" in adata.obsp:
        spatial_diag = spatial.graph_diagnostics(adata)
    else:
        spatial_diag = None

    # Domain summary
    domain_summary = None
    if "domain" in adata.obs:
        domain_counts = adata.obs["domain"].value_counts().to_dict()
        domain_summary = {"n_domains": len(domain_counts), "cells_per_domain": domain_counts}

    # Combine
    summary = {
        "qc": qc_summary,
        "spatial": spatial_diag,
        "domains": domain_summary,
    }

    # Save
    import json

    output_file = output or input_file.replace(".h5ad", "_summary.json")
    with open(output_file, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    logger.info(f"Summary saved to {output_file}")
    click.echo(f"Summary: {output_file}")


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--output-dir", "-o", type=click.Path(), required=True, help="Output directory")
@click.option("--embedding-format", type=click.Choice(["parquet", "csv"]), default="parquet")
def export_cmd(input_file, output_dir, embedding_format):
    """
    Export R-friendly outputs (domains, embeddings, manifest).

    INPUT_FILE: Path to H5AD file with Novae results
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Exporting {input_file} to {output_dir}")

    # Load
    adata = io.load_h5ad(input_file)

    # Export all
    exported = export.export_all(
        adata,
        output_dir=output_dir,
        embedding_format=embedding_format,
    )

    # Create manifest
    manifest = export.create_manifest(
        adata,
        input_files=[input_file],
        parameters=adata.uns.get("novae_params", {}),
    )

    # Add diagnostics if available
    if "spatial_connectivities" in adata.obsp:
        spatial_diag = spatial.graph_diagnostics(adata)
        manifest = export.add_spatial_summary_to_manifest(manifest, spatial_diag)

    # Save manifest
    manifest_file = Path(output_dir) / "run_manifest.json"
    export.save_manifest(manifest, str(manifest_file))

    click.echo(f"\n=== Exported Files ===")
    for key, path in exported.items():
        click.echo(f"  {key}: {path}")
    click.echo(f"  manifest: {manifest_file}")


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--outdir", "-o", type=click.Path(), required=True, help="Output directory")
@click.option("--assay", default="RNA", help="Seurat assay name [default: RNA]")
@click.option("--x-col", default="auto", help="X coordinate column [default: auto-detect]")
@click.option("--y-col", default="auto", help="Y coordinate column [default: auto-detect]")
@click.option("--sample-id-col", default="auto", help="Sample ID column [default: auto-detect]")
@click.option("--cell-id-col", default=None, help="Cell ID column [default: use rownames]")
@click.option("--counts-slot", default="counts", type=click.Choice(["counts", "data"]), 
              help="Seurat slot for counts [default: counts]")
@click.option("--keep-meta-regex", default=None, help="Regex to filter metadata columns")
@click.option("--overwrite", is_flag=True, help="Overwrite existing files")
@click.option("--no-cache", is_flag=True, help="Disable caching of conversion")
def convert(input_file, outdir, assay, x_col, y_col, sample_id_col, cell_id_col, 
            counts_slot, keep_meta_regex, overwrite, no_cache):
    """
    Convert Seurat .rds file to H5AD format.

    INPUT_FILE: Path to .rds file containing Seurat object
    """
    import json
    
    logger = logging.getLogger(__name__)
    logger.info(f"Converting {input_file} to H5AD")
    
    # Check R availability
    r_available, r_msg = io.check_r_available()
    if not r_available:
        click.echo(f"ERROR: {r_msg}", err=True)
        click.echo("\nPlease install R (>=4.2) and required packages:")
        click.echo("  R -e \"install.packages(c('Seurat', 'hdf5r', 'optparse'))\"")
        click.echo("  R -e \"remotes::install_github('mojaveazure/seurat-disk')\"")
        sys.exit(1)
    
    click.echo(f"✓ {r_msg}")
    
    # Check R packages
    packages_ok, packages_msg, missing = io.check_r_packages()
    if not packages_ok:
        click.echo(f"ERROR: {packages_msg}", err=True)
        sys.exit(1)
    
    click.echo("✓ All required R packages installed")
    
    # Run conversion
    click.echo(f"\nConverting {input_file}...")
    click.echo(f"  Assay: {assay}")
    click.echo(f"  Counts slot: {counts_slot}")
    click.echo(f"  Output directory: {outdir}")
    
    try:
        output_path, adata, info = io.convert_rds_to_h5ad_with_validation(
            input_rds=input_file,
            output_dir=outdir,
            assay=assay,
            counts_slot=counts_slot,
            x_col=x_col,
            y_col=y_col,
            sample_id_col=sample_id_col,
            cell_id_col=cell_id_col,
            keep_meta_regex=keep_meta_regex,
            overwrite=overwrite,
            use_cache=not no_cache,
        )
        
        click.echo(f"\n✓ Conversion successful!")
        click.echo(f"\n=== Conversion Summary ===")
        click.echo(f"Output file: {output_path}")
        click.echo(f"Cells: {info['n_cells']}")
        click.echo(f"Features: {info['n_features']}")
        click.echo(f"Cached: {info['cached']}")
        click.echo(f"Spatial coords: {info['has_spatial']}")
        click.echo(f"Cell ID: {info['has_cell_id']}")
        click.echo(f"Sample ID: {info['has_sample_id']}")
        
        # Save manifest
        manifest_path = Path(outdir) / "conversion_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(info, f, indent=2, default=str)
        
        click.echo(f"\nManifest saved to: {manifest_path}")
        
    except Exception as e:
        click.echo(f"\nERROR: Conversion failed", err=True)
        click.echo(str(e), err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
