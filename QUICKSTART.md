# Quick Reference Guide

## Prerequisites

**Python Version**: Requires Python **3.11, 3.12, or 3.13**

```bash
python --version  # Should be 3.11.x, 3.12.x, or 3.13.x
```

## Installation

### Quick Install (uv - Recommended)

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and setup
git clone https://github.com/eisascience/XSpatialNovae.git
cd XSpatialNovae

# Create environment and install
uv venv --python 3.11
source .venv/bin/activate
uv pip install -r requirements-uv.txt
uv pip install -e .
```

### Quick Install (pip - Traditional)

```bash
git clone https://github.com/eisascience/XSpatialNovae.git
cd XSpatialNovae

# Create environment and install
python3.11 -m venv venv
source venv/bin/activate  # or `venv\Scripts\activate` on Windows
pip install -r requirements.txt
pip install -e .
```

### Download Models (Optional but Recommended)

```bash
# Prefetch Novae model weights for offline use
python scripts/download_models.py

# Or for a specific model
python scripts/download_models.py --novae-model-id MICS-Lab/novae-human-0
```

### Optional: Histology Features

**Only needed for whole-slide image analysis. Skip if you're just using spatial coordinates + expression matrices.**

```bash
# macOS
brew install openslide
pip install -r requirements-histology.txt

# Ubuntu/Debian
sudo apt-get install openslide-tools
pip install -r requirements-histology.txt

# Windows: Download from https://openslide.org/download/ and add to PATH
pip install -r requirements-histology.txt
```

### Optional: R Dependencies for .rds Upload

**Only needed if you want to upload Seurat .rds files directly. Skip if you're using H5AD files.**

#### Install R (>= 4.2)
```bash
# macOS
brew install r

# Ubuntu/Debian
sudo apt-get install r-base r-base-dev

# Windows: Download from https://cran.r-project.org/
```

#### Install Required R Packages
```bash
# Run from terminal
Rscript -e "install.packages(c('optparse', 'Seurat', 'hdf5r'))"
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes'); remotes::install_github('mojaveazure/seurat-disk')"
```

#### Verify R Installation
```bash
# Basic verification (optparse, Seurat, hdf5r)
Rscript -e "library(Seurat); library(hdf5r); library(optparse); cat('✓ Core packages ready\n')"

# Full verification (including SeuratDisk)
Rscript -e "library(Seurat); library(SeuratDisk); library(hdf5r); library(optparse); cat('✓ All R packages ready\n')"
```

**Note:** This tool supports both SeuratObject v4 (slot-based) and v5 (layer-based) automatically.

## Usage

### 1. Streamlit GUI (Recommended for Interactive Use)

```bash
streamlit run app.py
```

Navigate through the tabs:
1. **Load Data** → Upload H5AD file
2. **QC Filtering** → Filter cells by metadata
3. **Spatial Neighbors** → Build neighbor graph
4. **Run Novae** → Preprocess and run model
5. **Results** → Visualize domains and embeddings
6. **Export** → Generate R-friendly outputs

### 2. Command Line Interface

```bash
# Validate H5AD file
novae-seurat-gui validate data.h5ad

# Preprocess (normalize, PCA, neighbors)
novae-seurat-gui preprocess data.h5ad \
    --neighbors-radius 150 \
    --n-top-genes 2000 \
    --output preprocessed.h5ad

# Run Novae model
novae-seurat-gui run-novae preprocessed.h5ad \
    --model MICS-Lab/novae-human-0 \
    --n-domains 10 \
    --output novae_results.h5ad

# Generate summaries
novae-seurat-gui summarize novae_results.h5ad \
    --output summary.json

# Export for R
novae-seurat-gui export novae_results.h5ad \
    --output-dir results/
```

### 3. Python API

```python
from novae_seurat_gui import io, qc, spatial, modeling, export

# Load
adata = io.load_h5ad("data.h5ad")

# QC
mask = qc.create_filter_mask(adata, {"nCount_RNA": (10, None)})
adata = adata[mask, :].copy()

# Spatial neighbors
adata = spatial.compute_neighbors(adata, method="radius", radius=150.0)

# Preprocess & Novae
adata = modeling.preprocess_for_novae(adata, n_top_genes=2000)
adata = modeling.run_novae_zeroshot(adata, n_domains=10)

# Export
export.export_all(adata, output_dir="results/")
```

### 4. Docker

```bash
# Build and run
docker-compose up

# Access at http://localhost:8501
```

## Common Workflows

### CosMx Transcriptomics (Zero-Shot)

```bash
# Preprocess with transcriptomics defaults
novae-seurat-gui preprocess data.h5ad \
    --neighbors-radius 150 \
    --n-top-genes 2000 \
    --output preprocessed.h5ad

# Run with pretrained model
novae-seurat-gui run-novae preprocessed.h5ad \
    --model MICS-Lab/novae-human-0 \
    --n-domains 10

# Export
novae-seurat-gui export preprocessed.h5ad --output-dir results/
```

### PhenoCycler Proteomics (Train from Scratch)

In Python:
```python
from novae_seurat_gui import modeling

# Preprocess with proteomics settings (quantile scaling)
adata = modeling.preprocess_for_novae(
    adata,
    skip_normalize=True,  # Don't use standard normalization
    skip_log=True,
)
adata = modeling.apply_quantile_scaling(adata)

# Train from scratch
adata = modeling.run_novae_training(
    adata,
    n_epochs=100,
    n_domains=10,
)
```

## Configuration Files

### Default Mappings (`configs/default_mapping.yaml`)

Customize column detection:
```yaml
coordinates:
  x_column: "x_slide_mm"
  y_column: "y_slide_mm"
  units: "mm"

sample_id:
  column: "SampleID"

cell_type:
  column: "cell_type_1"
  confidence_column: "posterior_probability"
```

### Example Run (`configs/example_run.yaml`)

Full pipeline configuration:
```yaml
input:
  file: "data/seurat_export.h5ad"

qc:
  filters:
    nCount_RNA: {min: 10, max: null}
    nFeature_RNA: {min: 5, max: 1000}

spatial:
  method: "radius"
  radius: 150.0

novae:
  n_domains: 10
  model_name: "MICS-Lab/novae-human-0"
```

## Troubleshooting

### Import Errors
```bash
# Reinstall dependencies
pip install -e ".[dev]"
```

### Spatial Coordinates Not Found
```python
# Manually specify
from novae_seurat_gui.io import ensure_spatial_coords

adata = ensure_spatial_coords(
    adata,
    x_col="your_x_column",
    y_col="your_y_column"
)
```

### Too Many/Few Neighbors
```python
# Adjust radius or k
adata = spatial.compute_neighbors(
    adata,
    method="radius",
    radius=200.0  # Increase for more neighbors
)
```

## Output Files

After running the complete pipeline, expect:

```
results/
├── domains.csv              # Domain assignments per cell
├── embeddings.parquet       # Novae embeddings
├── filtered_cell_ids.txt    # Cells passing QC
├── run_manifest.json        # Complete run metadata
└── import_to_seurat.R       # R script for importing
```

## R Integration

```r
# In R
source("results/import_to_seurat.R")

# Or manually:
library(Seurat)
library(arrow)

# Load domains
domains <- read.csv("results/domains.csv", row.names = "cell_id")
seurat_obj <- AddMetaData(seurat_obj, domains)

# Load embeddings
embeddings <- arrow::read_parquet("results/embeddings.parquet")
embeddings_mat <- as.matrix(embeddings[, -1])
rownames(embeddings_mat) <- embeddings$cell_id

novae_dr <- CreateDimReducObject(
  embeddings = embeddings_mat,
  key = "NOVAe_",
  assay = DefaultAssay(seurat_obj)
)
seurat_obj[["novae"]] <- novae_dr

# Visualize
DimPlot(seurat_obj, reduction = "novae", group.by = "domain_level_0")
```

## Troubleshooting

### .rds Upload Issues

**R not found:**
```bash
# Install R
brew install r  # macOS
sudo apt-get install r-base  # Ubuntu
```

**Missing R packages:**
```bash
Rscript -e "install.packages(c('optparse', 'Seurat', 'hdf5r'))"
Rscript -e "remotes::install_github('mojaveazure/seurat-disk')"
```

**SeuratObject v5 "slot defunct" error:**
- This tool automatically supports v4 and v5
- Ensure you're using the latest version
- Use `--verbose` flag to debug

**Coordinate detection fails:**
```bash
# Specify manually
novae-seurat-gui convert data.rds --outdir ./out \
  --x-col "x_slide_mm" --y-col "y_slide_mm"
```

### Other Issues

**Import errors:**
```bash
pip install -e .
```

**Streamlit blank page:**
```bash
streamlit run app.py --server.port 8502
```

**Model not found:**
```bash
python scripts/download_models.py
```

For detailed troubleshooting, see [README.md](README.md#troubleshooting)

## Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=novae_seurat_gui

# Run specific module
pytest tests/test_io.py -v
```

## Getting Help

- **Documentation**: See [README.md](README.md) for detailed information
- **Issues**: https://github.com/eisascience/XSpatialNovae/issues
- **Contributing**: See [CONTRIBUTING.md](CONTRIBUTING.md)

## Quick Tips

1. **Start with GUI** for exploration, then use CLI for batch processing
2. **Save preprocessed H5AD** to avoid re-running expensive steps
3. **Adjust neighbor radius** based on your tissue type and coordinate units
4. **Use QC filters** to remove low-quality cells before Novae
5. **Check manifests** for complete run documentation
6. **Test with small datasets** first before running on large data
