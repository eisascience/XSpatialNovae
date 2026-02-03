# Novae-Seurat-GUI

Python-first workflow to run the [Novae spatial foundation model](https://www.nature.com/articles/s41592-024-02465-w) on Seurat objects (R) exported to AnnData (Python), with an interactive GUI for quality control, parameter tuning, model runs, visualization, and export back to R/Seurat.

## Overview

This toolkit bridges the gap between R/Seurat spatial analysis workflows and Python-based deep learning models for spatial biology. It provides:

- **Seamless R ↔ Python interoperability**: Load Seurat objects exported as H5AD, process with Novae, export results back to R
- **Interactive Streamlit GUI**: Visual quality control, parameter tuning, and results exploration
- **Command-line interface**: Scriptable workflows for batch processing and reproducibility
- **Dual modality support**: 
  - CosMx (transcriptomics): Zero-shot inference with pretrained models
  - PhenoCycler (proteomics): Train-from-scratch workflow
- **Spatial analysis**: Build spatial neighbor graphs, assign domains, compute embeddings
- **R-friendly exports**: CSV domains, Parquet embeddings, JSON manifests, and R import snippets

## Features

### Input/Output
- Import Seurat objects via H5Seurat/H5AD interchange format
- Auto-detect spatial coordinates, sample identifiers, and cell metadata
- Export domain assignments, embeddings, and QC-filtered cell lists
- Generate R code snippets for seamless Seurat integration

### Quality Control
- Interactive filtering based on metadata (nCount_RNA, nFeature_RNA, Area, QC flags, etc.)
- Real-time spatial visualizations showing kept vs. filtered cells
- Summary statistics and plots

### Spatial Analysis
- Build spatial neighbor graphs (radius-based or k-nearest neighbors)
- Graph diagnostics: degree distribution, connected components, edge visualization
- Downsample for efficient visualization of large datasets

### Novae Model
- **CosMx mode**: Load pretrained models (e.g., MICS-Lab/novae-human-0) for zero-shot inference
- **PhenoCycler mode**: Train from scratch with appropriate preprocessing
- Configurable parameters: neighbors, domains, embedding dimensions, epochs
- Intelligent caching for responsive GUI experience

### Visualization
- Spatial scatter plots colored by domains, niches, QC flags, cell types
- Embedding plots (PCA, UMAP) colored by metadata
- Per-domain summaries: cell counts, composition, marker features

### Niche Analysis (Optional)
- Compute neighborhood composition features using cell type labels
- Cluster neighborhoods into niche states
- Confidence-weighted composition using posterior probabilities
- Export niche assignments and summaries

## Installation

### From source

```bash
git clone https://github.com/eisascience/XSpatialNovae.git
cd XSpatialNovae
pip install -e .
```

### For development

```bash
pip install -e ".[dev]"
```

## Quick Start

### GUI Mode

```bash
streamlit run app.py
```

Then navigate to `http://localhost:8501` in your browser.

### CLI Mode

```bash
# Validate H5AD schema
novae-seurat-gui validate data.h5ad

# Preprocess: normalize, PCA, neighbors
novae-seurat-gui preprocess data.h5ad --output processed.h5ad

# Run Novae (zero-shot)
novae-seurat-gui run-novae processed.h5ad --model MICS-Lab/novae-human-0

# Generate summaries
novae-seurat-gui summarize processed.h5ad --output summary.json

# Export R-friendly outputs
novae-seurat-gui export processed.h5ad --output-dir results/
```

## Usage

### 1. Prepare Data from R/Seurat

```r
# In R: Export Seurat object to H5AD
library(Seurat)
library(SeuratDisk)

# Convert Seurat to H5AD
SaveH5Seurat(seurat_obj, filename = "data.h5Seurat")
Convert("data.h5Seurat", dest = "h5ad")
```

### 2. Run in Python (GUI or CLI)

**Option A: Streamlit GUI**
- Load H5AD file
- Review and adjust detected mappings (coordinates, sample ID, cell types)
- Apply QC filters with live preview
- Build spatial neighbors and review diagnostics
- Run Novae model (zero-shot or train-from-scratch)
- Visualize results (domains, embeddings, compositions)
- Export R-friendly outputs

**Option B: Command Line**
```bash
# Full pipeline
novae-seurat-gui validate data.h5ad
novae-seurat-gui preprocess data.h5ad --neighbors-radius 150
novae-seurat-gui run-novae data.h5ad --model MICS-Lab/novae-human-0 --n-domains 10
novae-seurat-gui export data.h5ad --output-dir results/
```

### 3. Import Results back to R/Seurat

```r
# In R: Merge domains and embeddings
source("results/import_to_seurat.R")

# Or manually:
domains <- read.csv("results/domains.csv", row.names = "cell_id")
seurat_obj <- AddMetaData(seurat_obj, domains)

embeddings <- arrow::read_parquet("results/embeddings.parquet")
novae_dr <- CreateDimReducObject(
  embeddings = as.matrix(embeddings[, -1]), 
  key = "NOVAe_", 
  assay = "RNA"
)
seurat_obj[["novae"]] <- novae_dr
```

## Data Requirements

The input H5AD must contain:

- **adata.obs**: 
  - `cell_id` (unique identifier)
  - `sample_id` (sample/tissue identifier)
  - Metadata columns (nCount_RNA, nFeature_RNA, etc.)
- **adata.obsm["spatial"]**: Nx2 array with x, y coordinates
- **adata.layers["counts"]** or **adata.X**: Raw counts matrix

The GUI auto-detects common column names but allows manual override:
- Coordinates: `x_slide_mm`/`y_slide_mm` (preferred), `x_FOV_px`/`y_FOV_px`
- Sample ID: `SampleID` (preferred), `Tissue`, or custom
- Cell types: `cell_type_1`, confidence: `posterior_probability`

## Configuration

Default mappings and parameters can be specified in YAML config files:

```yaml
# configs/default_mapping.yaml
coordinates:
  x_column: "x_slide_mm"
  y_column: "y_slide_mm"
  units: "mm"

sample_id:
  column: "SampleID"

cell_type:
  column: "cell_type_1"
  confidence_column: "posterior_probability"

qc_filters:
  nCount_RNA_min: 10
  nFeature_RNA_min: 5
  Area_um2_min: 50
```

## Project Structure

```
XSpatialNovae/
├── app.py                          # Streamlit GUI entry point
├── novae_seurat_gui/              # Main package
│   ├── __init__.py
│   ├── cli.py                     # Command-line interface
│   ├── io/                        # Data loading and validation
│   │   ├── loader.py
│   │   ├── validator.py
│   │   └── converter.py
│   ├── qc/                        # Quality control
│   │   ├── filters.py
│   │   └── summaries.py
│   ├── spatial/                   # Spatial analysis
│   │   ├── neighbors.py
│   │   └── diagnostics.py
│   ├── modeling/                  # Novae wrappers
│   │   ├── novae_runner.py
│   │   └── parameters.py
│   ├── niche/                     # Niche analysis (optional)
│   │   ├── composition.py
│   │   └── clustering.py
│   ├── viz/                       # Visualization utilities
│   │   ├── spatial_plots.py
│   │   └── embedding_plots.py
│   └── export/                    # R-friendly exports
│       ├── writers.py
│       └── manifest.py
├── configs/                       # Configuration templates
│   ├── default_mapping.yaml
│   └── example_run.yaml
├── tests/                         # Unit tests
│   ├── test_io.py
│   ├── test_qc.py
│   ├── test_spatial.py
│   └── test_export.py
├── examples/                      # Example scripts
│   ├── example_workflow.py
│   └── example_notebook.ipynb
├── pyproject.toml                # Package configuration
├── requirements.txt              # Dependencies
└── README.md                     # This file
```

## Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=novae_seurat_gui --cov-report=html

# Run specific test module
pytest tests/test_io.py
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes with tests
4. Run tests and linting (`pytest`, `black .`, `flake8`)
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## Citation

If you use this tool in your research, please cite:

- **Novae paper**: [Citation for Novae model paper]
- **This repository**: [Citation for this repo when published]

## License

MIT License - see LICENSE file for details

## Acknowledgments

- Novae spatial foundation model developers
- Seurat and Scanpy communities
- Contributors and users

## Support

- **Issues**: https://github.com/eisascience/XSpatialNovae/issues
- **Discussions**: https://github.com/eisascience/XSpatialNovae/discussions
- **Email**: info@eisascience.com
