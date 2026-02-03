# Novae-Seurat-GUI

Python-first workflow to run the [Novae spatial foundation model](https://www.nature.com/articles/s41592-024-02465-w) on Seurat objects (R) with an interactive GUI for quality control, parameter tuning, model runs, visualization, and export back to R/Seurat.

## Overview

This toolkit bridges the gap between R/Seurat spatial analysis workflows and Python-based deep learning models for spatial biology. It provides:

- **Seamless R ↔ Python interoperability**: Load Seurat objects directly from .rds or via H5AD, process with Novae, export results back to R
- **First-class .rds support**: Upload Seurat .rds files directly in the GUI with automatic conversion to H5AD
- **Interactive Streamlit GUI**: Visual quality control, parameter tuning, and results exploration
- **Command-line interface**: Scriptable workflows for batch processing and reproducibility
- **Dual modality support**: 
  - CosMx (transcriptomics): Zero-shot inference with pretrained models
  - PhenoCycler (proteomics): Train-from-scratch workflow
- **Spatial analysis**: Build spatial neighbor graphs, assign domains, compute embeddings
- **R-friendly exports**: CSV domains, Parquet embeddings, JSON manifests, and R import snippets

## Features

### Input/Output
- **Input formats**: 
  - **.rds** (Seurat objects): Direct upload with automatic conversion
  - **.h5ad** (AnnData): Direct upload, ready for analysis
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

### Prerequisites

**Python Version**: This package requires Python **3.11, 3.12, or 3.13**. Python 3.9 and 3.10 are no longer supported, and Python 3.14+ is not yet tested.

Check your Python version:
```bash
python --version
```

### Option A: Using uv (Recommended)

[uv](https://github.com/astral-sh/uv) is a fast Python package and virtual environment manager, especially recommended for macOS Apple Silicon users.

#### 1. Install uv

```bash
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or with pip
pip install uv
```

#### 2. Create Virtual Environment

```bash
# Clone the repository
git clone https://github.com/eisascience/XSpatialNovae.git
cd XSpatialNovae

# Create a virtual environment with Python 3.11
uv venv --python 3.11

# Activate the environment
source .venv/bin/activate  # macOS/Linux
# or
.venv\Scripts\activate     # Windows
```

#### 3. Install Dependencies

```bash
# Install core dependencies
uv pip install -r requirements-uv.txt

# Install the package in development mode
uv pip install -e .

# Optional: Install development tools
uv pip install -e ".[dev]"

# Optional: Install histology extras (see "Optional Histology Features" below)
uv pip install -e ".[histology]"
```

### Option B: Using pip (Traditional)

#### 1. Create Virtual Environment

```bash
# Clone the repository
git clone https://github.com/eisascience/XSpatialNovae.git
cd XSpatialNovae

# Create virtual environment
python3.11 -m venv venv

# Activate the environment
source venv/bin/activate  # macOS/Linux
# or
venv\Scripts\activate     # Windows
```

#### 2. Install Dependencies

```bash
# Install from requirements.txt
pip install -r requirements.txt

# Install the package
pip install -e .

# Optional: Install development dependencies
pip install -e ".[dev]"

# Optional: Install histology extras (see "Optional Histology Features" below)
pip install -e ".[histology]"
```

### Optional Histology Features

The histology/whole-slide image features (OpenSlide integration) are **completely optional** and only needed if you're working with H&E or other histology images. Most users working with spatial transcriptomics coordinates and expression matrices **do not need these**.

#### When to Install Histology Dependencies

Install histology extras if you need to:
- Load and process whole-slide images (WSI)
- Work with `.svs`, `.ndpi`, or other histology formats
- Perform tissue segmentation or image-based analyses

#### Platform-Specific Installation

**macOS:**
```bash
# Install OpenSlide library via Homebrew
brew install openslide

# Install Python bindings
pip install -r requirements-histology.txt
# or
pip install -e ".[histology]"
```

**Ubuntu/Debian:**
```bash
# Install OpenSlide library
sudo apt-get update
sudo apt-get install openslide-tools

# Install Python bindings
pip install -r requirements-histology.txt
# or
pip install -e ".[histology]"
```

**Windows:**
1. Download OpenSlide binaries from https://openslide.org/download/
2. Extract and add to PATH
3. Install Python bindings:
```bash
pip install -r requirements-histology.txt
# or
pip install -e ".[histology]"
```

**Note**: If you don't need histology features, you can safely skip this entire section. The core Novae workflow for spatial transcriptomics works without OpenSlide.

### R Dependencies (for .rds Seurat Support)

If you want to upload Seurat .rds files directly in the GUI, you need to install R and several R packages. **This is optional** - you can still use the tool with H5AD files without R.

#### When to Install R Dependencies

Install R dependencies if you:
- Want to upload Seurat .rds files directly in the GUI
- Need to use the CLI `convert` command for .rds to .h5ad conversion
- Prefer working with Seurat objects in R and want seamless integration

#### R Installation

**Recommended: R version >= 4.2**

**macOS:**
```bash
# Using Homebrew
brew install r

# Or download from CRAN
# https://cran.r-project.org/bin/macosx/
```

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install r-base r-base-dev
```

**Windows:**
1. Download R from https://cran.r-project.org/bin/windows/base/
2. Run the installer
3. Ensure R is added to your system PATH

#### Required R Packages

After installing R, install these packages:

```r
# Start R and run:
install.packages(c("optparse", "Seurat", "hdf5r"))

# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Install SeuratDisk from GitHub
remotes::install_github("mojaveazure/seurat-disk")
```

**Note on SeuratObject Versions:**

This tool supports both **SeuratObject v4** (slot-based) and **SeuratObject v5** (layer-based) seamlessly. The conversion script automatically detects your SeuratObject version and uses the appropriate API:
- **v5 (≥5.0.0)**: Uses `LayerData()` and layer-based access
- **v4 (<5.0.0)**: Uses `GetAssayData()` with slot-based access

If you're using SeuratObject v5 and encounter errors about "slot is defunct", ensure you're using the latest version of this tool which includes v5 compatibility fixes.

#### Verify Installation

Check that R and packages are correctly installed:

```bash
# Check R is available
Rscript --version

# Check packages (from terminal)
Rscript -e "library(Seurat); library(SeuratDisk); library(hdf5r); library(optparse); cat('All packages loaded successfully\n')"
```

Or use the built-in checker:
```bash
# From Python
python -c "from novae_seurat_gui.io import check_r_available, check_r_packages; print(check_r_available()); print(check_r_packages())"
```

**Note**: Without R dependencies, you can still use the GUI with H5AD files. The GUI will show a friendly error message if you try to upload .rds files without R installed.

### Download Models / Assets

Novae model weights are hosted on [Hugging Face](https://huggingface.co/MICS-Lab) and may be downloaded automatically when first used. However, for **offline/cluster environments** or to avoid delays during first run, you can prefetch and cache them.

#### Automatic Download (Online Use)

Models will be downloaded automatically on first use and cached in:
- `~/.cache/huggingface/hub` (default)
- Or the location specified by `HF_HOME`, `TRANSFORMERS_CACHE`, or `HUGGINGFACE_HUB_CACHE` environment variables

#### Manual Prefetch (Recommended for Offline/Cluster)

```bash
# Download default Novae model (novae-human-0)
python scripts/download_models.py

# Download a specific model version
python scripts/download_models.py --novae-model-id MICS-Lab/novae-human-1

# Use a custom cache directory
python scripts/download_models.py --cache-dir /scratch/models

# For offline use, set the cache location before running the app
export HF_HOME=/path/to/cache
streamlit run app.py
```

The script will print the cache location and provide instructions for configuring your environment to use the cached models.

**Available Models:**
- `MICS-Lab/novae-human-0` - Default pretrained model for human tissue
- `MICS-Lab/novae-human-1` - Alternative/updated version (check Hugging Face for availability)

**Note**: If you encounter network issues or "model not found" errors, ensure you have:
1. An active internet connection
2. Access to Hugging Face Hub (no firewall blocking)
3. Run `python scripts/download_models.py` to prefetch models

### Verify Installation

```bash
# Check package is installed
python -c "import novae_seurat_gui; print('✓ Package installed successfully')"

# Run tests (optional)
pytest

# Start the GUI
streamlit run app.py
```

## Quick Start

### GUI Mode

```bash
streamlit run app.py
```

Then navigate to `http://localhost:8501` in your browser.

### CLI Mode

```bash
# Convert Seurat .rds to H5AD
novae-seurat-gui convert mydata.rds --outdir ./converted --assay RNA

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

### 1. Prepare Data

**Option A: Direct .rds Upload (Recommended)**

No preparation needed! Simply save your Seurat object in R and upload the .rds file:

```r
# In R: Save Seurat object
saveRDS(seurat_obj, "mydata.rds")
```

Then upload `mydata.rds` directly in the GUI or use the CLI:

```bash
# Convert via CLI
novae-seurat-gui convert mydata.rds --outdir ./converted \
  --assay RNA \
  --x-col x_slide_mm \
  --y-col y_slide_mm \
  --sample-id-col SampleID
```

**Option B: Manual H5AD Export (Alternative)**

If you prefer, you can manually export to H5AD in R:

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

```bash
streamlit run app.py
```

Then in the browser:
- **Upload .rds or .h5ad file**: Drag and drop or browse
  - For .rds: Configure assay and coordinate mappings
  - For .h5ad: Ready to use immediately
- Review and adjust detected mappings (coordinates, sample ID, cell types)
- Apply QC filters with live preview
- Build spatial neighbors and review diagnostics
- Run Novae model (zero-shot or train-from-scratch)
- Visualize results (domains, embeddings, compositions)
- Export R-friendly outputs

**Option B: Command Line**
```bash
# Convert .rds to .h5ad (if not done already)
novae-seurat-gui convert mydata.rds --outdir ./data

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

### Input Formats

The tool accepts two input formats:

#### 1. Seurat .rds Files (Recommended for R Users)

Upload Seurat objects directly as .rds files. The tool will automatically:
- Extract the specified assay (default: "RNA")
- Extract spatial coordinates from metadata
- Extract cell metadata
- Convert to H5AD format with proper schema
- Cache the conversion for future use

**Seurat object should contain:**
- Expression data in an assay (e.g., "RNA")
- Spatial coordinates in `@meta.data` (e.g., `x_slide_mm`, `y_slide_mm`)
- Sample identifiers in `@meta.data` (e.g., `SampleID`, `Tissue`)
- Optional: Cell type annotations and other metadata

#### 2. AnnData .h5ad Files

Pre-converted H5AD files with the following structure:

- **adata.obs**: 
  - `cell_id` (unique identifier) - created automatically if missing
  - `sample_id` (sample/tissue identifier) - auto-detected or created
  - Metadata columns (nCount_RNA, nFeature_RNA, etc.)
- **adata.obsm["spatial"]**: Nx2 array with x, y coordinates - created from obs columns if missing
- **adata.layers["counts"]** or **adata.X**: Raw counts matrix

### Auto-Detection and Mapping

The GUI auto-detects common column names but allows manual override:
- **Coordinates**: `x_slide_mm`/`y_slide_mm` (preferred), `x_FOV_px`/`y_FOV_px`, `x`/`y`
- **Sample ID**: `SampleID` (preferred), `sample_id`, `Tissue`, `tissue`
- **Cell types**: `cell_type_1`, `cell_type`, `celltype`
- **Confidence**: `posterior_probability`, `confidence`

For .rds files, you can specify these mappings during upload, or use "auto" to let the tool detect them.

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
│   │   ├── converter.py
│   │   └── convert.py             # RDS to H5AD conversion
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
├── scripts/                       # Utility scripts
│   ├── download_models.py
│   └── convert_seurat_rds_to_h5ad.R  # R conversion script
├── configs/                       # Configuration templates
│   ├── default_mapping.yaml
│   └── example_run.yaml
├── tests/                         # Unit tests
│   ├── test_io.py
│   ├── test_convert.py            # Conversion tests
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

## Troubleshooting

### .rds Conversion Issues

#### "Rscript not found in PATH"
**Problem**: R is not installed or not in your system PATH.

**Solution**:
1. Install R (>= 4.2) from https://cran.r-project.org/
2. Ensure R is in your PATH:
   ```bash
   # Test if R is accessible
   Rscript --version
   ```
3. If still not working on macOS: `brew install r`
4. If still not working on Ubuntu: `sudo apt-get install r-base r-base-dev`

#### "optparse package is required but not installed"
**Problem**: The optparse R package is missing.

**Solution**:
```bash
Rscript -e "install.packages('optparse')"
```

#### "GetAssayData slot is defunct" or "slot not found"
**Problem**: You're using SeuratObject v5, which changed from slot-based to layer-based access.

**Solution**: This tool automatically detects and handles both v4 and v5. If you see this error:
1. Ensure you're using the latest version of this tool
2. Run with `--verbose` flag to see version detection:
   ```bash
   novae-seurat-gui convert mydata.rds --outdir ./output --verbose
   ```
3. If the error persists, your Seurat object may have unusual layer structure. Try:
   ```r
   # In R: Check your Seurat object structure
   library(Seurat)
   obj <- readRDS("mydata.rds")
   print(packageVersion("SeuratObject"))
   print(Assays(obj))
   print(Layers(obj[["RNA"]]))  # For v5
   ```

#### "Missing R packages: Seurat, SeuratDisk, hdf5r"
**Problem**: Required R packages are not installed.

**Solution**:
```bash
# Install core packages
Rscript -e "install.packages(c('Seurat', 'hdf5r'))"

# Install remotes if needed
Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')"

# Install SeuratDisk from GitHub
Rscript -e "remotes::install_github('mojaveazure/seurat-disk')"
```

#### "Conversion timed out after 10 minutes"
**Problem**: Very large Seurat object or slow system.

**Solution**:
1. Check object size: `object.size(seurat_obj)` in R
2. Consider subsetting the data first in R before conversion
3. Close other applications to free up memory
4. The conversion is cached - subsequent runs with same parameters will be instant

#### Conversion succeeds but spatial coordinates missing
**Problem**: Auto-detection couldn't find coordinate columns.

**Solution**:
```bash
# Manually specify coordinate columns
novae-seurat-gui convert mydata.rds --outdir ./output \
  --x-col "your_x_column_name" \
  --y-col "your_y_column_name"
```

Or in the GUI: Adjust the detected mappings before clicking "Load Dataset"

#### "Error: Assay RNA not found"
**Problem**: Your Seurat object doesn't have the default "RNA" assay.

**Solution**:
```r
# In R: Check available assays
Assays(seurat_obj)
```

Then specify the correct assay:
```bash
novae-seurat-gui convert mydata.rds --outdir ./output --assay "Spatial"
```

### Other Common Issues

#### "No module named 'novae_seurat_gui'"
**Problem**: Package not installed.

**Solution**:
```bash
pip install -e .
```

#### GUI shows blank page
**Problem**: Streamlit not started correctly or port conflict.

**Solution**:
```bash
# Try a different port
streamlit run app.py --server.port 8502

# Or check if another process is using 8501
lsof -ti:8501 | xargs kill -9  # Kill process on port 8501 (macOS/Linux)
```

#### "Model not found" errors
**Problem**: Novae model weights not downloaded.

**Solution**:
```bash
# Prefetch the model
python scripts/download_models.py

# Or set cache directory
export HF_HOME=/path/to/cache
python scripts/download_models.py --cache-dir /path/to/cache
```

#### Import errors with scanpy/anndata
**Problem**: Version conflicts or missing dependencies.

**Solution**:
```bash
# Reinstall with clean environment
pip uninstall scanpy anndata -y
pip install -e .
```

For more help, please open an issue at: https://github.com/eisascience/XSpatialNovae/issues

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
