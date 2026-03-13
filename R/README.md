# R/README.md

# R Helper Layer for XSpatialNovae

This folder contains `runNovae.R`, a single-file R helper that lets you run
the XSpatialNovae/Novae spatial domain-assignment pipeline directly on an
in-memory Seurat object **without writing any expression matrices to disk**.

---

## Prerequisites

| Requirement | Minimum version | Notes |
|---|---|---|
| R | >= 4.2 | |
| SeuratObject | >= 5.0.0 | Older versions used deprecated `slot=` argument |
| Seurat | >= 5.0.0 | |
| reticulate | any | CRAN: `install.packages("reticulate")` |
| Matrix | any | Usually installed with R base |

Install or update R packages if needed:

```r
install.packages(c("Seurat", "SeuratObject", "reticulate", "Matrix"))
```

---

## Python environment setup (uv-based)

XSpatialNovae uses a **uv** virtual environment. Set it up once:

```bash
# From the XSpatialNovae repository root
uv venv --python 3.11        # or 3.12 / 3.13
source .venv/bin/activate    # macOS/Linux
# or: .venv\Scripts\activate  (Windows)

uv pip install -r requirements-uv.txt
uv pip install -e .
```

That gives you a Python binary at `.venv/bin/python` (`.venv\Scripts\python.exe`
on Windows). Pass that path to `runNovae()` via the `python` argument.

---

## Quick start

```r
# 1. Source the helper (or put R/ on your source path)
source("R/runNovae.R")

# 2. Load your Seurat object
SerObj <- readRDS("my_seurat.rds")

# 3. Run Novae
SerObj <- runNovae(
  SerObj,
  python    = ".venv/bin/python",   # path to your uv venv python
  n_domains = 10,
  verbose   = TRUE
)

# 4. Inspect results
table(SerObj$novae_domain)
head(SerObj@meta.data[, "novae_domain", drop = FALSE])
```

---

## Function signature

```r
runNovae(
  SerObj,
  assay            = NULL,               # defaults to DefaultAssay(SerObj)
  layer            = "counts",           # SeuratObject v5 layer name
  coords_from      = c("centroids",      # @images[[1]]$centroids@coords
                       "coordinates"),   # @images[[1]]@coordinates
  model            = "MICS-Lab/novae-human-0",
  n_domains        = 10L,
  neighbors_radius = NULL,
  neighbors_method = c("radius", "knn", "none"),
  random_state     = 42L,
  python           = NULL,               # path to python binary
  condaenv         = NULL,               # conda env name (alternative)
  virtualenv       = NULL,               # virtualenv name/path (alternative)
  reduction_name   = "novae",
  embedding_key    = NULL,
  add_embedding    = TRUE,
  verbose          = FALSE
)
```

### Key parameters

| Parameter | Description |
|---|---|
| `assay` | Seurat assay to use (default: `DefaultAssay()`). |
| `layer` | Layer holding raw counts inside the assay. SeuratObject v5 uses `layer=` not `slot=`. |
| `coords_from` | `"centroids"` (preferred for Xenium/CosMx) or `"coordinates"` (Visium-style). |
| `model` | HuggingFace model ID. If the `novae` package is not installed, a UMAP + Leiden proxy is used automatically. |
| `n_domains` | Number of spatial domains to assign via Leiden clustering. |
| `neighbors_radius` | When `neighbors_method = "radius"`, this value is used as `n_neighbors` (rounded) for `sc.pp.neighbors`. Scanpy does not support a true radius neighbor graph natively; set `neighbors_method = "none"` and build the graph externally for radius-based graphs. |
| `neighbors_method` | `"radius"` or `"knn"` to run an explicit spatial neighbor step; `"none"` to let the pipeline choose. |
| `python` | Full path to the Python binary in your uv venv. Use one of `python`, `condaenv`, or `virtualenv`. |
| `add_embedding` | If `TRUE`, adds the Novae embedding as a `DimReduc` object named `reduction_name`. |
| `verbose` | Print progress messages. |

### Return value

The function returns the input `SerObj` with:
- `meta.data$novae_domain` - character vector of domain labels per cell.
- (optional) A `DimReduc` slot named `reduction_name` (default `"novae"`)
  containing the embedding matrix.

---

## Neighbor computation and `neighbors_radius`

Scanpy's `pp.neighbors` (used internally) does **not** support a radius-based
neighbor graph natively - it uses kNN. When `neighbors_method = "radius"` and
`neighbors_radius` is set, `neighbors_radius` is rounded to an integer and
used as `n_neighbors`.

For a true radius-based neighbor graph (e.g., with `squidpy`), set
`neighbors_method = "none"`, compute the graph yourself, store it in
`adata.obsp["connectivities"]`, and then pass the object to the Python API
directly.

In most cases, leaving `neighbors_method = "none"` (the default) works well:
`preprocess_for_novae` computes PCA and `run_novae_zeroshot` builds a kNN
graph on the PCA embedding automatically.

---

## Full end-to-end example

```r
# --- Python setup (run once in terminal) ---
# uv venv --python 3.11
# source .venv/bin/activate
# uv pip install -r requirements-uv.txt
# uv pip install -e .

# --- R ---
library(Seurat)
library(reticulate)

source("R/runNovae.R")

SerObj <- readRDS("my_xenium_seurat.rds")

# Inspect coordinate structure
head(SerObj@images[[1]]$centroids@coords)

SerObj <- runNovae(
  SerObj,
  assay            = "RNA",
  layer            = "counts",
  coords_from      = "centroids",
  model            = "MICS-Lab/novae-human-0",
  n_domains        = 10,
  neighbors_method = "none",      # let pipeline decide
  random_state     = 42,
  python           = ".venv/bin/python",
  add_embedding    = TRUE,
  reduction_name   = "novae",
  verbose          = TRUE
)

# Domain labels
print(table(SerObj$novae_domain))

# Visualize (requires Seurat >= 5)
DimPlot(SerObj, reduction = "novae", group.by = "novae_domain")
SpatialDimPlot(SerObj, group.by = "novae_domain")
```

---

## Troubleshooting

**`slot` argument defunct error**

This means your Seurat/SeuratObject is < v5 or you are calling
`GetAssayData(..., slot=)` directly. `runNovae()` uses `LayerData()` internally
and should not trigger this error. Update to SeuratObject >= 5.0.0.

**`Python module 'novae_seurat_gui' is not available`**

Install the package in your venv:
```bash
source .venv/bin/activate
uv pip install -e .
```

**`Python module 'anndata' is not available`**

Install Python dependencies:
```bash
uv pip install -r requirements-uv.txt
```

**Coordinate mismatch error**

If the number of coordinate rows does not match the number of cells, provide
an object whose image coordinates have rownames equal to the Seurat cell IDs.

**Too many or too few spatial domains**

Adjust `n_domains`. The underlying Leiden resolution is tuned towards
`n_domains`; very small or very large values may produce fewer clusters than
requested. The Leiden algorithm treats `n_domains` as a target, not a
guarantee.
