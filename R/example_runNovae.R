# R/example_runNovae.R
# Minimal example demonstrating runNovae() on a small synthetic Seurat object.
# This script can be run interactively or in CI to validate the R helper.
#
# Requirements:
#   R packages: Seurat (>= 5), SeuratObject (>= 5), reticulate, Matrix
#   Python env: uv venv with novae_seurat_gui installed
#
# Usage (from repo root):
#   Rscript R/example_runNovae.R --python .venv/bin/python

# --- parse optional --python argument -------------------------------------
args <- commandArgs(trailingOnly = TRUE)
python_path <- NULL
idx <- which(args == "--python")
if (length(idx) > 0L && idx < length(args)) {
  python_path <- args[idx + 1L]
}

# --- load dependencies ----------------------------------------------------
for (pkg in c("Seurat", "SeuratObject", "reticulate", "Matrix")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required package '", pkg, "' is not installed. ",
      "Run: install.packages('", pkg, "')"
    )
  }
}

library(Matrix)

source(file.path(dirname(sys.frame(1L)$ofile), "runNovae.R"))

# --- build synthetic Seurat object ----------------------------------------
set.seed(42L)

n_cells <- 200L
n_genes <- 100L

cell_ids <- paste0("cell_", seq_len(n_cells))
gene_ids <- paste0("gene_", seq_len(n_genes))

# Sparse count matrix: genes x cells
counts <- Matrix::rsparsematrix(
  nrow = n_genes,
  ncol = n_cells,
  density = 0.3,
  rand.x = function(n) rpois(n, lambda = 5) + 1L
)
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids

# Create Seurat object (v5 API)
SerObj <- SeuratObject::CreateSeuratObject(counts = counts, assay = "RNA")

# --- attach synthetic spatial coordinates ---------------------------------
# Simulate a centroids-like structure using a simple list with a $coords slot
x_coords <- runif(n_cells, min = 0, max = 100)
y_coords <- runif(n_cells, min = 0, max = 100)
coord_mat <- matrix(
  c(x_coords, y_coords),
  ncol = 2L,
  dimnames = list(cell_ids, c("x", "y"))
)

# Minimal stand-in for the SeuratObject Centroids class structure.
# In real use, obj@images[[1]]$centroids is a proper SeuratObject class.
# Here we mimic it with a list so .get_spatial_coords() can extract coords.
centroids_obj <- structure(
  list(coords = coord_mat),
  class = "FakeCentroids"
)
fake_image <- structure(
  list(centroids = centroids_obj),
  class = "FakeImage"
)
SerObj@images[["slice1"]] <- fake_image

# --- verify coordinate extraction helper ----------------------------------
message("--- testing .get_spatial_coords() ...")
coords_out <- .get_spatial_coords(SerObj, img_index = 1L)
stopifnot(is.matrix(coords_out))
stopifnot(nrow(coords_out) == n_cells)
stopifnot(ncol(coords_out) == 2L)
message("    ok: ", nrow(coords_out), " x ", ncol(coords_out))

# --- verify matrix extraction helper --------------------------------------
message("--- testing .get_counts_matrix_seurat5() ...")
expr_out <- .get_counts_matrix_seurat5(SerObj, assay = "RNA", layer = "counts")
stopifnot(inherits(expr_out, "dgCMatrix"))
stopifnot(ncol(expr_out) == n_cells)
stopifnot(nrow(expr_out) == n_genes)
message("    ok: ", nrow(expr_out), " genes x ", ncol(expr_out), " cells")

# --- verify alignment helper ----------------------------------------------
message("--- testing .align_coords_to_cells() ...")
shuffled_ids <- sample(cell_ids)
aligned <- .align_coords_to_cells(coords_out, shuffled_ids)
stopifnot(identical(rownames(aligned), shuffled_ids))
message("    ok: alignment by rownames works")

# --- run full pipeline if python path provided ----------------------------
if (!is.null(python_path)) {
  message("--- running full runNovae() with python = ", python_path, " ...")
  SerObj <- runNovae(
    SerObj,
    assay            = "RNA",
    layer            = "counts",
    coords_from      = "centroids",
    model            = "MICS-Lab/novae-human-0",
    n_domains        = 5L,
    neighbors_method = "none",
    random_state     = 42L,
    python           = python_path,
    add_embedding    = TRUE,
    reduction_name   = "novae",
    verbose          = TRUE
  )

  stopifnot("novae_domain" %in% colnames(SerObj@meta.data))
  message("Domain distribution:")
  print(table(SerObj$novae_domain))

  if ("novae" %in% names(SerObj@reductions)) {
    message("Embedding dimensions: ", paste(dim(SerObj[["novae"]]), collapse = " x "))
  }
  message("--- runNovae() test passed.")
} else {
  message("Skipping full pipeline test (no --python path provided).")
  message("To run end-to-end: Rscript R/example_runNovae.R --python .venv/bin/python")
}

message("All helper tests passed.")
