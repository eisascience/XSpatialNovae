#!/usr/bin/env Rscript
# Convert Seurat .rds object to H5AD format via H5Seurat
# 
# Required R packages:
#   - Seurat
#   - SeuratDisk
#   - hdf5r
#   - optparse

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Error: optparse package is required but not installed.\n",
         "Install with: R -e \"install.packages('optparse')\"")
  }
  library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("--input_rds"), type = "character", default = NULL,
              help = "Path to input .rds file (Seurat object)", metavar = "FILE"),
  make_option(c("--output_h5ad"), type = "character", default = NULL,
              help = "Path to output .h5ad file", metavar = "FILE"),
  make_option(c("--assay"), type = "character", default = "RNA",
              help = "Assay name to extract [default: %default]", metavar = "NAME"),
  make_option(c("--counts_slot"), type = "character", default = "counts",
              help = "Slot to use for counts data: 'counts' or 'data' [default: %default]", metavar = "SLOT"),
  make_option(c("--x_col"), type = "character", default = "auto",
              help = "Column name for x coordinates [default: auto-detect]", metavar = "COL"),
  make_option(c("--y_col"), type = "character", default = "auto",
              help = "Column name for y coordinates [default: auto-detect]", metavar = "COL"),
  make_option(c("--sample_id_col"), type = "character", default = "auto",
              help = "Column name for sample ID [default: auto-detect]", metavar = "COL"),
  make_option(c("--cell_id_col"), type = "character", default = NULL,
              help = "Column name for cell ID [default: use rownames]", metavar = "COL"),
  make_option(c("--keep_meta_regex"), type = "character", default = NULL,
              help = "Regex pattern to filter metadata columns to keep", metavar = "REGEX"),
  make_option(c("--overwrite"), type = "logical", default = FALSE,
              help = "Overwrite existing output files [default: %default]"),
  make_option(c("--verbose"), type = "logical", default = FALSE,
              help = "Print verbose output including session info [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list,
                          description = "\nConvert Seurat .rds object to H5AD format")
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input_rds)) {
  stop("Error: --input_rds is required")
}
if (is.null(opt$output_h5ad)) {
  stop("Error: --output_h5ad is required")
}

# Check if input file exists
if (!file.exists(opt$input_rds)) {
  stop(paste("Error: Input file does not exist:", opt$input_rds))
}

# Check if output file exists and overwrite is FALSE
if (file.exists(opt$output_h5ad) && !opt$overwrite) {
  stop(paste("Error: Output file already exists:", opt$output_h5ad, 
             "\nUse --overwrite=TRUE to overwrite"))
}

cat("=== Seurat RDS to H5AD Conversion ===\n")
cat("Input RDS:", opt$input_rds, "\n")
cat("Output H5AD:", opt$output_h5ad, "\n")
cat("Assay:", opt$assay, "\n")
cat("Counts slot:", opt$counts_slot, "\n")

if (opt$verbose) {
  cat("\n=== Session Info ===\n")
  print(sessionInfo())
  cat("\n")
}

# Load required libraries
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Error: Seurat package is not installed.\nInstall with: install.packages('Seurat')")
  }
  if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    stop("Error: SeuratDisk package is not installed.\nInstall with: remotes::install_github('mojaveazure/seurat-disk')")
  }
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Error: hdf5r package is not installed.\nInstall with: install.packages('hdf5r')")
  }
  
  library(Seurat)
  library(SeuratDisk)
  library(hdf5r)
})

# Detect SeuratObject version for compatibility
seurat_obj_version <- try(packageVersion("SeuratObject"), silent = TRUE)
is_seurat_v5 <- FALSE
if (!inherits(seurat_obj_version, "try-error")) {
  is_seurat_v5 <- seurat_obj_version >= "5.0.0"
  if (opt$verbose) {
    cat("SeuratObject version:", as.character(seurat_obj_version), 
        ifelse(is_seurat_v5, "(v5 layer-based)", "(v4 slot-based)"), "\n")
  }
}

# Helper function to safely extract data from Seurat assay (v4/v5 compatible)
get_assay_data <- function(object, assay_name, layer_or_slot = "counts") {
  assay_obj <- object[[assay_name]]
  
  if (is_seurat_v5) {
    # SeuratObject v5: Use layer-based access
    if (inherits(assay_obj, "Assay5")) {
      # Get available layers
      available_layers <- Layers(assay_obj)
      
      if (opt$verbose) {
        cat("  Assay5 detected. Available layers:", paste(available_layers, collapse = ", "), "\n")
      }
      
      # Try to get the requested layer
      if (layer_or_slot == "data") {
        layer_name <- paste0("data.", DefaultLayer(assay_obj))
      } else {
        layer_name <- paste0("counts.", DefaultLayer(assay_obj))
      }
      
      # Check various layer naming patterns
      if (layer_or_slot %in% available_layers) {
        data <- LayerData(assay_obj, layer = layer_or_slot)
      } else if (layer_name %in% available_layers) {
        data <- LayerData(assay_obj, layer = layer_name)
      } else if (length(available_layers) > 0) {
        # Fall back to first available layer with appropriate prefix
        matching_layers <- grep(paste0("^", layer_or_slot), available_layers, value = TRUE)
        if (length(matching_layers) > 0) {
          data <- LayerData(assay_obj, layer = matching_layers[1])
          warning("Using layer '", matching_layers[1], "' as '", layer_or_slot, "' not found")
        } else {
          # Use default layer
          default_layer <- DefaultLayer(assay_obj)
          data <- LayerData(assay_obj, layer = paste0(layer_or_slot, ".", default_layer))
        }
      } else {
        stop("No layers available in Assay5 object")
      }
    } else {
      # Legacy Assay in v5
      data <- GetAssayData(assay_obj, layer = layer_or_slot)
    }
  } else {
    # SeuratObject v4 or earlier: Use slot-based access
    data <- GetAssayData(assay_obj, slot = layer_or_slot)
  }
  
  return(data)
}

# Helper function to check if data exists
has_assay_data <- function(object, assay_name, layer_or_slot = "counts") {
  tryCatch({
    data <- get_assay_data(object, assay_name, layer_or_slot)
    return(!is.null(data) && length(data) > 0)
  }, error = function(e) {
    return(FALSE)
  })
}

# Load Seurat object
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(opt$input_rds)

if (!inherits(seurat_obj, "Seurat")) {
  stop("Error: Input file does not contain a Seurat object")
}

cat("Loaded Seurat object with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "features\n")

# Check if assay exists
available_assays <- Assays(seurat_obj)
cat("Available assays:", paste(available_assays, collapse = ", "), "\n")

if (!(opt$assay %in% available_assays)) {
  stop(paste("Error: Assay", opt$assay, "not found in Seurat object.\n",
             "Available assays:", paste(available_assays, collapse = ", ")))
}

# Set default assay
DefaultAssay(seurat_obj) <- opt$assay
cat("Using assay:", opt$assay, "\n")

# Check counts slot using version-compatible helper
has_counts <- has_assay_data(seurat_obj, opt$assay, "counts")

if (opt$counts_slot == "counts" && !has_counts) {
  warning("Counts data is empty or not available. Falling back to 'data' slot/layer.")
  opt$counts_slot <- "data"
}

cat("Using counts slot/layer:", opt$counts_slot, "\n")

# Filter metadata columns if regex is provided
if (!is.null(opt$keep_meta_regex)) {
  cat("Filtering metadata columns with regex:", opt$keep_meta_regex, "\n")
  meta_cols <- colnames(seurat_obj@meta.data)
  keep_cols <- grep(opt$keep_meta_regex, meta_cols, value = TRUE)
  cat("Keeping", length(keep_cols), "out of", length(meta_cols), "metadata columns\n")
  seurat_obj@meta.data <- seurat_obj@meta.data[, keep_cols, drop = FALSE]
}

# Auto-detect or validate coordinate columns
meta_cols <- colnames(seurat_obj@meta.data)

detect_column <- function(candidates, available_cols, col_type) {
  for (candidate in candidates) {
    if (candidate %in% available_cols) {
      cat("Auto-detected", col_type, "column:", candidate, "\n")
      return(candidate)
    }
  }
  return(NULL)
}

# Detect x column
if (opt$x_col == "auto") {
  x_candidates <- c("x_slide_mm", "x_location", "x_centroid", "x_FOV_px", "x", "X")
  opt$x_col <- detect_column(x_candidates, meta_cols, "x coordinate")
  if (is.null(opt$x_col)) {
    warning("Could not auto-detect x coordinate column. Will check obsm['spatial'] after conversion.")
  }
} else if (!(opt$x_col %in% meta_cols)) {
  stop(paste("Error: x column", opt$x_col, "not found in metadata"))
}

# Detect y column
if (opt$y_col == "auto") {
  y_candidates <- c("y_slide_mm", "y_location", "y_centroid", "y_FOV_px", "y", "Y")
  opt$y_col <- detect_column(y_candidates, meta_cols, "y coordinate")
  if (is.null(opt$y_col)) {
    warning("Could not auto-detect y coordinate column. Will check obsm['spatial'] after conversion.")
  }
} else if (!(opt$y_col %in% meta_cols)) {
  stop(paste("Error: y column", opt$y_col, "not found in metadata"))
}

# Detect sample_id column
if (opt$sample_id_col == "auto") {
  sample_candidates <- c("SampleID", "sample_id", "Tissue", "tissue", "sample")
  opt$sample_id_col <- detect_column(sample_candidates, meta_cols, "sample_id")
  if (is.null(opt$sample_id_col)) {
    warning("Could not auto-detect sample_id column. Will create default sample_id after conversion.")
  }
} else if (!(opt$sample_id_col %in% meta_cols)) {
  stop(paste("Error: sample_id column", opt$sample_id_col, "not found in metadata"))
}

# Convert to H5Seurat
h5seurat_path <- sub("\\.h5ad$", ".h5Seurat", opt$output_h5ad)
cat("\nConverting to H5Seurat format...\n")
cat("Intermediate H5Seurat file:", h5seurat_path, "\n")

tryCatch({
  SaveH5Seurat(seurat_obj, filename = h5seurat_path, overwrite = opt$overwrite)
  cat("Successfully saved H5Seurat file\n")
}, error = function(e) {
  stop(paste("Error saving H5Seurat file:", e$message))
})

# Convert H5Seurat to H5AD
cat("\nConverting H5Seurat to H5AD...\n")
tryCatch({
  Convert(h5seurat_path, dest = "h5ad", overwrite = opt$overwrite, assay = opt$assay)
  cat("Successfully converted to H5AD\n")
}, error = function(e) {
  # Clean up H5Seurat file
  if (file.exists(h5seurat_path)) {
    file.remove(h5seurat_path)
  }
  stop(paste("Error converting to H5AD:", e$message))
})

# Clean up intermediate H5Seurat file
if (file.exists(h5seurat_path)) {
  file.remove(h5seurat_path)
  cat("Removed intermediate H5Seurat file\n")
}

# Output conversion summary
cat("\n=== Conversion Summary ===\n")
cat("Output file:", opt$output_h5ad, "\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of features:", nrow(seurat_obj), "\n")
cat("Assay used:", opt$assay, "\n")

if (!is.null(opt$x_col) && !is.null(opt$y_col)) {
  cat("Coordinate columns:", opt$x_col, ",", opt$y_col, "\n")
} else {
  cat("Coordinate columns: Not detected in metadata (check obsm['spatial'])\n")
}

if (!is.null(opt$sample_id_col)) {
  cat("Sample ID column:", opt$sample_id_col, "\n")
  n_samples <- length(unique(seurat_obj@meta.data[[opt$sample_id_col]]))
  cat("Number of samples:", n_samples, "\n")
} else {
  cat("Sample ID column: Not detected (will create default)\n")
}

cat("\n✓ Conversion complete!\n")
cat("Note: Python post-processing will ensure correct AnnData schema.\n")
