#!/usr/bin/env Rscript
# Convert Seurat .rds with Seurat v5 compatibility (monkey-patch + RNA3 trick)

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--input_rds"), type = "character", default = NULL, help = "Input .rds file"),
  make_option(c("--output_h5ad"), type = "character", default = NULL, help = "Output .h5ad file"),
  make_option(c("--out_dir"), type = "character", default = NULL, help = "Output directory"),
  make_option(c("--assay"), type = "character", default = "RNA", help = "Assay name"),
  make_option(c("--counts_slot"), type = "character", default = "counts", help = "counts or data"),
  make_option(c("--x_col"), type = "character", default = "auto", help = "X coordinate column"),
  make_option(c("--y_col"), type = "character", default = "auto", help = "Y coordinate column"),
  make_option(c("--sample_id_col"), type = "character", default = "auto", help = "Sample ID column"),
  make_option(c("--cell_id_col"), type = "character", default = NULL, help = "Cell ID column"),
  make_option(c("--keep_meta_regex"), type = "character", default = NULL, help = "Metadata filter regex"),
  make_option(c("--overwrite"), type = "logical", default = FALSE, help = "Overwrite existing"),
  make_option(c("--verbose"), type = "logical", default = FALSE, help = "Verbose output"),
  make_option(c("--use_legacy_conversion"), type = "logical", default = TRUE, help = "Convert Assay5 to legacy Assay format for v5 compatibility")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds)) stop("--input_rds required")
if (is.null(opt$output_h5ad) && is.null(opt$out_dir)) stop("--output_h5ad or --out_dir required")

if (is.null(opt$output_h5ad)) {
  dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
  opt$output_h5ad <- file.path(opt$out_dir, paste0(tools::file_path_sans_ext(basename(opt$input_rds)), ".h5ad"))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(hdf5r)
  library(Matrix)
  library(methods)
})

is_v5 <- packageVersion("SeuratObject") >= "5.0.0"

# Monkey-patch for v5
if (is_v5) {
  setMethod(
    f = "WriteH5Group",
    signature = signature(x = "Assay"),
    definition = function(x, name, hgroup, verbose = TRUE) {
      get_mat <- function(obj, what) {
        m <- tryCatch(SeuratObject::LayerData(obj, layer = what), error = function(e) NULL)
        if (!is.null(m)) return(m)
        m <- tryCatch(slot(obj, what), error = function(e) NULL)
        if (!is.null(m)) return(m)
        Matrix::Matrix(0, 0, 0, sparse = TRUE)
      }
      
      xgroup <- hgroup$create_group(name = name)
      
      for (i in c("counts", "data", "scale.data")) {
        dat <- get_mat(x, i)
        if (inherits(dat, "matrix")) dat <- Matrix::Matrix(dat, sparse = TRUE)
        if (inherits(dat, "dgTMatrix")) dat <- as(dat, "dgCMatrix")
        if (!SeuratDisk::IsMatrixEmpty(dat)) {
          if (verbose) message("Adding ", i, " for ", name)
          SeuratDisk::WriteH5Group(dat, i, xgroup, verbose)
        }
        if (i == "scale.data" && !SeuratDisk::IsMatrixEmpty(dat)) {
          SeuratDisk::WriteH5Group(rownames(dat), "scaled.features", xgroup, verbose)
        }
      }
      
      SeuratDisk::WriteH5Group(rownames(x), "features", xgroup, verbose)
      # Note: Using ::: to access unexported GuessDType function from SeuratDisk
      # This is required by SeuratDisk's HDF5 attribute API and has been stable across versions
      xgroup$create_attr("key", Seurat::Key(x), SeuratDisk:::GuessDType(Seurat::Key(x)))
      
      if (length(Seurat::VariableFeatures(x))) {
        if (verbose) message("Adding variable features for ", name)
        SeuratDisk::WriteH5Group(Seurat::VariableFeatures(x), "variable.features", xgroup, verbose)
      }
      
      if (ncol(x[[]])) {
        if (verbose) message("Adding meta.features for ", name)
        SeuratDisk::WriteH5Group(x[[]], "meta.features", xgroup, verbose)
      }
      
      SeuratDisk::WriteH5Group(Seurat::Misc(x), "misc", xgroup, verbose)
      invisible(NULL)
    }
  )
}

# RNA3 trick
convert_rna <- function(seu, assay = "RNA") {
  if (inherits(seu@assays[[assay]], "Assay5")) {
    tmp <- paste0(assay, "3")
    seu@assays[[tmp]] <- as(seu@assays[[assay]], "Assay")
    seu@assays[[assay]] <- NULL
    names(seu@assays)[names(seu@assays) == tmp] <- assay
    DefaultAssay(seu) <- assay
  }
  seu
}

seu <- readRDS(opt$input_rds)
DefaultAssay(seu) <- opt$assay

if (opt$use_legacy_conversion && is_v5) {
  seu <- convert_rna(seu, opt$assay)
}

h5s <- sub("\\.h5ad$", ".h5seurat", opt$output_h5ad)
SaveH5Seurat(seu, h5s, overwrite = opt$overwrite)
Convert(h5s, "h5ad", overwrite = opt$overwrite, assay = opt$assay)
file.remove(h5s)

cat("✓ Conversion complete:", opt$output_h5ad, "\n")
