# runNovae.R
# R helper layer for running Novae on an in-memory Seurat object via reticulate.
# Requires: reticulate, Matrix, Seurat (>= 5.0), SeuratObject (>= 5.0)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Extract the counts matrix from a Seurat v5 or v4 assay.
#'
#' Tries SeuratObject::LayerData() first (v5), falls back to
#' Seurat::GetAssayData() with the layer= argument (also v5-compatible).
#'
#' @param SerObj  A Seurat object.
#' @param assay   Assay name (character).
#' @param layer   Layer name, e.g. "counts" (character).
#' @return A sparse matrix (genes x cells), coerced to dgCMatrix.
#' @keywords internal
.get_counts_matrix_seurat5 <- function(SerObj, assay, layer) {
  expr <- tryCatch(
    SeuratObject::LayerData(SerObj[[assay]], layer = layer),
    error = function(e) NULL
  )

  if (is.null(expr)) {
    expr <- tryCatch(
      Seurat::GetAssayData(SerObj, assay = assay, layer = layer),
      error = function(e) NULL
    )
  }

  if (is.null(expr)) {
    stop(
      "Could not retrieve expression data from assay '", assay,
      "', layer '", layer, "'. ",
      "Ensure the Seurat object is v5-compatible and the layer name is correct."
    )
  }

  if (!inherits(expr, "dgCMatrix")) {
    expr <- methods::as(expr, "dgCMatrix")
  }

  expr
}


#' Extract spatial coordinates from a Seurat image slot.
#'
#' Tries (in order):
#'   1. obj@images[[img_index]]$centroids@coords  (preferred, Xenium/CosMx)
#'   2. obj@images[[img_index]]@coordinates       (Visium-style)
#'
#' @param SerObj    A Seurat object.
#' @param img_index Integer index or character name of the image to use.
#' @return A numeric matrix with columns x and y, rownames = cell IDs.
#' @keywords internal
.get_spatial_coords <- function(SerObj, img_index = 1L) {
  if (length(SerObj@images) == 0L) {
    stop("No images found in Seurat object. Spatial coordinates are required.")
  }

  img <- SerObj@images[[img_index]]

  coords <- tryCatch(
    {
      raw <- img$centroids@coords
      if (!is.matrix(raw)) raw <- as.matrix(raw)
      raw
    },
    error = function(e) NULL
  )

  if (is.null(coords)) {
    coords <- tryCatch(
      {
        raw <- img@coordinates
        if (!is.matrix(raw)) raw <- as.matrix(raw)
        raw
      },
      error = function(e) NULL
    )
  }

  if (is.null(coords)) {
    stop(
      "Could not extract spatial coordinates from Seurat images slot. ",
      "Tried obj@images[[", img_index, "]]$centroids@coords ",
      "and obj@images[[", img_index, "]]@coordinates. ",
      "Ensure the Seurat object contains valid spatial coordinate data."
    )
  }

  if (ncol(coords) < 2L) {
    stop("Spatial coordinates matrix must have at least 2 columns (x and y).")
  }

  coords[, 1:2, drop = FALSE]
}


#' Align spatial coordinates to expression cell IDs.
#'
#' If coords has rownames and all cell IDs are present in those rownames,
#' reorders rows to match cell_ids. Otherwise returns coords as-is (assuming
#' rows already correspond to the columns of the expression matrix).
#'
#' @param coords   Numeric matrix, rows = cells.
#' @param cell_ids Character vector of cell IDs from the expression matrix.
#' @return Aligned coordinate matrix.
#' @keywords internal
.align_coords_to_cells <- function(coords, cell_ids) {
  rn <- rownames(coords)
  if (!is.null(rn) && length(rn) > 0L && all(cell_ids %in% rn)) {
    return(coords[cell_ids, , drop = FALSE])
  }
  if (nrow(coords) != length(cell_ids)) {
    stop(
      "Number of spatial coordinate rows (", nrow(coords), ") does not match ",
      "number of cells in expression matrix (", length(cell_ids), "). ",
      "Provide coordinate rownames matching cell IDs for unambiguous alignment."
    )
  }
  coords
}


# ---------------------------------------------------------------------------
# Main user-facing function
# ---------------------------------------------------------------------------

#' Run Novae spatial domain assignment on a Seurat object
#'
#' Extracts expression and spatial coordinates from a Seurat object (v5
#' SeuratObject API), builds an in-memory AnnData via reticulate, runs the
#' XSpatialNovae/novae_seurat_gui Python pipeline, and returns domain labels
#' (and optionally embeddings) back into the Seurat object.
#'
#' No intermediate files are written to disk.
#'
#' @param SerObj         A Seurat object (SeuratObject >= 5.0).
#' @param assay          Assay to use. Defaults to \code{DefaultAssay(SerObj)}.
#' @param layer          Layer within the assay containing raw counts.
#'                       Defaults to \code{"counts"}.
#' @param coords_from    Where to look for spatial coordinates.
#'                       \code{"centroids"} uses
#'                       \code{obj@images[[1]]$centroids@coords};
#'                       \code{"coordinates"} uses
#'                       \code{obj@images[[1]]@coordinates}.
#' @param model          Pretrained Novae model identifier.
#'                       Default: \code{"MICS-Lab/novae-human-0"}.
#' @param n_domains      Number of spatial domains to assign. Default: 10.
#' @param neighbors_radius  Spatial neighbor radius in coordinate units.
#'                          Passed to \code{sc.pp.neighbors()} when
#'                          \code{neighbors_method} is \code{"radius"}.
#'                          Set to \code{NULL} to use the scanpy default.
#' @param neighbors_method  Neighbor computation strategy: \code{"radius"},
#'                          \code{"knn"}, or \code{"none"} to skip explicit
#'                          neighbor computation (let the pipeline decide).
#' @param random_state   Random seed. Default: 42.
#' @param python         Path to the Python binary in your uv virtual
#'                       environment, e.g. \code{".venv/bin/python"}.
#'                       Passed to \code{reticulate::use_python()}.
#' @param condaenv       Name of a conda environment to activate via
#'                       \code{reticulate::use_condaenv()}.
#'                       Use either \code{python} or \code{condaenv}, not both.
#' @param virtualenv     Name or path of a virtualenv to activate via
#'                       \code{reticulate::use_virtualenv()}.
#'                       Use either \code{python} or \code{virtualenv}, not both.
#' @param reduction_name Name for the Seurat DimReduc slot that receives the
#'                       Novae embedding. Default: \code{"novae"}.
#' @param embedding_key  Key to look up in \code{adata.obsm} for the embedding.
#'                       When \code{NULL} the function tries in order:
#'                       \code{X_novae}, \code{X_emb}, \code{X_umap},
#'                       \code{X_pca}.
#' @param add_embedding  If \code{TRUE} (default) and an embedding is found,
#'                       add it as a DimReduc to \code{SerObj}.
#' @param verbose        If \code{TRUE}, print progress messages.
#'
#' @return The input \code{SerObj} with:
#'   \itemize{
#'     \item \code{meta.data$novae_domain} containing the domain labels.
#'     \item (optional) A DimReduc named \code{reduction_name} containing the
#'           Novae embedding.
#'   }
#'
#' @examples
#' \dontrun{
#' SerObj <- runNovae(
#'   SerObj,
#'   python = ".venv/bin/python",
#'   n_domains = 10
#' )
#' table(SerObj$novae_domain)
#' }
#'
#' @export
runNovae <- function(
  SerObj,
  assay = NULL,
  layer = "counts",
  coords_from = c("centroids", "coordinates"),
  model = "MICS-Lab/novae-human-0",
  n_domains = 10L,
  neighbors_radius = NULL,
  neighbors_method = c("radius", "knn", "none"),
  random_state = 42L,
  python = NULL,
  condaenv = NULL,
  virtualenv = NULL,
  reduction_name = "novae",
  embedding_key = NULL,
  add_embedding = TRUE,
  verbose = FALSE
) {
  # --- prerequisite packages -------------------------------------------
  for (pkg in c("reticulate", "Matrix")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", pkg, "' is required. Install it with: ",
        "install.packages('", pkg, "')"
      )
    }
  }
  for (pkg in c("Seurat", "SeuratObject")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", pkg, "' is required. Install it with: ",
        "install.packages('", pkg, "')"
      )
    }
  }

  coords_from    <- match.arg(coords_from)
  neighbors_method <- match.arg(neighbors_method)

  # --- assay ---------------------------------------------------------------
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(SerObj)
  }

  if (verbose) message("Using assay: ", assay, ", layer: ", layer)

  # --- expression matrix ---------------------------------------------------
  if (verbose) message("Extracting expression matrix ...")
  expr <- .get_counts_matrix_seurat5(SerObj, assay, layer)

  cell_ids <- colnames(expr)
  gene_ids <- rownames(expr)

  if (is.null(cell_ids) || length(cell_ids) == 0L) {
    stop("Expression matrix has no column names (cell IDs).")
  }
  if (is.null(gene_ids) || length(gene_ids) == 0L) {
    stop("Expression matrix has no row names (gene/feature IDs).")
  }

  if (verbose) {
    message(
      "Expression matrix: ", nrow(expr), " genes x ", ncol(expr), " cells"
    )
  }

  # --- spatial coordinates -------------------------------------------------
  if (verbose) message("Extracting spatial coordinates ...")
  coords <- .get_spatial_coords(SerObj, img_index = 1L)
  coords <- .align_coords_to_cells(coords, cell_ids)

  if (verbose) {
    message(
      "Coordinate range: x [", round(min(coords[, 1]), 3),
      ", ", round(max(coords[, 1]), 3), "]",
      "  y [", round(min(coords[, 2]), 3),
      ", ", round(max(coords[, 2]), 3), "]"
    )
  }

  # --- configure Python environment ----------------------------------------
  n_env_args <- sum(c(!is.null(python), !is.null(condaenv), !is.null(virtualenv)))
  if (n_env_args > 1L) {
    stop(
      "Provide at most one of 'python', 'condaenv', or 'virtualenv', not multiple."
    )
  }

  if (!is.null(python)) {
    reticulate::use_python(python, required = TRUE)
  } else if (!is.null(condaenv)) {
    reticulate::use_condaenv(condaenv, required = TRUE)
  } else if (!is.null(virtualenv)) {
    reticulate::use_virtualenv(virtualenv, required = TRUE)
  }

  # --- import Python packages ----------------------------------------------
  if (verbose) message("Importing Python packages ...")

  for (mod in c("anndata", "numpy", "scipy.sparse", "pandas")) {
    if (!reticulate::py_module_available(mod)) {
      stop(
        "Python module '", mod, "' is not available. ",
        "Install it in your environment, e.g.: uv pip install ", mod
      )
    }
  }
  if (!reticulate::py_module_available("novae_seurat_gui")) {
    stop(
      "Python module 'novae_seurat_gui' is not available. ",
      "Install XSpatialNovae in your environment: uv pip install -e ."
    )
  }

  anndata_mod <- reticulate::import("anndata")
  np          <- reticulate::import("numpy")
  sp          <- reticulate::import("scipy.sparse")
  pd          <- reticulate::import("pandas")
  nsg         <- reticulate::import("novae_seurat_gui")

  # --- build in-memory AnnData ---------------------------------------------
  if (verbose) message("Building in-memory AnnData ...")

  # dgCMatrix -> scipy.sparse.csc_matrix
  # We pass the three CSC components explicitly so that the conversion is
  # unambiguous regardless of the reticulate version.
  expr_data    <- expr@x
  expr_indices <- expr@i          # 0-based row indices (CSC)
  expr_indptr  <- expr@p          # column pointers
  n_genes      <- nrow(expr)
  n_cells      <- ncol(expr)

  expr_csc <- sp$csc_matrix(
    reticulate::tuple(
      np$array(expr_data,    dtype = np$float64),
      np$array(expr_indices, dtype = np$int32),
      np$array(expr_indptr,  dtype = np$int32)
    ),
    shape = reticulate::tuple(n_genes, n_cells)
  )
  # Transpose to cells x genes as AnnData expects X to be obs x var
  expr_csc_t <- expr_csc$T$tocsc()

  obs_df <- pd$DataFrame(
    list(cell_id = cell_ids),
    index = cell_ids
  )
  var_df <- pd$DataFrame(
    list(gene_id = gene_ids),
    index = gene_ids
  )

  adata <- anndata_mod$AnnData(
    X   = expr_csc_t,
    obs = obs_df,
    var = var_df
  )

  adata$obsm[["spatial"]] <- np$array(coords, dtype = np$float64)

  if (verbose) {
    message(
      "AnnData built: ", reticulate::py_to_r(adata$n_obs), " obs x ",
      reticulate::py_to_r(adata$n_vars), " vars"
    )
  }

  # --- spatial neighbors (optional) ----------------------------------------
  # Note: scanpy's pp.neighbors does not support a radius parameter natively;
  # it uses kNN. When neighbors_method = "radius" and neighbors_radius is set,
  # neighbors_radius is used as n_neighbors (rounded to an integer) as a rough
  # approximation. For true radius-based graphs, set neighbors_method = "none"
  # and build the graph externally before calling runNovae(), or use squidpy.
  if (neighbors_method != "none") {
    sc_mod <- reticulate::import("scanpy")
    n_neighbors_arg <- if (
      neighbors_method == "radius" && !is.null(neighbors_radius)
    ) {
      as.integer(round(neighbors_radius))
    } else {
      15L  # scanpy default
    }
    if (verbose) {
      message(
        "Computing spatial kNN neighbors (n_neighbors = ", n_neighbors_arg,
        ", use_rep = spatial) ..."
      )
    }
    sc_mod$pp$neighbors(
      adata,
      n_neighbors  = n_neighbors_arg,
      n_pcs        = 0L,
      use_rep      = "spatial",
      random_state = as.integer(random_state)
    )
  }

  # --- preprocess ----------------------------------------------------------
  if (verbose) message("Preprocessing for Novae ...")
  adata <- nsg$modeling$preprocess_for_novae(
    adata,
    random_state = as.integer(random_state)
  )

  # --- run Novae -----------------------------------------------------------
  if (verbose) message("Running Novae zero-shot model ...")
  adata <- nsg$modeling$run_novae_zeroshot(
    adata,
    model_name   = model,
    n_domains    = as.integer(n_domains),
    random_state = as.integer(random_state)
  )

  # --- pull results back into Seurat ---------------------------------------
  if (verbose) message("Retrieving domain labels ...")

  domain_py <- adata$obs[["domain"]]
  domain_r  <- reticulate::py_to_r(domain_py)

  if (is.null(names(domain_r))) {
    names(domain_r) <- cell_ids
  } else {
    # Reorder to match cell_ids in case Python returned a different order
    domain_r <- domain_r[cell_ids]
  }

  SerObj <- Seurat::AddMetaData(
    SerObj,
    metadata = domain_r,
    col.name = "novae_domain"
  )

  if (verbose) {
    tab <- table(domain_r)
    message(
      "Domain assignment complete. ", length(tab), " domain(s) found."
    )
  }

  # --- optional embedding --------------------------------------------------
  if (add_embedding) {
    obsm_keys <- tryCatch(
      reticulate::py_to_r(adata$obsm$keys()),
      error = function(e) character(0L)
    )

    candidate_keys <- if (!is.null(embedding_key)) {
      embedding_key
    } else {
      c("X_novae", "X_emb", "X_umap", "X_pca")
    }

    found_key <- NULL
    for (k in candidate_keys) {
      if (k %in% obsm_keys) {
        found_key <- k
        break
      }
    }

    if (!is.null(found_key)) {
      if (verbose) message("Adding embedding '", found_key, "' as DimReduc '", reduction_name, "' ...")
      emb <- reticulate::py_to_r(adata$obsm[[found_key]])
      if (!is.matrix(emb)) emb <- as.matrix(emb)
      rownames(emb) <- cell_ids
      colnames(emb) <- paste0(toupper(reduction_name), "_", seq_len(ncol(emb)))

      SerObj[[reduction_name]] <- Seurat::CreateDimReducObject(
        embeddings = emb,
        key        = paste0(toupper(reduction_name), "_"),
        assay      = assay
      )
    } else if (verbose) {
      message(
        "No embedding found in adata.obsm (tried: ",
        paste(candidate_keys, collapse = ", "), "). ",
        "DimReduc will not be added."
      )
    }
  }

  if (verbose) message("runNovae complete.")

  SerObj
}
