## 可复用的 10x/Seurat RNA QC 与 DoubletFinder 函数。

read_sigcell_rna <- function(input_file, sample_id) {
  old_object <- readRDS(input_file)
  if (!inherits(old_object, "Seurat")) {
    stop("Input is not a Seurat object: ", input_file)
  }
  if (!"RNA" %in% SeuratObject::Assays(old_object) ||
      !"counts" %in% SeuratObject::Layers(old_object[["RNA"]])) {
    stop("Input requires RNA counts: ", input_file)
  }

  object <- Seurat::CreateSeuratObject(
    counts = SeuratObject::LayerData(old_object, assay = "RNA", layer = "counts"),
    project = sample_id,
    min.cells = 3,
    min.features = 0
  )
  object$sample_id <- sample_id

  mitochondrial_genes <- grep("^(mt-|MT-)", rownames(object), value = TRUE)
  if (length(mitochondrial_genes) == 0L) {
    warning("No mitochondrial genes matched in ", sample_id)
    object$percent.mt <- 0
  } else {
    object[["percent.mt"]] <- Seurat::PercentageFeatureSet(
      object,
      features = mitochondrial_genes
    )
  }
  object
}

filter_sigcell_qc <- function(object, min_features, max_features, max_percent_mt) {
  keep_cells <-
    object$nFeature_RNA >= min_features &
    object$nFeature_RNA <= max_features &
    object$percent.mt < max_percent_mt
  subset(object, cells = colnames(object)[keep_cells])
}

sigcell_df_function <- function(candidates) {
  package_namespace <- asNamespace("DoubletFinder")
  available <- candidates[vapply(
    candidates,
    exists,
    logical(1),
    envir = package_namespace,
    inherits = FALSE
  )]
  if (length(available) == 0L) {
    stop("DoubletFinder does not provide: ", paste(candidates, collapse = ", "))
  }
  get(available[[1L]], envir = package_namespace, inherits = FALSE)
}

estimate_sigcell_doublet_rate <- function(n_cells, rate_per_1000, max_rate) {
  min(max_rate, rate_per_1000 * n_cells / 1000)
}

run_sigcell_doubletfinder <- function(
  object,
  sample_id,
  normalization_method,
  scale_factor,
  variable_feature_method,
  n_variable_features,
  n_pcs,
  dims_to_use,
  precluster_resolution,
  doublet_rate_per_1000,
  doublet_rate_max,
  doublet_pN
) {
  if (ncol(object) < 100L) {
    stop(sample_id, " has fewer than 100 QC-passed cells; DoubletFinder is not suitable.")
  }

  # 这些预处理仅用于 DoubletFinder 判定 doublet，不作为最终对象保存。
  object <- Seurat::NormalizeData(
    object,
    normalization.method = normalization_method,
    scale.factor = scale_factor,
    verbose = FALSE
  )
  object <- Seurat::FindVariableFeatures(
    object,
    selection.method = variable_feature_method,
    nfeatures = n_variable_features,
    verbose = FALSE
  )
  object <- Seurat::ScaleData(
    object,
    features = Seurat::VariableFeatures(object),
    verbose = FALSE
  )
  object <- Seurat::RunPCA(
    object,
    features = Seurat::VariableFeatures(object),
    npcs = n_pcs,
    verbose = FALSE
  )

  available_pcs <- ncol(Seurat::Embeddings(object, "pca"))
  if (is.null(dims_to_use) || !length(dims_to_use)) stop("dims_to_use must be explicit.")
  df_dims <- seq_len(min(max(dims_to_use), available_pcs))
  selected_pc <- max(df_dims)
  if (length(df_dims) < 2L) {
    stop(sample_id, " does not have enough PCA dimensions for DoubletFinder.")
  }
  object <- Seurat::FindNeighbors(
    object,
    reduction = "pca",
    dims = df_dims,
    verbose = FALSE
  )
  object <- Seurat::FindClusters(
    object,
    resolution = precluster_resolution,
    cluster.name = "doublet_preclusters",
    verbose = FALSE
  )

  param_sweep <- sigcell_df_function(c("paramSweep", "paramSweep_v3"))
  summarize_sweep <- sigcell_df_function(c("summarizeSweep", "summarizeSweep_v3"))
  find_pk <- sigcell_df_function("find.pK")
  model_homotypic <- sigcell_df_function("modelHomotypic")
  call_doublets <- sigcell_df_function(c("doubletFinder", "doubletFinder_v3"))

  sweep_results <- param_sweep(object, PCs = df_dims, sct = FALSE)
  sweep_statistics <- summarize_sweep(sweep_results, GT = FALSE)
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  pk_statistics <- find_pk(sweep_statistics)
  valid_pk <- which(
    is.finite(pk_statistics$BCmetric) &
      !is.na(suppressWarnings(as.numeric(as.character(pk_statistics$pK))))
  )
  if (length(valid_pk) == 0L) {
    stop("No usable pK was found for ", sample_id)
  }
  best_index <- valid_pk[which.max(pk_statistics$BCmetric[valid_pk])]
  best_pk <- as.numeric(as.character(pk_statistics$pK[best_index]))

  expected_rate <- estimate_sigcell_doublet_rate(
    ncol(object),
    rate_per_1000 = doublet_rate_per_1000,
    max_rate = doublet_rate_max
  )
  expected_doublets <- max(1L, as.integer(round(ncol(object) * expected_rate)))
  homotypic_fraction <- model_homotypic(object$doublet_preclusters)
  expected_heterotypic_doublets <- max(
    1L,
    as.integer(round(expected_doublets * (1 - homotypic_fraction)))
  )

  object <- call_doublets(
    object,
    PCs = df_dims,
    pN = doublet_pN,
    pK = best_pk,
    nExp = expected_heterotypic_doublets,
    reuse.pANN = formals(call_doublets)$reuse.pANN,
    sct = FALSE
  )
  classification_columns <- grep("^DF.classifications", colnames(object[[]]), value = TRUE)
  if (length(classification_columns) == 0L) {
    stop("DoubletFinder returned no classification for ", sample_id)
  }
  classification_column <- classification_columns[[length(classification_columns)]]
  object$doublet_call <- as.character(object[[]][[classification_column]])
  object$doublet_call[object$doublet_call != "Doublet"] <- "Singlet"

  list(
    object = object,
    summary = data.frame(
      sample_id = sample_id,
      cells_entered_doubletfinder = ncol(object),
      estimated_doublet_rate = expected_rate,
      selected_pK = best_pk,
      expected_heterotypic_doublets = expected_heterotypic_doublets,
      doublets_called = sum(object$doublet_call == "Doublet"),
      singlets_retained = sum(object$doublet_call == "Singlet"),
      selected_pcs = selected_pc,
      dimensions_used = paste0("1:", selected_pc),
      stringsAsFactors = FALSE
    )
  )
}

write_sigcell_keep_list <- function(object, output_file) {
  # 一行一个 barcode、无表头；readLines() 读取后可直接用于 subset()。
  writeLines(rownames(object[[]]), con = output_file)
}

sigcell_scdblfinder_processing <- function(e, dims, normalization_method, scale_factor) {
  if (!requireNamespace("Matrix", quietly = TRUE) ||
      !requireNamespace("rsvd", quietly = TRUE)) {
    stop("scDblFinder processing requires Matrix and rsvd.")
  }
  if (!identical(normalization_method, "LogNormalize")) {
    stop("sigcell_scdblfinder_processing currently supports normalization_method = `LogNormalize` only.")
  }
  library_size <- Matrix::colSums(e)
  if (any(!is.finite(library_size) | library_size <= 0)) {
    stop("scDblFinder requires positive library sizes for every cell.")
  }
  normalized <- e %*% Matrix::Diagonal(x = scale_factor / library_size)
  if (inherits(normalized, "sparseMatrix")) {
    normalized@x <- log1p(normalized@x)
  } else {
    normalized <- log1p(normalized)
  }
  k <- min(as.integer(dims), nrow(normalized) - 1L, ncol(normalized) - 1L)
  if (k < 2L) stop("scDblFinder processing needs at least two PCA dimensions.")
  pca <- rsvd::rpca(Matrix::t(normalized), k = k, center = TRUE, scale = FALSE)$x
  rownames(pca) <- colnames(e)
  pca
}

estimate_sigcell_scdblfinder_rate <- function(n_cells, rate_per_1000, max_rate) {
  min(max_rate, rate_per_1000 * n_cells / 1000)
}

sample_sigcell_seed <- function(sample_id, random_seed) {
  code_points <- utf8ToInt(sample_id)
  as.integer((random_seed + sum(code_points * seq_along(code_points))) %% .Machine$integer.max)
}

run_sigcell_scdblfinder <- function(
  counts, sample_id, n_variable_features, doublet_dims, pca_plot_dims,
  precluster_resolution, doublet_rate_per_1000, doublet_rate_max, random_seed,
  normalization_method, scale_factor, variable_feature_method
) {
  if (ncol(counts) < 100L) {
    stop(sample_id, " has fewer than 100 QC-passed cells; scDblFinder is not suitable.")
  }
  clustering_object <- Seurat::CreateSeuratObject(
    counts = counts, project = sample_id, min.cells = 0, min.features = 0
  )
  clustering_object <- Seurat::NormalizeData(
    clustering_object, normalization.method = normalization_method,
    scale.factor = scale_factor, verbose = FALSE
  )
  clustering_object <- Seurat::FindVariableFeatures(
    clustering_object, selection.method = variable_feature_method,
    nfeatures = min(n_variable_features, nrow(clustering_object)), verbose = FALSE
  )
  clustering_object <- Seurat::ScaleData(
    clustering_object, features = Seurat::VariableFeatures(clustering_object), verbose = FALSE
  )
  clustering_object <- Seurat::RunPCA(
    clustering_object, features = Seurat::VariableFeatures(clustering_object),
    npcs = min(pca_plot_dims, ncol(clustering_object) - 1L), verbose = FALSE
  )
  available_pcs <- ncol(Seurat::Embeddings(clustering_object, "pca"))
  clustering_dims <- seq_len(min(doublet_dims, available_pcs))
  if (length(clustering_dims) < 2L) stop(sample_id, " does not have enough PCA dimensions for clustering.")
  clustering_object <- Seurat::FindNeighbors(
    clustering_object, reduction = "pca", dims = clustering_dims, verbose = FALSE
  )
  clustering_object <- Seurat::FindClusters(
    clustering_object, resolution = precluster_resolution,
    cluster.name = "scDblFinder_preclusters", verbose = FALSE
  )
  clusters <- clustering_object$scDblFinder_preclusters
  names(clusters) <- colnames(clustering_object)
  pca_stdev <- clustering_object[["pca"]]@stdev
  rm(clustering_object)
  gc(verbose = FALSE)

  precluster_mode <- "cluster_based"
  n_preclusters <- length(unique(as.character(clusters)))
  if (n_preclusters < 2L) {
    warning("Only one precluster was found for ", sample_id, "; using random artificial doublets.")
    clusters <- NULL
    precluster_mode <- "random_fallback"
  }
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  dbr <- estimate_sigcell_scdblfinder_rate(ncol(sce), doublet_rate_per_1000, doublet_rate_max)
  bp <- BiocParallel::SerialParam(RNGseed = sample_sigcell_seed(sample_id, random_seed))
  processing <- function(e, dims) {
    sigcell_scdblfinder_processing(e, dims, normalization_method, scale_factor)
  }
  result <- scDblFinder::scDblFinder(
    sce,
    clusters = if (is.null(clusters)) NULL else clusters[colnames(sce)],
    dbr = dbr,
    nfeatures = min(n_variable_features, nrow(sce)),
    dims = min(doublet_dims, ncol(sce) - 1L),
    processing = processing,
    BPPARAM = bp,
    verbose = FALSE
  )
  classification <- as.character(result$scDblFinder.class)
  score <- as.numeric(result$scDblFinder.score)
  if (length(classification) != ncol(sce) || length(score) != ncol(sce) || anyNA(classification)) {
    stop("scDblFinder returned an invalid result for ", sample_id)
  }
  list(
    dbr = dbr,
    clusters = clusters,
    classification = classification,
    score = score,
    pca_stdev = pca_stdev,
    clustering_dims = clustering_dims,
    n_preclusters = n_preclusters,
    precluster_mode = precluster_mode
  )
}
