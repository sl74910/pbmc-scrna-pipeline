# Small, explicit helpers shared by the PBMC scripts.
# Every argument is required so that a changed analysis parameter is visible.
# Main contracts: sample_files -> sorted file table; read_cells -> barcodes;
# read_seurat_counts -> a Seurat object containing only requested cells;
# preprocess -> normalized/scaled PCA object; merge_samples -> merged object;
# qc_filter -> QC-passed object; save_* -> a written PDF; make_*_plots -> PDFs.

pbmc_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "))
  }
  invisible(packages)
}

pbmc_validate_parameters <- function(parameters) {
  if (!is.list(parameters) || is.null(names(parameters))) {
    stop("parameters must be a named list.")
  }
  unset <- names(parameters)[vapply(parameters, function(value) {
    length(value) == 0L || anyNA(value)
  }, logical(1))]
  if (length(unset) > 0L) {
    stop("Parameter(s) are NA or empty: ", paste(unset, collapse = ", "))
  }
  invisible(parameters)
}

pbmc_make_dirs <- function(paths) {
  if (!is.character(paths) || anyNA(paths) || any(!nzchar(paths))) {
    stop("paths must contain non-empty directory names.")
  }
  for (path in paths) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path) || file.access(path, 2L) != 0) {
      stop("Directory is not writable: ", path)
    }
  }
  invisible(paths)
}

pbmc_apply_samples <- function(indices, FUN, workers, packages, seed) {
  if (length(indices) == 0L) return(list())
  if (workers <= 1L) return(lapply(indices, FUN))
  pbmc_require_packages(c("future", "future.apply", packages))
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = min(as.integer(workers), length(indices)))
  future.apply::future_lapply(indices, FUN, future.seed = seed, future.packages = packages)
}

pbmc_sample_files <- function(input_dir, pattern) {
  if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No input files found in: ", input_dir)

  filename <- basename(files)
  sample_id <- sub("_Seurat\\.rds$", "", filename)
  day <- suppressWarnings(as.integer(sub("^PBMC_([0-9]+)day_.*$", "\\1", sample_id)))
  measurement <- suppressWarnings(as.integer(sub("^PBMC_[0-9]+day_([0-9]+)$", "\\1", sample_id)))
  order_index <- order(is.na(day), day, is.na(measurement), measurement, sample_id)
  data.frame(
    file = files[order_index],
    filename = filename[order_index],
    sample_id = sample_id[order_index],
    day = day[order_index],
    measurement = measurement[order_index],
    stringsAsFactors = FALSE
  )
}

pbmc_make_manifest <- function(input_dir, keep_dir, days, measurements, mouse_id, keep_suffix) {
  if (!length(days) || !length(measurements) || anyNA(days) || anyNA(measurements)) {
    stop("days and measurements cannot be empty or NA.")
  }
  if (length(mouse_id) != 1L || is.na(mouse_id) || !nzchar(mouse_id)) {
    stop("mouse_id must be one non-empty string.")
  }
  day_levels <- paste0("Day", days)
  grid <- expand.grid(
    day = days,
    measurement = measurements,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$day, grid$measurement), , drop = FALSE]
  sample_id <- paste0("PBMC_", grid$day, "day_", grid$measurement)
  data.frame(
    file = file.path(input_dir, paste0(sample_id, "_Seurat.rds")),
    filename = paste0(sample_id, "_Seurat.rds"),
    sample_id = sample_id,
    keep_file = file.path(keep_dir, paste0(sample_id, keep_suffix)),
    day = grid$day,
    day_group = factor(paste0("Day", grid$day), levels = day_levels),
    mouse_id = mouse_id,
    mouse_timepoint_id = paste0(mouse_id, "_Day", grid$day),
    measurement = grid$measurement,
    measurement_id = paste0("R", grid$measurement),
    stringsAsFactors = FALSE
  )
}

pbmc_check_manifest <- function(manifest, check_keep_files) {
  missing_inputs <- manifest$filename[!file.exists(manifest$file)]
  if (length(missing_inputs)) {
    stop("Missing input RDS: ", paste(missing_inputs, collapse = ", "))
  }
  if (isTRUE(check_keep_files)) {
    missing_keep <- basename(manifest$keep_file[!file.exists(manifest$keep_file)])
    if (length(missing_keep)) {
      stop("Missing keep lists: ", paste(missing_keep, collapse = ", "))
    }
  }
  invisible(manifest)
}

pbmc_read_cells <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
  cells <- trimws(readLines(path, warn = FALSE))
  cells <- unique(cells[nzchar(cells)])
  if (!length(cells)) stop("No cell barcodes found in ", label, ": ", path)
  cells
}

pbmc_read_seurat_counts <- function(input_file, keep_cells, sample_id, metadata) {
  object <- readRDS(input_file)
  if (!inherits(object, "Seurat")) stop("Input is not a Seurat object: ", input_file)
  if (!"RNA" %in% SeuratObject::Assays(object) ||
      !"counts" %in% SeuratObject::Layers(object[["RNA"]])) {
    stop("Input requires an RNA counts layer: ", input_file)
  }
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  keep <- intersect(unique(keep_cells), colnames(counts))
  if (!length(keep)) stop("No requested barcodes match ", sample_id)
  result <- Seurat::CreateSeuratObject(
    counts = counts[, keep, drop = FALSE],
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
  if (!is.null(metadata)) {
    if (!is.list(metadata) || is.null(names(metadata))) stop("metadata must be a named list.")
    for (name in names(metadata)) result[[name]] <- metadata[[name]]
  }
  result
}

pbmc_add_manifest_metadata <- function(object, manifest_row, day_levels) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.")
  n <- ncol(object)
  object$sample_id <- rep(as.character(manifest_row$sample_id[[1L]]), n)
  object$tissue <- rep("PBMC", n)
  object$day <- rep(manifest_row$day[[1L]], n)
  object$day_group <- factor(rep(as.character(manifest_row$day_group[[1L]]), n), levels = day_levels)
  object$measurement <- rep(manifest_row$measurement[[1L]], n)
  object$measurement_id <- rep(as.character(manifest_row$measurement_id[[1L]]), n)
  object$mouse_id <- rep(as.character(manifest_row$mouse_id[[1L]]), n)
  object$mouse_timepoint_id <- rep(as.character(manifest_row$mouse_timepoint_id[[1L]]), n)
  object
}

pbmc_merge_samples <- function(objects, project, split_column) {
  if (!is.list(objects) || !length(objects)) stop("objects cannot be empty.")
  if (length(objects) == 1L) result <- objects[[1L]] else {
    result <- merge(
      x = objects[[1L]],
      y = objects[-1L],
      add.cell.ids = names(objects),
      project = project
    )
  }
  if (length(unique(as.character(result[[split_column]][, 1L]))) > 1L) {
    count_layers <- grep("^counts", SeuratObject::Layers(result[["RNA"]]), value = TRUE)
    if (length(count_layers) == 1L) {
      result[["RNA"]] <- split(result[["RNA"]], f = result[[split_column]][, 1L])
    }
  }
  result
}

pbmc_preprocess <- function(object, normalization_method, scale_factor,
                            variable_method, n_variable_features, n_pcs, random_seed) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.")
  if (ncol(object) < 3L) stop("At least three cells are needed for PCA.")
  n_variable_features <- min(as.integer(n_variable_features), nrow(object))
  n_pcs <- min(as.integer(n_pcs), ncol(object) - 1L, nrow(object) - 1L)
  if (n_variable_features < 2L || n_pcs < 2L) stop("Not enough genes/cells for preprocessing.")
  object <- Seurat::NormalizeData(object, normalization.method = normalization_method,
                                  scale.factor = scale_factor, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, selection.method = variable_method,
                                         nfeatures = n_variable_features, verbose = FALSE)
  object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
  Seurat::RunPCA(object, features = Seurat::VariableFeatures(object), npcs = n_pcs,
                 seed.use = random_seed, verbose = FALSE)
}

pbmc_log_normalize_features <- function(object, assay, features, normalization_method, scale_factor) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.")
  if (!identical(normalization_method, "LogNormalize")) {
    stop("pbmc_log_normalize_features currently supports normalization_method = `LogNormalize` only.")
  }
  counts <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  features <- intersect(features, rownames(counts))
  if (!length(features)) stop("No requested features are present in the counts layer.")
  library_size <- Matrix::colSums(counts)
  if (any(!is.finite(library_size) | library_size <= 0)) {
    stop("Log-normalization requires positive library sizes for every cell.")
  }
  normalized <- counts[features, , drop = FALSE] %*%
    Matrix::Diagonal(x = scale_factor / library_size)
  if (inherits(normalized, "sparseMatrix")) {
    normalized@x <- log1p(normalized@x)
  } else {
    normalized <- log1p(normalized)
  }
  normalized
}

pbmc_add_qc_metadata <- function(object, mitochondrial_pattern, hemoglobin_genes) {
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object.")
  object[["percent_mt"]] <- Seurat::PercentageFeatureSet(object, pattern = mitochondrial_pattern)
  object[["percent_hb"]] <- Seurat::PercentageFeatureSet(object, features = hemoglobin_genes)
  object
}

pbmc_filter_qc <- function(object, min_features, max_features, max_counts,
                           max_percent_mt, max_percent_hb) {
  metadata <- object[[]]
  keep <- metadata$nFeature_RNA >= min_features &
    metadata$nFeature_RNA <= max_features &
    metadata$nCount_RNA <= max_counts &
    metadata$percent_mt <= max_percent_mt &
    metadata$percent_hb <= max_percent_hb
  keep[is.na(keep)] <- FALSE
  result <- object[, keep]
  if (!ncol(result)) stop("QC removed every cell.")
  result
}

pbmc_qc_plot <- function(object, sample_name, thresholds, point_fraction) {
  required <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_hb")
  missing <- setdiff(required, colnames(object[[]]))
  if (length(missing)) stop("Missing QC metadata: ", paste(missing, collapse = ", "))
  if (!is.list(thresholds) || !all(c("min_features", "max_features", "max_counts",
                                    "max_percent_mt", "max_percent_hb") %in% names(thresholds))) {
    stop("thresholds must contain max_features, max_counts, max_percent_mt and max_percent_hb.")
  }
  metadata <- object[[]]
  values <- do.call(rbind, lapply(required, function(feature) {
    data.frame(feature = feature, value = as.numeric(metadata[[feature]]), sample = sample_name)
  }))
  values$feature <- factor(values$feature, levels = required)
  points <- do.call(rbind, lapply(split(values, values$feature), function(data) {
    n <- floor(nrow(data) * point_fraction)
    if (n == 0L) data[FALSE, , drop = FALSE] else data[sample.int(nrow(data), n), , drop = FALSE]
  }))
  limits <- data.frame(
    feature = factor(c("nFeature_RNA", "nFeature_RNA", "nCount_RNA", "percent_mt", "percent_hb"), levels = required),
    value = c(thresholds$min_features, thresholds$max_features, thresholds$max_counts,
              thresholds$max_percent_mt, thresholds$max_percent_hb)
  )
  ggplot2::ggplot(values, ggplot2::aes(sample, value)) +
    ggplot2::geom_point(data = points, position = ggplot2::position_jitter(width = 0.12),
                        size = 0.55, alpha = 0.34, colour = "#1F5F8B") +
    ggplot2::geom_hline(data = limits, ggplot2::aes(yintercept = value), inherit.aes = FALSE,
                        colour = "#C96B6B", linetype = "dashed", linewidth = 0.45) +
    ggplot2::geom_boxplot(width = 0.30, outlier.shape = NA, fill = "white",
                          colour = "#123B5D", linewidth = 0.45) +
    ggplot2::stat_summary(fun = median, geom = "point", shape = 23, size = 1.8,
                          fill = "#D97706", colour = "#123B5D", stroke = 0.35) +
    ggplot2::facet_wrap(~feature, ncol = 4, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL, title = sub("_Seurat\\.rds$", "", sample_name)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey88"),
                   plot.title = ggplot2::element_text(hjust = 0.5))
}

pbmc_save_pdf <- function(plot, file, width, height) {
  grDevices::pdf(file, width = width, height = height, compress = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(file)
}

pbmc_save_plot_grid <- function(plots, file, ncol, nrow, width, height) {
  if (!length(plots)) stop("plots cannot be empty.")
  pbmc_save_pdf(patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow, byrow = TRUE),
                file, width, height)
}

pbmc_save_elbow_plot <- function(object, selected_pcs, file, n_pcs) {
  available <- ncol(Seurat::Embeddings(object, "pca"))
  pbmc_save_pdf(Seurat::ElbowPlot(object, ndims = min(n_pcs, available)) +
                  ggplot2::geom_vline(xintercept = selected_pcs, colour = "red", linetype = "dashed"),
                file, 7, 5)
}

pbmc_save_scdblfinder_elbow <- function(stdev, selected_pcs, sample_id, file, n_pcs) {
  stdev <- stdev[seq_len(min(as.integer(n_pcs), length(stdev)))]
  pbmc_save_pdf({
    plot(seq_along(stdev), stdev, type = "b", pch = 16, cex = 0.55,
         xlab = "PC", ylab = "Standard deviation", main = sample_id)
    abline(v = selected_pcs, col = "red", lty = 2)
  }, file, 7, 5)
  invisible(file)
}

pbmc_require_metadata <- function(object, columns) {
  missing <- setdiff(columns, colnames(object[[]]))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "))
  invisible(object)
}

pbmc_require_reductions <- function(object, reductions) {
  missing <- setdiff(reductions, Seurat::Reductions(object))
  if (length(missing)) stop("Missing reductions: ", paste(missing, collapse = ", "))
  invisible(object)
}

pbmc_load_rdata_object <- function(file, object_name) {
  if (!file.exists(file)) stop("Missing RData file: ", file)
  environment <- new.env(parent = emptyenv())
  load(file, envir = environment)
  if (!exists(object_name, envir = environment, inherits = FALSE)) {
    stop("RData does not contain object `", object_name, "`.")
  }
  environment[[object_name]]
}

pbmc_order_day_values <- function(values, preferred_levels) {
  observed <- unique(as.character(values))
  if (anyNA(observed)) stop("Day values cannot be NA.")
  ordered <- c(intersect(preferred_levels, observed), setdiff(observed, preferred_levels))
  if (!length(ordered)) stop("No day values available.")
  ordered
}

pbmc_order_sample_levels <- function(object, sample_column, day_column) {
  metadata <- unique(object[[]][, c(sample_column, day_column), drop = FALSE])
  metadata[[sample_column]] <- as.character(metadata[[sample_column]])
  metadata[[day_column]] <- as.character(metadata[[day_column]])
  metadata$day_number <- suppressWarnings(as.integer(sub("^Day", "", metadata[[day_column]])))
  metadata$batch_number <- suppressWarnings(as.integer(sub(".*_(\\d+)$", "\\1", metadata[[sample_column]])))
  metadata <- metadata[order(is.na(metadata$day_number), metadata$day_number,
                             is.na(metadata$batch_number), metadata$batch_number,
                             metadata[[sample_column]]), , drop = FALSE]
  metadata[[sample_column]]
}

pbmc_make_umap_panels <- function(object, reduction, cluster_column, sample_column, day_column) {
  Seurat::DimPlot(object, reduction = reduction, group.by = sample_column) /
    Seurat::DimPlot(object, reduction = reduction, group.by = day_column) /
    Seurat::DimPlot(object, reduction = reduction, group.by = cluster_column, label = TRUE)
}

pbmc_make_annotation_plots <- function(object, reduction, label_column, figure_dir,
                                       day_column, sample_column, preferred_days, prefix) {
  pbmc_require_reductions(object, reduction)
  pbmc_require_metadata(object, c(label_column, day_column, sample_column))
  day_levels <- pbmc_order_day_values(object[[day_column]][, 1L], preferred_days)
  object[[day_column]] <- factor(as.character(object[[day_column]][, 1L]), levels = day_levels)
  object[[sample_column]] <- factor(as.character(object[[sample_column]][, 1L]),
                                    levels = pbmc_order_sample_levels(object, sample_column, day_column))
  by_day <- Seurat::DimPlot(object, reduction = reduction, group.by = label_column,
                            split.by = day_column, ncol = 1, combine = TRUE, raster = TRUE)
  by_sample <- Seurat::DimPlot(object, reduction = reduction, group.by = label_column,
                               split.by = sample_column, ncol = 5, combine = TRUE, raster = TRUE)
  pbmc_save_pdf(by_day, file.path(figure_dir, paste0(prefix, "_UMAP_by_day.pdf")), 8, 22)
  pbmc_save_pdf(by_sample, file.path(figure_dir, paste0(prefix, "_UMAP_by_batch_day.pdf")), 22, 16)
  invisible(object)
}

pbmc_score_modules <- function(expression, modules) {
  score <- sapply(modules, function(genes) {
    genes <- intersect(genes, rownames(expression))
    if (!length(genes)) rep(NA_real_, ncol(expression)) else colMeans(expression[genes, , drop = FALSE])
  })
  score <- as.matrix(score)
  rownames(score) <- colnames(expression)
  score
}

pbmc_best_module <- function(score_matrix) {
  if (!is.matrix(score_matrix) || ncol(score_matrix) == 0L) {
    stop("score_matrix must have at least one module column.")
  }
  safe <- score_matrix
  safe[is.na(safe)] <- -Inf
  best <- colnames(safe)[max.col(safe, ties.method = "first")]
  sorted <- t(apply(safe, 1L, sort, decreasing = TRUE))
  margin <- if (ncol(safe) == 1L) rep(0, nrow(safe)) else sorted[, 1L] - sorted[, 2L]
  margin[!is.finite(margin)] <- NA_real_
  score <- sorted[, 1L]
  score[is.infinite(score)] <- NA_real_
  data.frame(candidate = best, score = score, margin = margin,
             row.names = rownames(score_matrix), stringsAsFactors = FALSE)
}
