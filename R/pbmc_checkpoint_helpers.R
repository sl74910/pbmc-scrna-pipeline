# Compact RData checkpoints for the PBMC pipeline.

compact_pbmc_for_checkpoint <- function(
  object,
  layers,
  reductions
) {
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }

  if (!"RNA" %in% SeuratObject::Assays(object)) {
    stop("The checkpoint object must contain an RNA assay.")
  }

  available_layers <- SeuratObject::Layers(object[["RNA"]])
  requested_layers <- unlist(lapply(layers, function(layer_type) {
    SeuratObject::Layers(object[["RNA"]], search = layer_type)
  }), use.names = FALSE)
  requested_layers <- unique(requested_layers)

  if (length(requested_layers) == 0L) {
    stop(
      "None of the requested RNA layers are present. Requested: ",
      paste(layers, collapse = ", "), "; available: ",
      paste(available_layers, collapse = ", ")
    )
  }

  missing_reductions <- setdiff(reductions, Seurat::Reductions(object))
  if (length(missing_reductions) > 0L) {
    stop(
      "Requested reductions are absent: ",
      paste(missing_reductions, collapse = ", ")
    )
  }

  variable_features <- Seurat::VariableFeatures(object, assay = "RNA")
  compact_object <- Seurat::DietSeurat(
    object = object,
    assays = "RNA",
    layers = requested_layers,
    dimreducs = reductions,
    graphs = NULL,
    misc = FALSE
  )
  Seurat::VariableFeatures(compact_object, assay = "RNA") <- intersect(
    variable_features,
    rownames(compact_object)
  )
  if ("commands" %in% methods::slotNames(compact_object)) {
    compact_object@commands <- list()
  }

  compact_object
}

save_pbmc_rdata_checkpoint <- function(
  object,
  checkpoint_name,
  checkpoint_dir,
  analysis_parameters,
  layers,
  reductions,
  extra_metadata
) {
  if (!is.character(checkpoint_name) || length(checkpoint_name) != 1L ||
      !nzchar(checkpoint_name)) {
    stop("checkpoint_name must be one non-empty string.")
  }

  if (!is.list(analysis_parameters) || !is.list(extra_metadata)) {
    stop("analysis_parameters and extra_metadata must both be lists.")
  }

  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  pbmc <- compact_pbmc_for_checkpoint(
    object = object,
    layers = layers,
    reductions = reductions
  )
  checkpoint_metadata <- c(
    list(
      checkpoint_name = checkpoint_name,
      saved_at = Sys.time(),
      analysis_parameters = analysis_parameters,
      retained_rna_layers = SeuratObject::Layers(pbmc[["RNA"]]),
      retained_reductions = Seurat::Reductions(pbmc)
    ),
    extra_metadata
  )

  rdata_file <- file.path(checkpoint_dir, paste0(checkpoint_name, ".RData"))
  rdata_tmp <- paste0(rdata_file, ".tmp")
  checkpoint_environment <- list2env(
    list(pbmc = pbmc, checkpoint_metadata = checkpoint_metadata),
    parent = emptyenv()
  )

  # xz keeps the three requested checkpoints as small as possible on disk.
  save(
    list = c("pbmc", "checkpoint_metadata"),
    file = rdata_tmp,
    envir = checkpoint_environment,
    compress = "xz"
  )
  if (!file.rename(rdata_tmp, rdata_file)) {
    stop("Could not finalize checkpoint: ", rdata_file)
  }

  message("已保存精简 RData 检查点：", rdata_file)
  invisible(rdata_file)
}
