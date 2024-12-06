#' Export cytoprofiling cell table to Seurat object
#'
#' Export cytoprofiling cell table to Seurat object
#'
#' @param df Cell table
#' @return Seurat object
#' @export
cytoprofiling_to_seurat <- function(df) {
  min.features <- 5
  min.cells <- 3

  batches <- get_barcoding_batches(df)
  all_targets <- c()
  for (batch in batches) {
    all_targets <- c(all_targets, get_default_normalization_targets(df, batch))
  }

  genes <- lapply(strsplit(all_targets, ".", fixed = TRUE), function(x) {
    x[[1]]
  })

  all_targets <- all_targets[!duplicated(genes)]

  result <- Seurat::CreateSeuratObject(t(df[, all_targets]), min.features = min.features, min.cells = min.cells)
  result <- Seurat::AddMetaData(result, as.factor(df[, c("Well")]), col.name = "Well")
  result <- Seurat::AddMetaData(result, as.factor(df[, c("WellLabel")]), col.name = "WellLabel")
  Seurat::Idents(result) <- result$WellLabel
  for (batch in batches) {
    rownames(result[["RNA"]]) <- gsub(paste("\\.", batch, sep = ""), "", rownames(result[["RNA"]]))
  }
  return(result)
}
