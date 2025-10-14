library(Seurat)

#' Export cytoprofiling cell table to Seurat object
#'
#' @param df Cell table (see load_cytoprofiling)
#' @param panel_json JSON data for Panel.json (see load_panel)
#' @param use_gene_names If TRUE, use gene names as rownames, otherwise use target names. Duplicate gene names will be removed.
#' @param split_assays If TRUE, split the transcript and protein barcoding targets into separate Seurat assays ("Transcript" and "Protein", respectively).
#' @return Seurat object
#' @export
cytoprofiling_to_seurat <- function(df, panel_json, use_gene_names = FALSE, split_assays = FALSE) {
  barcoding_targets <- panel_json$BarcodingTargets

  if (split_assays) {
    transcript_targets <- barcoding_targets[sapply(barcoding_targets, function(x) x$TargetType == "Transcript")]
  } else {
    transcript_targets <- barcoding_targets
  }

  transcript_features <- sapply(transcript_targets, function(x) {
    paste(x$Target, x$BatchName, sep = ".")
  })

  if (use_gene_names) {
    transcript_genes <- lapply(strsplit(transcript_features, ".", fixed = TRUE), function(x) {
      x[[1]]
    })
    transcript_features <- transcript_features[!duplicated(transcript_genes)]
  }

  transcript_counts <- t(df[, transcript_features])
  colnames(transcript_counts) <- df[, "Cell"]

  assay.name <- if (split_assays) "RNA" else "Barcoding"

  result <- Seurat::CreateSeuratObject(transcript_counts, assay = assay.name)
  result <- Seurat::AddMetaData(result, as.factor(df[, c("Well")]), col.name = "Well")
  result <- Seurat::AddMetaData(result, as.factor(df[, c("WellLabel")]), col.name = "WellLabel")
  Seurat::Idents(result) <- result$WellLabel

  if (use_gene_names) {
    rownames(result[[assay.name]]) <- sapply(rownames(result[[assay.name]]), function(x) {
      strsplit(x, ".", fixed = TRUE)[[1]][1]
    })
  }

  if (split_assays) {
    protein_targets <- barcoding_targets[sapply(barcoding_targets, function(x) x$TargetType == "Protein")]
    protein_features <- sapply(protein_targets, function(x) {
      paste(x$Target, x$BatchName, sep = ".")
    })

    if (use_gene_names) {
      protein_genes <- lapply(strsplit(protein_features, ".", fixed = TRUE), function(x) {
        x[[1]]
      })
      protein_features <- protein_features[!duplicated(protein_genes)]
    }

    protein_counts <- t(df[, protein_features])
    colnames(protein_counts) <- df[, "Cell"]

    result[["Protein"]] <- Seurat::CreateAssayObject(counts = protein_counts, assay = "Protein")

    if (use_gene_names) {
      rownames(result[["Protein"]]) <- sapply(rownames(result[["Protein"]]), function(x) {
        strsplit(x, ".", fixed = TRUE)[[1]][1]
      })
    }
  }

  return(result)
}
