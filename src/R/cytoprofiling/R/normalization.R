#' Normalized cell table by aggregated counts
#'
#' Normalize each cell by number of counts in each batch. If list of normalization targets is
#' provided, normalize by the sum of counts in the targets. Otherwise, normalize by the sum of
#' default normalization targets (see get_default_normalization_targets).
#'
#' @param df Cell table
#' @param batch_names List of batch names
#' @param normalization_targets List of normalization targets
#' @return Normalized cell table
#' @export
normalize_cells_by_aggregated_counts <- function(df, batch_names = c(), normalization_targets = c()) {
  result <- df
  if (length(batch_names) == 0) {
    batch_names <- get_barcoding_batches(result)
  }
  for (batch in batch_names) {
    if (length(normalization_targets) == 0) {
      normalization_targets <- get_default_normalization_targets(result, batch)
    }
    normalization_targets <- get_default_normalization_targets(result, batch)
    norm_values <- rowSums(result[, normalization_targets], na.rm = TRUE)
    for (target in get_all_targets(result, batch)) {
      result[target] <- result[target] / norm_values
    }
  }
  return(result)
}

#' Normalize samples by median of ratios method
#'
#' Normalize cytoprofiling cell table such that most targets have smallest absolute log ratio
#' across provided batches and wells.
#'
#' @param df Cell table
#' @param batch_names List of batch names
#' @param well_names List of well names
#' @return Normalized cell table
#' @export
normalize_wells_by_median_of_ratios <- function(df, batch_names = c(), well_names = c()) {
  if (length(well_names) == 0) {
    well_names <- unique(df$Well)
  }

  if (length(batch_names) == 0) {
    batch_names <- get_barcoding_batches(df)
  }

  # filter df to only include wells in well_names
  result <- df[df$Well %in% well_names, ]
  for (batch_name in batch_names) {
    normalization_targets <- get_default_normalization_targets(df, batch_name)

    # aggregate by mean count per well, ignoring NAs
    well_counts <- aggregate(result[normalization_targets], by = list(result$Well), FUN = mean, na.rm = TRUE)

    # log transform well_counts
    well_counts[, -1] <- log(well_counts[, -1])

    # center each column of well_counts by mean
    well_counts[, -1] <- sweep(well_counts[, -1], 2, colMeans(well_counts[, -1], na.rm = TRUE))

    # find the median of each row of well_counts
    median_counts <- apply(well_counts[, -1], 1, median, na.rm = TRUE)

    # calculate the scaling factor as the expontential of the median log ratio
    scaling_factors <- exp(median_counts)

    # for each group of rows for each well in the original df, divide by the scaling factor
    all_targets <- get_all_targets(result, batch_name)
    for (well_idx in seq_along(well_names)) {
      well <- well_names[well_idx]
      well_rows <- df$Well == well
      result[well_rows, all_targets] <- result[well_rows, all_targets] / scaling_factors[well_idx]
    }
  }
  return(result)
}

#' Perform default normalization for cells and samples
#'
#' Normalize cells by assigned counts per batch and normalize samples by median of ratios.
#'
#' @param df Cell table
#' @param batch_names List of batch names
#' @param well_names List of well names
#' @return Normalized cell table
#' @export
normalize_cytoprofiling <- function(df, batch_names = c(), well_names = c()) {
  return(normalize_wells_by_median_of_ratios(
    normalize_cells_by_aggregated_counts(
      df,
      batch_names = batch_names
    ),
    batch_names = batch_names,
    well_names = well_names
  ))
}
