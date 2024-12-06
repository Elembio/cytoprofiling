#' Filter cells to remove outliers by area
#'
#' @param df Cell table
#' @param well_names List of well names
#' @param min_threshold Lower percentile threshold
#' @param max_threshold Upper percentile threshold
#' @return Filtered cell table
#' @export
filter_cells_by_area_percentile <- function(df, well_names = c(), min_threshold = 3, max_threshold = 97) {
  if (length(well_names) == 0) {
    well_names <- unique(df$Well)
  }
  filter_status <- rep(FALSE, dim(df)[[1]])
  for (well in well_names) {
    thresholds <- quantile(df[df$Well == well, ]$Area, c(min_threshold / 100.0, max_threshold / 100.0), na.rm = TRUE)
    filter_status <- filter_status | (df$Well == well) & (df$Area > thresholds[[1]]) & (df$Area < thresholds[[2]])
  }
  return(df[filter_status, ])
}

#' Filter cells to remove outliers by assigned counts
#'
#' Cells must fall within percentile range for all batches to be kept
#'
#' @param df Cell table
#' @param batch_names List of batch names
#' @param well_names List of well names
#' @param min_threshold Lower percentile threshold
#' @param max_threshold Upper percentile threshold
#' @return Filtered cell table
#' @export
filter_cells_by_assigned_counts_percentile <- function(df, batch_names = c(), well_names = c(), min_threshold = 3, max_threshold = 97) {
  if (length(batch_names) == 0) {
    batch_names <- get_barcoding_batches(df)
  }
  if (length(well_names) == 0) {
    well_names <- unique(df$Well)
  }

  final_filter_status <- rep(TRUE, dim(df)[[1]])

  for (batch in batch_names) {
    normalization_targets <- get_default_normalization_targets(df, batch)
    norm_values <- rowSums(df[, normalization_targets], na.rm = TRUE)

    current_filter_status <- rep(FALSE, dim(df)[[1]])
    for (well in well_names) {
      thresholds <- quantile(norm_values[(df$Well == well) & (norm_values > 0)], c(min_threshold / 100.0, max_threshold / 100.0), na.rm = TRUE)
      current_filter_status <- current_filter_status | (df$Well == well) & (norm_values > thresholds[[1]]) & (norm_values < thresholds[[2]])
    }
    final_filter_status <- final_filter_status & current_filter_status
  }
  return(df[final_filter_status, ])
}



#' Filter cells on assigned rate
#'
#' Cells must exceed the assignment rate threshold in all batches
#' to be kept
#'
#' @param df Cell table
#' @param batch_names List of batch names
#' @param well_names List of well names
#' @param assigned_rate_threshold Assigned rate threshold
#' @return Filtered cell table
#' @export
filter_cells_by_assigned_rate <- function(df, batch_names = c(), well_names = c(), assigned_rate_threshold = 0.5) {
  if (length(batch_names) == 0) {
    batch_names <- get_barcoding_batches(df)
  }
  if (length(well_names) == 0) {
    well_names <- unique(df$Well)
  }

  final_filter_status <- rep(TRUE, dim(df)[[1]])
  for (batch_name in batch_names) {
    batch_columns <- get_default_normalization_targets(df, batch_name)
    assigned_counts <- rowSums(df[, batch_columns], na.rm = TRUE)
    unassigned_counts <- df[[paste("Unassigned.", batch_name, sep = "")]]
    current_filter_status <- rep(FALSE, dim(df)[[1]])
    for (well_name in well_names) {
      current_filter_status <- current_filter_status | ((df$Well == well_name) & (assigned_counts > 0) & (assigned_counts / (unassigned_counts + assigned_counts) > assigned_rate_threshold))
    }
    final_filter_status <- final_filter_status & current_filter_status
  }
  return(df[final_filter_status, ])
}

#' Perform default filtering of cells (filter by area, assigned counts, and assigned rate)
#'
#' @param df Cell table
#' @return Filtered cell table
#' @export
filter_cells <- function(df, batch_names = c(), well_names = c()) {
  return(filter_cells_by_assigned_rate(
    filter_cells_by_area_percentile(
      filter_cells_by_assigned_counts_percentile(
        df,
        batch_names = batch_names,
        well_names = well_names
      ),
      well_names = well_names
    ),
    batch_names = batch_names,
    well_names = well_names
  ))
}
