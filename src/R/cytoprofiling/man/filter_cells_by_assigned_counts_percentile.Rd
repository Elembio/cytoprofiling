% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering.R
\name{filter_cells_by_assigned_counts_percentile}
\alias{filter_cells_by_assigned_counts_percentile}
\title{Filter cells to remove outliers by assigned counts}
\usage{
filter_cells_by_assigned_counts_percentile(
  df,
  batch_names = c(),
  well_names = c(),
  min_threshold = 3,
  max_threshold = 97
)
}
\arguments{
\item{df}{Cell table}

\item{batch_names}{List of batch names}

\item{well_names}{List of well names}

\item{min_threshold}{Lower percentile threshold}

\item{max_threshold}{Upper percentile threshold}
}
\value{
Filtered cell table
}
\description{
Cells must fall within percentile range for all batches to be kept
}
