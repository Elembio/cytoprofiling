% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{normalize_wells_by_median_of_ratios}
\alias{normalize_wells_by_median_of_ratios}
\title{Normalize samples by median of ratios method}
\usage{
normalize_wells_by_median_of_ratios(df, batch_names = c(), well_names = c())
}
\arguments{
\item{df}{Cell table}

\item{batch_names}{List of batch names}

\item{well_names}{List of well names}
}
\value{
Normalized cell table
}
\description{
Normalize cytoprofiling cell table such that most targets have smallest absolute log ratio
across provided batches and wells.
}
