library(arrow)
library(rjson)

#' Load a parquet file into a data frame
#'
#' Calls arrow::read_parquet to load cytoprofiling
#' parquet file into a dataframe
#'
#' @param input_filename Path to the parquet file
#' @return Data frame
#' @export
load_cytoprofiling <- function(input_filename) {
  return(as.data.frame(arrow::read_parquet(input_filename)))
}

#' Load a Panel.json file
#' 
#' Calls rjson::fromJSON to load a Panel.json file
#' 
#' @param input_filename Path to the Panel.json file
#' @return JSON data as a list
#' @export
load_panel <- function(input_filename) {
  return (rjson::fromJSON(file = input_filename))
}


#' Get a list of barcoding batches
#'
#' Get a list of barcoding batches present in the cell table
#'
#' @param df Cell table
#' @return List of barcoding batches
#' @export
get_barcoding_batches <- function(df) {
  return(unique(lapply(strsplit(colnames(df)[grep("\\.[BP]0", colnames(df), fixed = FALSE)], ".", fixed = TRUE), function(x) {
    x[[2]]
  })))
}

#' Get normalization targets for a batch
#'
#' Get a standard list of targets to use for normalization in a batch
#'
#' @param df Cell table
#' @param batch Batch name
#' @return List of normalization targets
#' @export
get_default_normalization_targets <- function(df, batch) {
  return(colnames(df)[endsWith(colnames(df), paste(".", batch, sep = "")) & (!grepl("_Nuclear.", colnames(df), fixed = TRUE)) & (!grepl("Unassigned", colnames(df), fixed = TRUE))])
}

#' Get all targets for a batch
#'
#' Get all targets for a batch
#' @param df Cell table
#' @param batch Batch name
#' @return List of all targets
#' @export
get_all_targets <- function(df, batch) {
  return(colnames(df)[endsWith(colnames(df), paste(".", batch, sep = ""))])
}
