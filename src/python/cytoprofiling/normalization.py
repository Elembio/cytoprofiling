
import numpy as np
import pandas as pd
from typing import Callable
from cytoprofiling.utility import get_all_targets, get_barcoding_batches, get_wells, get_default_normalization_targets


def normalize_cells_by_aggregated_counts(df : pd.DataFrame, batch_names : list[str] = None, well_names : list[str] = None, aggregation_func : Callable = np.nansum , normalization_targets : list[str] = []) -> pd.DataFrame:
    """Normalize counts in each cell by some aggregated counts per cell

    Normalized counts are rescaled to approximately maintain the original counts per cell
    in the raw data

    Args:
      df : Input dataframe with counts
      batch_names : List of batches to normalize
      well_names : List of wells to include in normalization
      aggregation_func : Function to use for counts aggregation (e.g., nansum or nanmedian), must accept axis arg
      batch_targets: List of target subset to use for calculating of aggregated counts

    Returns:
      Dataframe with normalized counts
    """
    if well_names is None:
        well_names = get_wells(df)
    if batch_names is None:
        batch_names = get_barcoding_batches(df)
    result = df.copy()
    for batch_name in batch_names:
      all_batch_targets = get_default_normalization_targets(df, batch_name)
      
      if len(normalization_targets) > 0:
          batch_norm_targets = []
          for target in normalization_targets:
              if f"{target}.{batch_name}" in all_batch_targets:
                  batch_norm_targets.append(f"{target}.{batch_name}")
      else:
          batch_norm_targets = all_batch_targets

      median_counts_by_well = []
      for well_name in well_names:
          median_counts_by_well.append(np.nanmedian(aggregation_func(df.loc[df["Well"] == well_name, batch_norm_targets], axis = 1)))
      
      #
      # this is the target median counts per cell for all wells
      #
      target_median = np.nanmedian(median_counts_by_well)

      norm_values = aggregation_func(df[batch_norm_targets], axis = 1)
      all_correction_targets = get_all_targets(df, batch_name)
      
      for correction_target in all_correction_targets:
          result[correction_target] = target_median * np.divide(df[correction_target], norm_values)

    return result

def normalize_wells_by_median_of_ratios(df : pd.DataFrame, batch_names : list[str] = None, well_names : list[str] = None) -> pd.DataFrame:
    """Normalize counts in each well by median of ratios method

    Args:
      df : Input dataframe with counts
      batch_names : List of batches to normalize
      well_names : List of wells to include in normalization

    Returns:
      Dataframe with normalized counts
    """
    if well_names is None:
        well_names = get_wells(df)
    if batch_names is None:
        batch_names = get_barcoding_batches(df)
    result = df.copy()
    for batch_name in batch_names:
      #
      # get data and targets for normalization
      #
      batch_targets = get_default_normalization_targets(df, batch_name)
      all_wells_df = df.loc[df["Well"].isin(well_names), ["Well", ] + batch_targets]

      #
      # Aggregate mean counts by well
      #
      grouped_all_wells_df = all_wells_df.groupby("Well").mean()

      #
      # log transformation of mean counts
      #
      grouped_all_wells_df[batch_targets] = np.log(grouped_all_wells_df[batch_targets])

      #
      # Generate the log transformed counts for the psuedo-reference sample
      #
      psuedo_reference = np.nanmean(grouped_all_wells_df, axis = 0)

      #
      # Scale counts by the psuedo-reference
      #
      grouped_all_wells_df[batch_targets] -= psuedo_reference

      #
      # Find the scaling factor for each well
      #
      well_log_ratios = np.nanmedian(grouped_all_wells_df, axis = 1)
      scaling_factors = np.exp(well_log_ratios)

      #
      # Apply the scaling factor
      #
      all_targets = get_all_targets(df, batch_name)
      

      for well_idx in range(len(grouped_all_wells_df.index)):
          well_name = grouped_all_wells_df.index[well_idx]
          result.loc[result["Well"] == well_name, all_targets] =  result.loc[result["Well"] == well_name, all_targets] / scaling_factors[well_idx]
    
    return result

def normalize_cytoprofiling(df : pd.DataFrame, batch_names : list[str] = None, well_names : list[str] = None) -> pd.DataFrame:
    """Perform default normalization for cells and samples

    Normalize cells by assigned counts per batch and normalize samples by median of ratios.

    Args:
      df : Input dataframe with counts
      batch_names : List of batch names to normalize (default is to normalize all barcoding batches)
      well_names : List of wells to include in normalization (default is all wells)

    Returns:
      Dataframe with normalized counts
    """
    if well_names is None:
        well_names = get_wells(df)
    if batch_names is None:
        batch_names = get_barcoding_batches(df)
    
    return normalize_wells_by_median_of_ratios(
        normalize_cells_by_aggregated_counts(df, batch_names, well_names),
        batch_names,
        well_names
    )