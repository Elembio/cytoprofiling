import numpy as np
import pandas as pd
from cytoprofiling.utility import get_wells, get_barcoding_batches, get_default_normalization_targets

def down_sample_cells(df : pd.DataFrame, target_count : int) -> pd.DataFrame:
    """Down sample to fixed number of cells per well
    
    Args:
        df : Input dataframe with cells
        target_count : Target number of cells per well

    Returns:
        Down sampled data frame
    """

    keep = np.zeros(df.shape[0])
    for well in get_wells(df):
        num_cells = np.sum(df["Well"] == well)
        cell_target = min(num_cells, target_count)
        cell_keep = ([1,] * cell_target) + ([0,] * (num_cells - cell_target))
        np.random.shuffle(cell_keep)
        keep[df["Well"] == well] = cell_keep
    
    result = df.loc[keep == 1,:]
    result.reset_index(drop=True, inplace=True)
    return result


def filter_cells_by_area_percentile(df : pd.DataFrame, well_names : list[str] = None, min_percentile : float = 3, max_percentile : float = 97) -> pd.DataFrame:
    """Filter cells to remove outliers by area

    Args:
      df : Input dataframe with cells
      well_names : List of wells to include in filtering
      min_percentile : Percentile to determine minimum value for thresholding
      max_percentile : Percentile to determine maximum value for thresholding

    Returns:
      Filtered data frame
    """
    if well_names is None:
        well_names = get_wells(df)
    filter = np.zeros(df.shape[0], dtype=np.int8)
    for well_name in well_names:
        well_areas = df[df["Well"] == well_name]["Area"]
        min_threshold = np.nanpercentile(well_areas, min_percentile) - 0.001
        max_threshold = np.nanpercentile(well_areas, max_percentile) + 0.001

        filter[(df["Well"] == well_name)&(well_areas >= min_threshold)&(well_areas <= max_threshold)] = 1
    
    result = df.loc[filter > 0,:]
    result.reset_index(drop=True, inplace=True)
    return result


def filter_cells_by_assigned_counts_percentile(df : pd.DataFrame, batch_names : list[str] = None, well_names : list[str] = None, min_percentile : float = 3, max_percentile : float = 97) -> pd.DataFrame:
    """Filter cells to remove outliers by assigned counts in a batch

    Cells must fall within percentile range for all batches to be kept

    Args:
      df : Input dataframe with cells
      batch_name : List of batches to include in filtering
      well_names : List of wells to include in filtering
      min_percentile : Percentile to determine minimum value for thresholding
      max_percentile : Percentile to determine maximum value for thresholding

    Returns:
      Filtered data frame
    """
    if well_names is None:
        well_names = get_wells(df)
    if batch_names is None:
        batch_names = get_barcoding_batches(df)

    filters = []
    for batch_name in batch_names:
        batch_columns = get_default_normalization_targets(df, batch_name)
        assigned_counts = np.nansum(df[batch_columns], axis = 1)
        filter = np.zeros(df.shape[0], dtype=np.uint8)

        for well_name in well_names:
            well_counts = assigned_counts[(df["Well"] == well_name)&(assigned_counts > 0)]
            min_threshold = np.nanpercentile(well_counts, min_percentile) - 0.001
            max_threshold = np.nanpercentile(well_counts, max_percentile) + 0.001
            filter[(df["Well"] == well_name)&(assigned_counts > min_threshold)&(assigned_counts<max_threshold)] = 1
        filters.append(filter)

    # compute the intersection of all filters
    filter = np.ones(df.shape[0], dtype=np.uint8)
    for f in filters:
        filter = filter & f

    result = df.loc[filter > 0,:]
    result.reset_index(drop=True, inplace=True)
    return result

def filter_cells_by_assigned_rate(df : pd.DataFrame, batch_names : list[str] = None, well_names : list[str] = None, assigned_rate : float = 0.5) -> pd.DataFrame:
    """Filter cells on assigned rate

    Cells must exceed the assignment rate threshold in all batches
    to be kept

    Args:
      df : Input dataframe with cells
      batch_name : List of batches to use for filtering
      well_names : List of wells to use in filtering
      assigned_rate : Assigned rate to use for filtering cells

    Returns:
      Filtered data frame
    """
    if well_names is None:
        well_names = get_wells(df)
    if batch_names is None:
        batch_names = get_barcoding_batches(df)

    filters = []
    for batch_name in batch_names:
        batch_columns = get_default_normalization_targets(df, batch_name)
        assigned_counts = np.nansum(df[batch_columns], axis = 1)
        unassigned_counts = df[f"Unassigned.{batch_name}"]
        filter = np.zeros(df.shape[0], dtype=np.uint8)

        for well_name in well_names:
            filter[(df["Well"] == well_name)&(assigned_counts > 0)&(assigned_counts / ((unassigned_counts + assigned_counts)) > assigned_rate)] = 1
        filters.append(filter)

    # compute the intersection of all filters
    filter = np.ones(df.shape[0], dtype=np.uint8)
    for f in filters:
        filter = filter & f

    result = df.loc[filter > 0,:]

    result.reset_index(drop=True, inplace=True)
    return result

def filter_cells(df : pd.DataFrame, batch_names : list[str] = None, well_names : list[str] = None, stats = None) -> pd.DataFrame:
    """Perform default filtering of cells (filter by area, assigned counts, and assigned rate)

    Args:
      df : Input dataframe with cells
      batch_names : List of batches to use for normalization
      well_names : List of wells to include in normalization
      stats : Dictionary to store statistics about filtering under the key 'filter_cells'
      
    Returns:
      Filtered data frame
    """
    if well_names is None:
        well_names = get_wells(df)
    if batch_names is None:
        batch_names = get_barcoding_batches(df)
    if stats is None:
        stats = {}
    
    result = filter_cells_by_area_percentile(df, well_names)
    for well in well_names:
        stats.setdefault(well, []).append(np.sum(df["Well"] == well) - np.sum(result["Well"] == well))
    df = result
    result = filter_cells_by_assigned_counts_percentile(df, batch_names, well_names)
    for well in well_names:
        stats.setdefault(well, []).append(np.sum(df["Well"] == well) - np.sum(result["Well"] == well))
    df = result
    result = filter_cells_by_assigned_rate(df, batch_names, well_names, 0.5)
    for well in well_names:
        stats.setdefault(well, []).append(np.sum(df["Well"] == well) - np.sum(result["Well"] == well))
    for well in well_names:
        stats.setdefault(well, []).append(np.sum(result["Well"] == well))

    stats_df = pd.DataFrame.from_dict(stats)
    stats_df.index = ["Area", "AssignedCounts", "AssignedRate", "Passing"]
    stats["filter_cells"] = stats_df
    return result
