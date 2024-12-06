import pandas as pd
try:
    import anndata as ad
except ImportError:
    _has_anndata = False
else:
    _has_anndata = True

def get_wells(df : pd.DataFrame) -> list[str]:
    """Get list of well names in input data frame

    Args:
      df : Input dataframe

    Returns:
      List of wells
    """
    return list(sorted(df["Well"].unique()))

def get_cells_per_well(df : pd.DataFrame) -> pd.core.series.Series:
    """Get number of cells for each well

    Args:
      df : Input dataframe

    Returns:
      Number of cells for each well
    """
    return df.groupby(["Well",]).size()

def get_barcoding_batches(df : pd.DataFrame) -> list[str]:
    """Get list of all barcoding batch names in data frame

    Args:
      df : Input dataframe

    Returns:
      List of barcoding batches
    """
    return [val.split(".")[1] for val in df.columns if "Nuclear" not in val and "Unassigned" in val]

def get_all_targets(df : pd.DataFrame, batch_name: str) -> list[str]:
    """Get list of all targets for a batch

    Args:
      df : Input dataframe

    Returns:
      List of all targets (column headers from data frame)
    """
    result = []
    for col in df.columns:
        if col.endswith(f".{batch_name}"):
            result.append(col)
    return result

def get_default_normalization_targets(df : pd.DataFrame, batch_name : str) -> list[str]:
    """Get default list of targets to use for normalization of batch (excluding nuclear, decoy, and unassigned)

    Args:
      df : Input dataframe with all targets
      batch_name : Name of batch to query targets

    Returns:
      List of targets (column headers of input data frame)
    """
    result = []
    for col in df.columns:
        if col.endswith(f".{batch_name}") and "_Decoy" not in col and not col.startswith("NSB") and not col.startswith("Unassigned") and not col.endswith(f"_Nuclear.{batch_name}"):
            result.append(col)
    return result

def cytoprofiling_to_anndata(df : pd.DataFrame):
    """Convert cytoprofiling dataframe to anndata object

    Args:
      df : Input dataframe

    Returns:
      Anndata object
    """ 
    if not _has_anndata:
        raise ImportError("anndata is required for cytoprofiling_to_anndata")

    # drop nans
    df = df.dropna(axis=1, how="all")
    df = df.dropna(axis=0, how="any")

    # take barcoding counts as features
    feature_names = [x for x in df.columns if any([y in x for y in get_barcoding_batches(df)]) and not "Unassigned" in x and not "Nuclear" in x]

    # cell names
    cell_names = [f"Cell_{x}" for x in df["Cell"].tolist()]

    counts = df[feature_names]

    feature_names_no_batch = [x.split(".")[0] for x in feature_names]

    adata = ad.AnnData(counts)
    adata.obs_names = cell_names
    adata.var_names = feature_names_no_batch

    # add additional information to carry it along
    for metric in [col for col in df.columns if not col in feature_names]:
      adata.obs[metric] = df[metric].tolist()

    # removed duplicated targets
    adata = adata[:, ~adata.var_names.duplicated()].copy()

    return adata