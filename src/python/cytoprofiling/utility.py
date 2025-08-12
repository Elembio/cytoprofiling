import pandas as pd
import numpy as np
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

def get_cell_profiler_modules() -> list[str]:
    """Get default list of cell profiler modules included in RawCellStats output
    
    Returns:
      List of cell profiler modules
    
    """
    return ["Texture", "Intensity", "Granularity", "RadialDistribution", "Location", "AreaShape"]

def get_imaging_targets_from_panel(panel_json : dict) -> list[str]:
    """Get list of imaging targets from panel json
    
    Args:
      panel_json : Panel as json dictionary object

    Returns:
      List of imaging target names
    """
    return list(target_data["Target"] for target_data in panel_json["ImagingTargets"])

def get_barcoding_batches_from_panel(panel_json : dict) -> list[str]: 
    """Get list of names of barcoding batches from panel json
    
    Args:
      panel_json : Panel as json dictionary object

    Returns:
      List of barcoding batch names
    """
    return list(batch_data["BatchName"] for batch_data in panel_json["BarcodingPrimerTubes"])

def get_imaging_batches_from_panel(panel_json : dict) -> list[str]:
    """Get list of names of imaging batches from panel json
    
    Args:
      panel_json : Panel as json dictionary object

    Returns:
      List of imaging batch names
    """
    return list(batch_data["BatchName"] for batch_data in panel_json["ImagingPrimerTubes"])

def get_targets_for_batch(panel_json: dict, batch_name: str) -> list[str]:
  """
  Get list of target names for a specific barcoding batch from the panel JSON.

  Args:
    panel_json: Panel configuration as a dictionary (as read by json.load)
    batch_name: Name of the barcoding batch.

  Returns:
    List of target names for the specified batch.
  """
  return [target["Target"] for target in panel_json["BarcodingTargets"] if target["BatchName"] == batch_name]

def get_target_name_for_barcode_index(
  barcode_index: int,
  panel_json: dict,
  batch_name: str
) -> str:
  """
  Returns the target name corresponding to a given barcode index for a specific batch.

  Args:
    barcode_index: The index of the barcode as read from barcodes parquet file. If 0, returns "Unassigned".
    panel_json: Panel configuration as a dictionary (as read by json.load)
    batch_name: The name of the batch to retrieve targets from.

  Returns:
    The target name corresponding to the barcode index, or "Unassigned" if index is 0.

  Raises:
    ValueError: If the barcode index is out of range for the specified batch.
  """
  if barcode_index == 0:
    return "Unassigned"
  targets = get_targets_for_batch(panel_json, batch_name)
  if barcode_index > len(targets):
    raise ValueError(f"Barcode index {barcode_index} is out of range for batch '{batch_name}' with {len(targets)} targets.")
  return targets[barcode_index - 1]
  

def cytoprofiling_to_anndata(df : pd.DataFrame, panel_json : dict = None, drop_na : bool = True) -> ad.AnnData:
  """Convert cytoprofiling dataframe to anndata object

  The following obs will be set from the dataframe:
    - WellLabel
    - Well
    - Tile
    - X
    - Y
  
  The following varm column annotations will be always be present:
    - is_nuclear
    - is_unassigned
    - cell_profiler_module
  
  If panel_json is not None, the following varm column annotations will be added:
    - batch
    - gene
    - imaging_target
    - measurement_type
    - control_type
    _ probe_concentration
  
  For all varm annotations, "NA" will be used if the annotation is unknown or not applicable. 

  Args:
    df : Input dataframe
    panel_json : Panel as json dictionary object (default: None). If string provided, will load json from named file. 
    drop_na : Drop rows and columns with any NA values (default: True)

  Returns:
    Anndata object
  """ 
  if not _has_anndata:
      raise ImportError("anndata is required for cytoprofiling_to_anndata")

  def assign_cell_profile_modules(data_columns, cell_profiler_modules):
      result = [""] * len(data_columns)
      for idx in range(len(data_columns)):
          assigned_category = "NA"
          for cell_profiler_module in cell_profiler_modules:
              if data_columns[idx].startswith(cell_profiler_module):
                  assigned_category = cell_profiler_module
                  break
          result[idx] = assigned_category
      return np.array(result)

  def assign_batches(data_columns, batches):
      result = [""] * len(data_columns)
      for idx in range(len(data_columns)):
          assigned_category = "NA"
          for batch in batches:
              if f".{batch}_" in data_columns[idx] or data_columns[idx].endswith(f".{batch}."):
                  assigned_category = batch
                  break
          result[idx] = assigned_category
      return np.array(result)

  def assign_imaging_target(data_columns, imaging_target):
      result = [""] * len(data_columns)
      for idx in range(len(data_columns)):
          assigned_category = "NA"
          for target in imaging_target:
              if f"_{target}." in data_columns[idx] or data_columns[idx].startswith(f"{target}.") or data_columns[idx].startswith(f"{target}_"):
                  assigned_category = target
                  break
          result[idx] = assigned_category
      return np.array(result)

  def assign_gene(data_columns, panel_json):
      result = [""] * len(data_columns)
      targets = [(target.get("Target", ""), target.get("BatchName", "")) for target in panel_json.get("BarcodingTargets",[])]
      target2gene = {}
      for target in targets:
          target2gene[target[0] + "." + target[1]] = target[0]
          target2gene[target[0] + "_Nuclear." + target[1]] = target[0]

      for idx in range(len(data_columns)):
          assigned_gene = "NA"
          if data_columns[idx] in target2gene:
              assigned_gene = target2gene[data_columns[idx]]
          result[idx] = assigned_gene
      return np.array(result)

  def assign_measurement_type(data_columns, panel_json):
      # get a list of protein targets from the panel
      protein_targets = set()
      for target in panel_json.get("BarcodingTargets",[]):
          if target.get("TargetType","") == "Protein":
              protein_targets.add(target.get("Target","") + "." + target.get("BatchName",""))
              protein_targets.add(target.get("Target","") + "_Nuclear." + target.get("BatchName",""))

      # get a list of transcript targets from the panel
      transcript_targets = set()
      for target in panel_json["BarcodingTargets"]:
          if target["TargetType"] == "Transcript":
              transcript_targets.add(target.get("Target","") + "." + target.get("BatchName",""))
              transcript_targets.add(target.get("Target","") + "_Nuclear." + target.get("BatchName",""))
  
      # get a list of imaging targets from the panel
      imaging_targets = get_imaging_targets_from_panel(panel_json)

      result = [""] * len(data_columns)
      for idx in range(len(data_columns)):
          assigned_category = "NA"
          if data_columns[idx] in protein_targets:
              assigned_category = "Protein"
          elif data_columns[idx] in transcript_targets:
              assigned_category = "RNA"
          elif data_columns[idx] in ["Area", "AreaUm", "NuclearArea", "NuclearAreaUm"]:
              assigned_category = "morphology"
          else:
              for target in imaging_targets:
                  if f"_{target}." in data_columns[idx] or data_columns[idx].startswith(f"{target}.") or data_columns[idx].startswith(f"{target}_"):
                      assigned_category = "morphology"
                      break
          result[idx] = assigned_category
      return np.array(result)

  def assign_control_type(data_columns, panel_json):
      # get a list of protein targets from the panel
      target2control_type = {}
      for target in panel_json.get("BarcodingTargets",[]) + panel_json.get("ImagingTargets",[]):
        target2control_type[target.get("Target","") + "." + target.get("BatchName","")] = target.get("ControlType", "")
        target2control_type[target.get("Target","") + "_Nuclear." + target.get("BatchName","")] = target.get("ControlType", "")
      
      result = [""] * len(data_columns)
      for idx in range(len(data_columns)):
          assigned_category = "NA"
          if data_columns[idx] in target2control_type:
              assigned_category = target2control_type[data_columns[idx]]
          result[idx] = assigned_category
      return np.array(result)
  
  def assign_probe_concentration(data_columns, panel_json):
      # get a list of protein targets from the panel
      target2probe_concentration = {}
      for target in panel_json.get("BarcodingTargets",[]) + panel_json.get("ImagingTargets",[]):
          target2probe_concentration[target.get("Target","") + "." + target.get("BatchName","")] = target.get("ProbeConcentration", float("nan"))
          target2probe_concentration[target.get("Target","") + "_Nuclear." + target.get("BatchName","")] = target.get("ProbeConcentration", float("nan"))
      
      result = [float("nan")] * len(data_columns)
      for idx in range(len(data_columns)):
          assigned_value = float("nan")
          if data_columns[idx] in target2probe_concentration:
              assigned_value = target2probe_concentration[data_columns[idx]]
          result[idx] = assigned_value
      return np.array(result)

  # if panel_json is a string, load json from file
  if isinstance(panel_json, str):
    import json
    with open(panel_json, "rb") as f:
      panel_json = json.load(f)

  # drop nans
  if drop_na:
    df = df.dropna(axis=1, how="all")
    df = df.dropna(axis=0, how="any")

  # Cell names
  cell_names = [f"Cell_{x}" for x in df["Cell"].tolist()]

  obs_columns = ["WellLabel", "Well", "Tile", "X", "Y"]

  data_columns = [col for col in df.columns if col not in ["Cell", ] + obs_columns]

  result = ad.AnnData(df[data_columns])
  result.obs_names = cell_names
  result.var_names = data_columns

  # Add additional information to carry it along
  for col in obs_columns:
      if col in df.columns:
          result.obs[col] = df[col].tolist()

  # nuclear
  result.var["is_nuclear"] = np.array(["_Nuclear." in col for col in data_columns], dtype=bool)
  # unassigned
  result.var["is_unassigned"] = np.array(["Unassigned." in col for col in data_columns], dtype=bool)
  # assign cell profiler module
  result.var["cell_profiler_module"] = assign_cell_profile_modules(data_columns, get_cell_profiler_modules())

  # batch
  if panel_json is not None:
    # batch
    result.var["batch"] = assign_batches(data_columns, get_barcoding_batches_from_panel(panel_json) + get_imaging_batches_from_panel(panel_json))
    # gene
    result.var["gene"] = assign_gene(data_columns, panel_json)
    # assign imaging target
    result.var["imaging_target"] = assign_imaging_target(data_columns, get_imaging_targets_from_panel(panel_json))
    # assign measurement type
    result.var["measurement_type"] = assign_measurement_type(data_columns, panel_json)
    # assign control type
    result.var["control_type"] = assign_control_type(data_columns, panel_json)
    # assign probe concentration
    result.var["probe_concentration"] = assign_probe_concentration(data_columns, panel_json)

  return result