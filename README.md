# Cytoprofiling
This repository contains methods and examples for processing data from the Element Biosciences Teton cytoprofiling assay. Both python and R packages are provided. In addition, Cellpose segmentation models trained on Cell Paint images from the Element Biosciences Teton cytoprofiling assay are included.

## Installation
After downloading the repository and changing the current directory to the downloaded location, package installation is performed in the typical fashion for each language. 

python
```
cd <path to downloaded location of repository>
pip3 install src/python

# if desired, can specify installation of additional dependences for integration with scanpy/anndata with the "scanpy" extra
pip3 install src/python[scanpy]
```

R
```
setwd("<path to downloaded location of repository>")
library("devtools")
devtools::install("src/R/cytoprofiling")
```
## Methods
Analogous methods for data processing with matching names are provided in the packages for each languages. The typical pre-processing workflow is to load the RawCellStats file into a dataframe, filter cells, and then normalize. 

### Loading

A key output of the Teton cytoprofiling assay is the RawCellStats.parquet file. This is a data table in parquet format with each row giving information about a single cell and each column representing a feature profiled across all cells. For each barcoded target "{target}" profiled in a batch "{batch}, the parquet table will contain a column "{target}.{batch}" giving the counts of that target in each cell. Each such target will also have an additional column, "{target}_Nuclear.{batch}" giving the subset of those counts that are nuclear. Examples of loading the parquet table are provided for each language below

python
```
import pandas as pd
df = pd.read_parquet("RawCellStats.parquet")
```

R
```
library(cytoprofiling)
df <- load_cytoprofiling("RawCellStats.parquet")
```

### Filtering 

#### filter_cells_by_area_percentile
Filter cells by area to remove unusually large or small objects that may represent erroneous segmentation events. 

#### filter_cells_by_assigned_counts_percentile
Filter cells by assigned counts in each batch. Cells with a very small number of assigned counts per batch may represent a localized issue with barcoding performance for that cell.  

#### filter_cells_by_assigned_rate
Filter cells by assigned rate in each batch. Cells with a low assignment rate in a batch may represent a localized issue with barcoding performance for that cell.  

#### filter_cells
Helper method to perform default filtering across all wells and batches. Performs filtering by
* Area percentile
* Assigned counts percentile
* Assigned rate
with default values. 

### Normalizing 
#### normalize_cells_by_aggregated_counts
To normalize across cells, the counts for each cells are normalized by the counts in that cell. By default, the counts will be normalized per batch by the total assigned counts in that batch; however, it is also possible to specify a custom set of targets for normalization if some invariant set of targets is known for a particular condition. 

#### normalize_wells_by_median_of_ratios
This method normalizes counts between wells to ensure that most targets show no fold change across wells

#### normalize_cytoprofiling
Helper method to perform default normalization. Normalizes
* Cells by total assigned counts per batch
* Wells by median of ratios
with default values. 

### Exporting
After pre-processing of the data, it is possible to export the data for use with popular 3rd party single cell analysis platforms. 

#### cytoprofiling_to_anndata
Export cytoprofiling python pandas dataframe to anndata object

#### cytoprofiling_to_seurat
Export cytoprofling R dataframe to Seurat object

## Examples

A simple example to execute the default pre-processing workflow is shown below

python
```
import pandas as pd
from cytoprofiling import normalize_cytoprofiling, filter_cells

df = pd.read_parquet("RawCellStats.parquet")
norm_df = normalize_cytoprofiling(filter_cells(df))
```

R
```
library(cytoprofiling)
df <- load_cytoprofiling("RawCellStats.parquet")
norm_df <- normalize_cytoprofiling(filter_cells(df))
```

Examples of downstream analysis workflows following pre-processing are supplied in the "examples" directory for each language. 