###################################################################################################
# LIVER TOXICITY ANALYSIS: Mini-panel run

################################################## Libraries ##################################################
# Load and install all required R packages. 
pkgs <- c("Seurat", "cytoprofiling", "dplyr", "ComplexHeatmap", "circlize", "ggplot2")

for (p in pkgs) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
library("BiocManager")

if (!require("glmGamPoi", quietly = TRUE)) {
  BiocManager::install("glmGamPoi")
}
library("glmGamPoi")

################################################## Color Palettes ##################################################
# Define custom color palettes for plotting each drug treatment. 
# These colors will be used in UMAPs and other figures for consistency.
drug_colors <- c(
  "DMSO"          = "#878D96",   # Vehicle control
  "Media"         = "#C1C7CD",   # No treatment control
  "Acetaminophen" = "#BC2e62",   # Classic hepatotoxin
  "Amiodarone"    = "#9c9aff",   # Anti-arrhythmic drug
  "Tunicamycin"   = "#0a7357",   # ER stress inducer
  "Diclofenac"    = "#a3edf7",   # NSAID
  "Doxorubicin"   = "#75b1f0"    # Positive control for apoptosis
)

################################################## Functions ##################################################
##### cytoprofiling_to_seurat: converts RawCellStats into a Seurat object
# Extracts assay targets per barcoding batch, unions them, and builds a cells×features matrix.
# Adds well metadata and sets identities to WellLabel for convenient grouping/plotting.
# min.features/min.cells gate out sparse features/cells.
cytoprofiling_to_seurat <- function(df, min.features = 25, min.cells = 20) {
  batches <- get_barcoding_batches(df)
  all_targets <- unique(unlist(lapply(batches, function(b) get_default_normalization_targets(df, b))))
  obj <- Seurat::CreateSeuratObject(t(df[, all_targets]),
                                    min.features = min.features, min.cells = min.cells)
  obj <- Seurat::AddMetaData(obj, as.factor(df[["Well"]]),      col.name = "Well")
  obj <- Seurat::AddMetaData(obj, as.factor(df[["WellLabel"]]), col.name = "WellLabel")
  Seurat::Idents(obj) <- obj$WellLabel
  obj
}

##### label_to_treatment: collapses triplicate wells into a single treatment condition
# Converts replicate-labeled wells (e.g., "ACE_10mM_Rep1") to a unified treatment key (e.g., "ACE-10").
# Useful for averaging/aggregation and simpler plotting facets/legends.

label_to_treatment <- function(x) {
  dplyr::case_when(
    x %in% c("AMO_30uM_Rep1","AMO_30uM_Rep2","AMO_30uM_Rep3") ~ "AMO-30",
    x %in% c("DOX_15uM_Rep1","DOX_15uM_Rep2","DOX_15uM_Rep3") ~ "DOX-10",
    x %in% c("ACE_10mM_Rep1","ACE_10mM_Rep2","ACE_10mM_Rep3") ~ "ACE-10",
    x %in% c("TUN_2.4uM_Rep1","TUN_2.4uM_Rep2","TUN_2.4uM_Rep3") ~ "TUN-2",
    x %in% c("DIC_200uM_Rep1","DIC_200uM_Rep2","DIC_200uM_Rep3") ~ "DIC-200",
    x %in% c("MEDIA_Rep1","MEDIA_Rep2","MEDIA_Rep3") ~ "MEDIA",
    x %in% c("DMSO_Rep1","DMSO_Rep2","DMSO_Rep3") ~ "DMSO",
    TRUE ~ NA_character_
  )
}

##### treatment_to_drug_group: maps treatment codes to drug labels
treatment_to_drug_group <- function(x) {
  dplyr::case_when(
    grepl("^ACE", x) ~ "Acetaminophen",
    grepl("^AMO", x) ~ "Amiodarone",
    grepl("^TUN", x) ~ "Tunicamycin",
    grepl("^DIC", x) ~ "Diclofenac",
    grepl("^DOX", x) ~ "Doxorubicin",
    x == "DMSO"      ~ "DMSO",
    x == "MEDIA"     ~ "Media",
    TRUE             ~ "Other"
  )
}

################################################## Load & Subset Wells ##################################################
##### Load RawCellStats file (per-cell matrix + metadata from CytoProfiling)
# Replace your filename here
mini_file <- "RawCellStats.parquet"
mini_raw  <- load_cytoprofiling(mini_file)

##### Define the wells to include downstream
# Keeps controls, four drugs, and doxorubicin; MEDIA_Rep3 was excluded due to QC.

wells_to_analyze <- c(
  "DMSO_Rep1","DMSO_Rep2","DMSO_Rep3",
  "ACE_10mM_Rep1","ACE_10mM_Rep2","ACE_10mM_Rep3",
  "AMO_30uM_Rep1","AMO_30uM_Rep2","AMO_30uM_Rep3",
  "DOX_15uM_Rep1","DOX_15uM_Rep2","DOX_15uM_Rep3",
  "TUN_2.4uM_Rep1","TUN_2.4uM_Rep2","TUN_2.4uM_Rep3",
  "DIC_200uM_Rep1","DIC_200uM_Rep2","DIC_200uM_Rep3",
  "MEDIA_Rep1","MEDIA_Rep2"
)

mini_raw <- subset(mini_raw, subset = WellLabel %in% wells_to_analyze)

################################################## Store morphology as metadata ##################################################
##### Build a morphology metadata frame (per cell) to attach to Seurat
# Remove protein features (suffix .B01–.B08) so only morphology + RNA remain here.
# Add a "treatment" column that collapses triplicates (e.g., "ACE-10.1" -> "ACE-10").
morph_meta_raw <- mini_raw %>%
  dplyr::select(!matches("\\.B0[1-8]$")) %>% 
  dplyr::mutate(treatment = sub("\\.[0-9]+$", "", WellLabel))

#### (Optional) Excludes doxorubicin from morphology data
morph_meta_filtered <- morph_meta_raw %>%
  dplyr::filter(!is.na(treatment), treatment != "DOX-10")

##### Saves morphology matrices for later analysis, if desired
#saveRDS(morph_meta_raw, "morph_meta_raw.rds")
#saveRDS(morph_meta_filtered, "morph_meta_filtered.rds")

################################################## SCTransform Workflow ##################################################
##### Cell-level QC/filtering using cytoprofiling defaults
# Removes segmentation outliers, low-assignment cells, and extreme total count outliers.
mini_filt <- filter_cells(mini_raw)

##### Create a Seurat object from filtered cells
sc <- cytoprofiling_to_seurat(mini_filt)

##### Normalize and variance-stabilize counts with SCTransform
# Replaces log-normalization; handles depth correction and residual variance.
sc <- SCTransform(sc, assay = "RNA", verbose = FALSE)

##### Keeps treatments in a particular order for downstream visualization (Negative controls left; treatments center; positive control right)
sc$treatment  <- label_to_treatment(sc$WellLabel)
sc$drug_group <- factor(
  treatment_to_drug_group(sc$treatment),
  levels = c("DMSO", "Media", "Acetaminophen", "Amiodarone", "Diclofenac", "Tunicamycin", "Doxorubicin")
)

##### Attaches morphology columns to Seurat metadata
sc <- AddMetaData(sc, metadata = morph_meta_raw[rownames(sc@meta.data), , drop = FALSE])

##### Dimensionality reduction: PCA followed by UMAP
# Use a fixed seed for reproducible embeddings.
sc <- RunPCA(sc, npcs = 30, features = rownames(sc))
sc <- RunUMAP(sc, dims = 1:30, seed.use = 1234)

##### Figure 1: UMAP by drug treatment
DimPlot(sc, group.by = "drug_group", cols = drug_colors) + ggtitle(NULL)

################################################## Target Groups & Heatmaps ##################################################
##### Define pathway/grouping blocks for heatmap row structure
# Lists of RNA/protein targets per biological theme
gene_groups <- list(
  "Apoptosis"                 = c("CASP3-CSA.B08","BAX-CSA.B08","phosS70-Bcl2.B01"),
  "Cell Cycle"                = c("RRM2.B03","DSCC1.B04","CDK1.B03","KIF20B.B02","CENPE.B04","DLGAP5.B03","TTK.B06","Ki67.B01"),
  "Mitochondrial Dysfunction" = c("CAPN2.B05","SIRT3-CM.B08","LY96.B05","FASTK-CM.B08","phosY705-STAT3.B01"),
  "Chemokines"                = c("CXCL8.B05","CXCL10.B06"),
  "Cytokines"                 = c("IL18.B06","IL-6.B01","IL8-InnImm.B08","TNF.B01","TGFB1.B06")
)
features_for_avg <- unique(unlist(gene_groups))

##### Compute average expression per drug group for selected features
# Aggregates SCTransformed expression by 'drug_group' for both RNA and protein features.
avg <- AverageExpression(
  sc,
  features = features_for_avg,
  group.by = "drug_group",
  assays = DefaultAssay(sc)
)[[DefaultAssay(sc)]]

##### Build a combined negative-control baseline and compute per-feature Z-scores
# Combine DMSO + Media as "Neg Control".
# Z-score features across groups, then center each feature on the Neg Control mean.
stopifnot(all(c("DMSO","Media") %in% colnames(avg)))
avg <- cbind(avg, "Neg Control" = rowMeans(avg[, c("DMSO","Media"), drop = FALSE], na.rm = TRUE))
z   <- t(scale(t(avg)))
zc  <- sweep(z, 1, z[, "Neg Control"], "-")
zplot <- zc[, setdiff(colnames(zc), c("Neg Control","DMSO","Media")), drop = FALSE]


##### Order rows by predefined groups and create split factors for block labels
# Keep only features present, preserve group order, and build a factor for group splits.
gene_groups <- lapply(gene_groups, function(g) intersect(g, rownames(zplot)))
gene_groups <- gene_groups[lengths(gene_groups) > 0]
row_order   <- unlist(gene_groups, use.names = FALSE)
zplot       <- zplot[row_order, , drop = FALSE]
row_split   <- factor(rep(names(gene_groups), lengths(gene_groups)), levels = names(gene_groups))

##### Create a symmetric color scale around 0 and clip extremes
# Uses the 98th percentile to avoid outliers dominating the color mapping.
lim <- as.numeric(quantile(abs(zc), 0.98, na.rm = TRUE))
col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("#3074E3","white","#BC2E62"))

################################################## Split RNA vs Protein Heatmaps ##################################################
##### Separate protein from RNA features
# Protein features carry .B01–.B08 suffixes; everything else is RNA
all_features     <- rownames(zplot)
protein_features <- grep("\\.(B01|B08)$", all_features, value = TRUE)
rna_features     <- setdiff(all_features, protein_features)
zplot_protein <- zplot[protein_features, , drop = FALSE]
zplot_rna     <- zplot[rna_features,     , drop = FALSE]

# Build split factors for each heatmap after subsetting
if (is.null(names(row_split))) names(row_split) <- rownames(zplot)
row_split_protein <- droplevels(row_split[match(rownames(zplot_protein), names(row_split))])
row_split_rna     <- droplevels(row_split[match(rownames(zplot_rna),     names(row_split))])

##### Clean row labels for display
# Remove trailing .Bxx suffixes and make names unique to avoid duplicate row labels.
clean_suffix <- function(x) make.unique(sub("\\.B\\d{2}$", "", x))
rownames(zplot_protein) <- clean_suffix(rownames(zplot_protein))
rownames(zplot_rna)     <- clean_suffix(rownames(zplot_rna))

##### Build heatmaps (protein and RNA)
# No column clustering: preserves treatment order set earlier.
# Row splits show pathway blocks; dendrograms within blocks can reveal feature substructure
ht_prot <- ComplexHeatmap::Heatmap(
  zplot_protein, name = "Z", col = col_fun,
  cluster_rows = TRUE, cluster_columns = FALSE,
  row_split = row_split_protein, show_row_dend = FALSE,
  row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 11),
  row_names_side = "left", row_names_gp = grid::gpar(fontsize = 10),
  column_title = NULL,
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
)
ht_rna <- ComplexHeatmap::Heatmap(
  zplot_rna, name = "Z", col = col_fun,
  cluster_rows = TRUE, cluster_columns = FALSE,
  row_split = row_split_rna, show_row_dend = FALSE,
  row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 11),
  row_names_side = "left", row_names_gp = grid::gpar(fontsize = 10),
  column_title = NULL,
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
)


##### Draw heatmaps to the active device
# Render protein first, then RNA, so they appear one after another.
draw(ht_prot)
draw(ht_rna)