###################################################################################################
##### Data Spotlight: T Cell Activation

################################################## Libraries ##################################################
# Load all required R packages.
suppressPackageStartupMessages({
  library(Seurat)           # For single-cell data integration, normalization, clustering, and visualization
  library(cytoprofiling)    # For loading and preprocessing Teton CytoProfiling data
  library(dplyr)            # For general data wrangling
  library(ggplot2)          # For plotting
})

# Use Arial on Windows
if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

################################################## Globals ##################################################
# Define custom color palettes for plotting each treatment.
# These colors will be used in UMAPs and other figures for consistency.
drug_colors <- c(
  "Media + DMSO"     = "#878D96",
  "Media"            = "#C1C7CD",
  "CD3/CD28 30min"   = "#0a7537",
  "CD3/CD28 6h"      = "#9c9aff",
  "CD3/CD28 24h"     = "#0f41a1",
  "Ionomycin + PMA"  = "#BC2e62"
)

# Common ordering of treatments for legends and facets
treat_levels_keep <- c("Media + DMSO", "CD3/CD28 30min", "CD3/CD28 6h", "CD3/CD28 24h", "Ionomycin + PMA")

# Map treatments to numeric time for kinetics
time_map <- c(
  "Media + DMSO"    = 0,
  "CD3/CD28 30min"  = 0.5,
  "CD3/CD28 6h"     = 6,
  "CD3/CD28 24h"    = 24
)

# labels for the time axis
x_labels <- c("0h", "30m", "6h", "24h")

# Plot defaults
base_font_family <- "Arial"
base_font_size   <- 10
legend_point_sz  <- 2

# Reusable ggplot theme with Arial and uniform sizing
theme_aviti <- function() {
  theme_classic(base_size = base_font_size, base_family = base_font_family) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      axis.text        = element_text(family = base_font_family, size = base_font_size),
      axis.title       = element_text(family = base_font_family, size = base_font_size),
      plot.title       = element_text(family = base_font_family, size = base_font_size, face = "bold"),
      legend.position  = "top"
    )
}

################################################## Functions ##################################################
##### cytoprofiling_to_seurat: converts RawCellStats into a Seurat object
# Extracts assay targets per barcoding batch, unions them, and builds a cells×features matrix.
# Adds well metadata and sets identities to WellLabel for convenient grouping/plotting.
# min.features/min.cells gate out sparse features/cells.
cytoprofiling_to_seurat <- function(df, min.features = 10, min.cells = 20) {
  batches <- get_barcoding_batches(df)
  all_targets <- unique(unlist(lapply(batches, function(b) get_default_normalization_targets(df, b))))
  obj <- Seurat::CreateSeuratObject(
    counts = t(df[, all_targets]),
    min.features = min.features,
    min.cells = min.cells
  )
  obj <- Seurat::AddMetaData(obj, as.factor(df[["Well"]]),      col.name = "Well")
  obj <- Seurat::AddMetaData(obj, as.factor(df[["WellLabel"]]), col.name = "WellLabel")
  Seurat::Idents(obj) <- obj$WellLabel
  obj
}

##### label_to_treatment: collapses duplicate wells into a single treatment condition
label_to_treatment <- function(x) {
  dplyr::case_when(
    x %in% "Ionomycin_PMA"                ~ "Ionomycin + PMA",
    x %in% "Media_DMSO"                   ~ "Media + DMSO",
    x %in% "Media"                        ~ "Media",
    x %in% "Anti-CD3_Anti-CD28_30min"     ~ "CD3/CD28 30min",
    x %in% "Anti-CD3_Anti-CD28_6hr"       ~ "CD3/CD28 6h",
    x %in% "Anti-CD3_Anti-CD28_24hr"      ~ "CD3/CD28 24h",
    TRUE                                  ~ NA_character_
  )
}

################################################## Load & Subset Wells ##################################################
# Load RawCellStats file (per-cell matrix + metadata from CytoProfiling)
#### IMPORTANT: replace the file name below with the location of your parquet file
parquet_file  <- "CPKOL-0053_RawCellStats.parquet"
cyto   <- load_cytoprofiling(parquet_file)
cyto  <- filter_cells(cyto)
t.cell.activation <- cytoprofiling_to_seurat(cyto)  # ~48k cells

# Drop feature columns (e.g., .B01-.B08) from metadata and attach remaining metadata
cyto_metadata <- cyto %>% select(-matches("\\.B0[1-8]$"))
t.cell.activation    <- AddMetaData(t.cell.activation, metadata = cyto_metadata)

# Basic cell filtering and downsampling
t.cell.activation          <- subset(t.cell.activation, subset = nFeature_RNA > 50 & nCount_RNA > 100)
t.cell.activation <- subset(t.cell.activation, downsample = 3000)

# Add treatment label
t.cell.activation$treatment <- label_to_treatment(t.cell.activation$WellLabel)

# Normalize and embed
t.cell.activation <- SCTransform(t.cell.activation, assay = "RNA", verbose = FALSE,
                         new.assay.name = "SCT", vst.flavor = "v2", min_cells = 1)
t.cell.activation <- RunPCA(t.cell.activation, npcs = 30, features = rownames(t.cell.activation))
t.cell.activation <- RunUMAP(t.cell.activation, dims = 1:30, seed.use = 1234)

####################################################### Figure 1A: Merged UMAP #######################################################
# Use a consistent legend order and palette
DimPlot(t.cell.activation, reduction = "umap", group.by = "treatment") +
  ggtitle(NULL) +
  scale_color_manual(values = drug_colors, breaks = treat_levels_keep, labels = treat_levels_keep) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = legend_point_sz, shape = 16))) +
  theme_aviti() +
  theme(strip.text = element_blank())

####################################################### Figure 1B: Split UMAP #######################################################
t.sub <- subset(t.cell.activation, subset = treatment %in% treat_levels_keep & !is.na(treatment))
t.sub$treatment <- factor(t.sub$treatment, levels = treat_levels_keep)

DimPlot(t.sub, reduction = "umap", split.by = "treatment", group.by = "treatment") +
  ggtitle(NULL) +
  scale_color_manual(values = drug_colors, breaks = treat_levels_keep, labels = treat_levels_keep) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = legend_point_sz, shape = 16))) +
  theme_aviti() +
  theme(strip.text = element_blank())

####################################################### Fig 2A: Marker Gene Identification and Heatmap #######################################################
# Drop unwanted treatments and find top markers per cluster
drop_treatments <- c("Media", "Ionomycin + PMA")  # <-- edit these if needed
t.sub.mark      <- subset(t.cell.activation, subset = !treatment %in% drop_treatments & !is.na(treatment))

table(t.sub.mark$treatment)

markers <- FindAllMarkers(t.sub.mark, only.pos = TRUE)
top10 <- markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.1) %>%
  slice_head(n = 10) %>%
  ungroup()

DoHeatmap(t.sub.mark, features = top10$gene, assay = "SCT", slot = "scale.data") + NoLegend()

############################################# Early phospho kinetics & Figure 2B ########################################
# Define features once
marks <- c("phosS240-S244-RPS6.B01", "CD69.B04", "NR4A1.B07", "IL2RA.B04", "CD28.B04")

# Build long-form data, filter to time-mapped treatments, and compute summary stats
dat <- FetchData(t.cell.activation, vars = c("treatment", marks), assay = "protein", layer = "data") %>%
  tidyr::pivot_longer(cols = all_of(marks), names_to = "feature", values_to = "value") %>%
  filter(treatment %in% names(time_map)) %>%
  mutate(time_hr = unname(time_map[treatment]))

stats <- dat %>%
  group_by(time_hr, treatment, feature) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    se   = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop_last"
  ) %>%
  ungroup()

# Normalize to 0h within each feature (fold-change) and keep SE on same scale
stats_norm <- stats %>%
  group_by(feature) %>%
  mutate(
    mean0     = mean[time_hr == 0],
    mean_norm = mean / mean0,
    se_norm   = se   / mean0
  ) %>%
  ungroup()

# Figure 2B: timeline of normalized expression
ggplot(stats_norm, aes(x = time_hr, y = mean_norm, color = feature, group = feature)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_norm - 1.96 * se_norm,
                    ymax = mean_norm + 1.96 * se_norm),
                width = 0.3, linewidth = 0.6) +
  scale_x_continuous(
    breaks = c(0, 0.5, 6, 24),
    labels = x_labels,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    x = "Time after stimulation",
    y = "Fold-change vs 0h",
    title = "Normalized kinetics (0h = 1)",
    color = "Feature"
  ) +
  theme_aviti()
