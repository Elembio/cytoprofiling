library(rjson)
library(arrow)
library(Seurat)
library(cytoprofiling)

panel_json_filename <- ""
parquet_filename <- ""

panel_json <- fromJSON(file = panel_json_filename)
teton_df <- load_cytoprofiling(parquet_filename)
teton_df <- normalize_cytoprofiling(filter_cells(teton_df))
seurat_data <- cytoprofiling_to_seurat(teton_df, panel_json, use_gene_names = FALSE, split_assays = FALSE)

Seurat::DefaultAssay(seurat_data) <- "Barcoding"
seurat_data <- NormalizeData(seurat_data, normalization.method = "LogNormalize", scale.factor = 5000)
Seurat::VariableFeatures(seurat_data) <- rownames(seurat_data)
all.genes <- rownames(seurat_data)
seurat_data <- ScaleData(seurat_data, features = all.genes)
seurat_data <- RunPCA(seurat_data, features = VariableFeatures(object = seurat_data), npcs = 20, reduction.name = "barcoding_pca")
seurat_data <- RunUMAP(seurat_data, dims = 1:10, reduction = "barcoding_pca", reduction.name = "umap.barcoding", reduction.key = "barcodingUMAP_")
Seurat::DimPlot(seurat_data, reduction = "umap.barcoding")