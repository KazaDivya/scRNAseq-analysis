# ============================================================
# 02_normalization.R
# Normalize, find variable features, scale, PCA
# ============================================================

library(Seurat)
library(ggplot2)



# ---- Load filtered object 
seurat_obj <- readRDS("GSE67835/seurat_filtered.rds")
cat("Loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes\n")

# ---- Normalize ----
# Divides each cell's counts by total counts, multiplies by 10,000, log-transforms
seurat_obj <- NormalizeData(seurat_obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
cat("Normalization done\n")

# ---- Find Variable Features ----
# Finds the 2000 most variable genes across cells (most informative for clustering)
seurat_obj <- FindVariableFeatures(seurat_obj,
                                   selection.method = "vst",
                                   nfeatures = 2000)

# Plot top variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
cat("Top 10 most variable genes:\n")
print(top10)

plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# Save variable features plot
png("figures/02_variable_features.png", width = 1000, height = 700, res = 150)
print(plot2)
dev.off()
cat("Saved: figures/02_variable_features.png\n")


# ---- Scale Data ----
# Zero-centers and scales each gene — required before PCA
seurat_obj <- ScaleData(seurat_obj)
cat("Scaling done\n")

# ---- PCA ----
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
cat("PCA done\n")

# Visualize top genes driving each PC
print(seurat_obj[["pca"]], dims = 1:3, nfeatures = 5)

# Elbow plot — helps decide how many PCs to use for UMAP
# Save elbow plot
png("figures/02_elbow_plot.png", width = 800, height = 600, res = 150)
ElbowPlot(seurat_obj, ndims = 30)
dev.off()
cat("Saved: figures/02_elbow_plot.png\n")


# ---- Save ----
saveRDS(seurat_obj, file = "GSE67835/seurat_normalized.rds")
cat("Saved normalized object!")