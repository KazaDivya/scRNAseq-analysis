# ============================================================
# 03_UMAP_clustering.R
# Find neighbours, cluster cells, run UMAP
# ============================================================
rm(list = ls())

library(Seurat)
library(ggplot2)



# ---- Load from Step 02 output ----
seurat_obj <- readRDS("GSE67835/seurat_normalized.rds")
cat("Loaded:", ncol(seurat_obj), "cells\n")

# ---- Step 1: Find Neighbours ----
# Builds a graph of similar cells using the top 15 PCs
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
cat("Neighbours found\n")

# ---- Step 2: Find Clusters ----
# Resolution controls number of clusters
# Lower = fewer clusters, Higher = more clusters
# 0.5 is standard starting point for ~500 cells
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
cat("Clusters found\n")
cat("Number of clusters:", length(unique(seurat_obj$seurat_clusters)), "\n")
print(table(seurat_obj$seurat_clusters))

# ---- Step 3: Run UMAP ----
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)
cat("UMAP done\n")

# ---- Plot 1: UMAP coloured by cluster ----
p1 <- DimPlot(seurat_obj,
              reduction = "umap",
              label = TRUE,
              label.size = 5,
              pt.size = 1.5) +
  ggtitle("UMAP - Clusters") +
  theme_classic()

print(p1)

# Save
png("figures/03_UMAP_clusters.png", width = 900, height = 700, res = 150)
print(p1)
dev.off()
cat("Saved: figures/03_UMAP_clusters.png\n")

# ---- Plot 2: UMAP coloured by known marker genes ----
# These are the brain cell type markers we saw in PCA
marker_genes <- c("AQP4",    # Astrocytes
                  "MAP1B",   # Neurons
                  "PLP1",    # Oligodendrocytes
                  "CD74",    # Microglia
                  "CLDN5")   # Endothelial cells

p2 <- FeaturePlot(seurat_obj,
                  features = marker_genes,
                  ncol = 3,
                  pt.size = 1) &
  theme_classic()

print(p2)

# Save
png("figures/03_UMAP_markers.png", width = 1200, height = 800, res = 150)
print(p2)
dev.off()
cat("Saved: figures/03_UMAP_markers.png\n")

# ---- Save clustered object ----
saveRDS(seurat_obj, file = "GSE67835/seurat_clustered.rds")
cat("Saved: seurat_clustered.rds\n")