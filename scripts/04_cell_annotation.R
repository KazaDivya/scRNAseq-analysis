# ============================================================
# 04_cell_annotation.R
# Annotate clusters with brain cell type labels
# ============================================================
rm(list = ls())

library(Seurat)
library(ggplot2)


# Load clustered object
seurat_obj <- readRDS("GSE67835/seurat_clustered.rds")
cat("Loaded:", ncol(seurat_obj), "cells\n")

# Assign cell type names to clusters
new_labels <- c(
  "0" = "Neurons",
  "1" = "Fetal_Neurons",
  "2" = "Fetal_Neurons",
  "3" = "Astrocytes",
  "4" = "OPC",
  "5" = "Oligodendrocytes",
  "6" = "Microglia"
)

seurat_obj <- RenameIdents(seurat_obj, new_labels)

# Store labels in metadata
seurat_obj$cell_type <- Idents(seurat_obj)

# Check counts per cell type
cat("Cells per type:\n")
print(table(seurat_obj$cell_type))

# Plot annotated UMAP
p <- DimPlot(seurat_obj,
             reduction = "umap",
             label = TRUE,
             label.size = 5,
             pt.size = 1.5,
             repel = TRUE) +
  ggtitle("UMAP - Brain Cell Types") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14))

print(p)

png("figures/04_UMAP_annotated.png", width = 1000, height = 700, res = 150)
print(p)
dev.off()
cat("Saved: figures/04_UMAP_annotated.png\n")

# Save annotated object
saveRDS(seurat_obj, file = "GSE67835/seurat_annotated.rds")
cat("Saved: seurat_annotated.rds\n")