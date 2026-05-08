# ============================================================
# 01_load_data_QC_filtering.R
# create Seurat object, QC & filter
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

# ---- Load file list ----
files <- list.files("GSE67835/raw_cells",
                    pattern = "\\.csv\\.gz$",
                    full.names = TRUE)
cat("Files found:", length(files), "\n")

# ---- Build count matrix ----
cat("Building count matrix...\n")
count_list <- lapply(files, function(f) {
  df <- read.table(gzfile(f), header = FALSE, sep = "\t",
                   col.names = c("Gene", "Count"),
                   strip.white = TRUE)
  return(df$Count)
})

gene_names <- read.table(gzfile(files[1]), header = FALSE,
                         sep = "\t", strip.white = TRUE)[, 1]

count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- gsub("_.*", "", basename(files))
cat("Count matrix:", nrow(count_matrix), "genes x", 
    ncol(count_matrix), "cells\n")

# ---- Create Seurat Object ----
seurat_obj <- CreateSeuratObject(counts = count_matrix,
                                 project = "GSE67835_Brain",
                                 min.cells = 3,
                                 min.features = 200)
cat("Seurat object created\n")
print(seurat_obj)

# ---- QC Metrics ----
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj,
                                                   pattern = "^MT-")
cat("\n QC Summary:\n")
cat("--- Genes per cell (nFeature_RNA) ---\n")
print(summary(seurat_obj$nFeature_RNA))
cat("--- Counts per cell (nCount_RNA) ---\n")
print(summary(seurat_obj$nCount_RNA))
cat("--- % Mitochondrial (percent.mt) ---\n")
print(summary(seurat_obj$percent.mt))

# ---- QC Violin Plot ----
png("figures/01_QC_violin_beforeFilter.png", 
    width = 1200, height = 600, res = 150)
VlnPlot(seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0.5)
dev.off()
cat("Saved: figures/01_QC_violin_beforeFilter.png\n")

# Filter out low-quality cells
# We will remove:
# - Cells with fewer than 500 genes    (empty droplets / dead cells)
# - Cells with more than 10,000 genes  (doublets — two cells in one)

# Before filtering
cat("Before filtering:\n")
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

# Apply filters
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 500 & 
                       nFeature_RNA < 10000)

# After filtering
cat("\nAfter filtering:\n")
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

# Full Seurat object summary
cat("Seurat Object Summary:")
print(seurat_obj)

# ---- Save filtered Seurat object ----
saveRDS(seurat_obj, file = "GSE67835/seurat_filtered.rds")
cat("Filtered Seurat object saved!\n")

