# ============================================================
# 06_figures.R
# Violin plots and Heatmap of top DEGs
# ============================================================
rm(list = ls())

library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)


# Load objects
seurat_obj  <- readRDS("GSE67835/seurat_annotated.rds")
deg_results <- read.csv("GSE67835/DEGs_Neurons_vs_Astrocytes.csv")

Idents(seurat_obj) <- "cell_type"

cat("Loaded successfully\n")

# ---- Figure 1: Violin plots of top DEGs ----
# Pick top 3 up in each direction
top_neurons    <- deg_results %>%
  filter(direction == "Up in Neurons") %>%
  arrange(p_val_adj) %>%
  head(3) %>%
  pull(gene)

top_astrocytes <- deg_results %>%
  filter(direction == "Up in Astrocytes") %>%
  arrange(p_val_adj) %>%
  head(3) %>%
  pull(gene)

top_genes <- c(top_neurons, top_astrocytes)
cat("Top genes for violin plot:", top_genes, "\n")

# Plot only Neurons and Astrocytes for clarity
subset_obj <- subset(seurat_obj,
                     idents = c("Neurons", "Astrocytes"))

p_violin <- VlnPlot(subset_obj,
                    features = top_genes,
                    ncol = 3,
                    pt.size = 0.5,
                    cols = c("Neurons"    = "#E41A1C",
                             "Astrocytes" = "#4DAF4A")) &
  theme(plot.title   = element_text(size = 11, face = "bold"),
        axis.title.x = element_blank())

print(p_violin)

png("figures/06_violin_top_DEGs.png",
    width = 1400, height = 900, res = 150)
print(p_violin)
dev.off()
cat("Saved: figures/06_violin_top_DEGs.png\n")

# ---- Figure 2: Heatmap of top 40 DEGs ----
# Get top 20 from each direction
top20_neurons    <- deg_results %>%
  filter(direction == "Up in Neurons") %>%
  arrange(p_val_adj) %>%
  head(20) %>%
  pull(gene)

top20_astrocytes <- deg_results %>%
  filter(direction == "Up in Astrocytes") %>%
  arrange(p_val_adj) %>%
  head(20) %>%
  pull(gene)

heatmap_genes <- c(top20_neurons, top20_astrocytes)

# Use Seurat's DoHeatmap on Neurons + Astrocytes only
p_heatmap <- DoHeatmap(subset_obj,
                       features = heatmap_genes,
                       group.colors = c("Neurons"    = "#E41A1C",
                                        "Astrocytes" = "#4DAF4A"),
                       label = TRUE,
                       size = 4) +
  scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C")) +
  labs(title = "Top 20 DEGs: Neurons vs Astrocytes") +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 13))

print(p_heatmap)

png("figures/06_heatmap_top_DEGs.png",
    width = 1200, height = 1400, res = 150)
print(p_heatmap)
dev.off()
cat("Saved: figures/06_heatmap_top_DEGs.png\n")

cat("\nAll figures saved!\n")