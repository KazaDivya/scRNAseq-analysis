# ============================================================
# 05_DESeq2.R
# Differential Expression Analysis - Neurons vs Astrocytes
# ============================================================
rm(list = ls())

library(Seurat)
library(DESeq2)
library(ggplot2)
library(dplyr)



# Load annotated object
seurat_obj <- readRDS("GSE67835/seurat_annotated.rds")
cat("Loaded:", ncol(seurat_obj), "cells\n")
print(table(seurat_obj$cell_type))

# ---- Set cell type as active identity ----
Idents(seurat_obj) <- "cell_type"

# ---- Run DESeq2: Neurons vs Astrocytes ----
cat("\nRunning DESeq2: Neurons vs Astrocytes...\n")

deg_results <- FindMarkers(seurat_obj,
                           ident.1 = "Neurons",
                           ident.2 = "Astrocytes",
                           test.use = "DESeq2",
                           min.pct = 0.25,
                           logfc.threshold = 0.25)

cat("Done!\n")
cat("Total DEGs found:", nrow(deg_results), "\n")

# ---- Clean up results ----
deg_results$gene <- rownames(deg_results)
deg_results <- deg_results %>%
  arrange(p_val_adj) %>%
  mutate(direction = case_when(
    avg_log2FC > 0 & p_val_adj < 0.05 ~ "Up in Neurons",
    avg_log2FC < 0 & p_val_adj < 0.05 ~ "Up in Astrocytes",
    TRUE ~ "Not Significant"
  ))

# Print top results
cat("\nTop 10 significant DEGs:\n")
print(head(deg_results, 10))

# Count up/down
cat("\nSummary:\n")
print(table(deg_results$direction))

# Save results
write.csv(deg_results,
          "GSE67835/DEGs_Neurons_vs_Astrocytes.csv",
          row.names = FALSE)
cat("Saved: DEGs_Neurons_vs_Astrocytes.csv\n")

# ---- Volcano Plot ----
p_volcano <- ggplot(deg_results,
                    aes(x = avg_log2FC,
                        y = -log10(p_val_adj),
                        color = direction)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(
    "Up in Neurons"    = "#E41A1C",
    "Up in Astrocytes" = "#377EB8",
    "Not Significant"  = "grey70"
  )) +
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "black") +
  labs(title = "DEGs: Neurons vs Astrocytes",
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)",
       color = "") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14))

print(p_volcano)

png("figures/05_volcano_plot.png", width = 900, height = 700, res = 150)
print(p_volcano)
dev.off()
cat("Saved: figures/05_volcano_plot.png\n")

# Save object
saveRDS(seurat_obj, "GSE67835/seurat_DEG.rds")