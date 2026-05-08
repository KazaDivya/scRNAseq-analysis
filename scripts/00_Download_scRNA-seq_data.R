# ============================================================
# scRNA-seq Analysis: GSE67835 Human Brain Cells
# Darmanis et al., 2015 - PNAS
# ============================================================

# ---- Step 01: Load Libraries ----
library(Seurat)
library(GEOquery)
library(dplyr)
library(ggplot2)

# ---- Step 01: Download Data from GEO ----
getGEOSuppFiles("GSE67835")
untar("GSE67835/GSE67835_RAW.tar", exdir = "GSE67835/raw_cells")

# Confirm download
files <- list.files("GSE67835/raw_cells",
                    pattern = "\\.csv\\.gz$",
                    full.names = TRUE)
cat("Total cells downloaded:", length(files), "\n")

# ---- Step 01: Download Metadata ----
gse <- getGEO("GSE67835", GSEMatrix = TRUE, AnnotGPL = FALSE)
# Combine both platform metadata
metadata1 <- pData(gse[[1]])  # GPL15520 - 138 cells
metadata2 <- pData(gse[[2]])  # GPL18573 - remaining cells

cat("Platform 1 rows:", nrow(metadata1), "\n")
cat("Platform 2 rows:", nrow(metadata2), "\n")

# Combine them
metadata <- rbind(metadata1, metadata2)
cat("Total metadata rows:", nrow(metadata), "\n")

# Save combined metadata
write.csv(metadata, "GSE67835/metadata_raw.csv", row.names = TRUE)

# ---- Step 01: Peek at one cell CSV ----
one_cell <- read.table(gzfile(files[1]), 
                       header = FALSE, 
                       sep = "\t",
                       col.names = c("Gene", "Count"),
                       strip.white = TRUE)  # removes extra spaces

cat("Genes measured:", nrow(one_cell), "\n")
head(one_cell, 10)



