# scRNA-seq Analysis: Human Brain Cell Types
### Differential Expression Analysis of GSE67835 (Darmanis et al., 2015)

---

## Overview

This project performs a complete single-cell RNA sequencing (scRNA-seq) analysis pipeline on publicly available human brain cell data. Starting from raw count matrices downloaded from GEO, the pipeline covers quality control, normalization, dimensionality reduction, clustering, cell type annotation, differential expression analysis, and figure generation.

**Key comparison:** Neurons vs Astrocytes differential gene expression

---

## Dataset

| Field | Details |
|-------|---------|
| **GEO Accession** | [GSE67835](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835) |
| **Paper** | Darmanis et al., *"A survey of human brain transcriptome diversity at the single cell level"*, PNAS 2015 |
| **Technology** | SMART-seq (full-length scRNA-seq) |
| **Organism** | Homo sapiens |
| **Tissue** | Human cerebral cortex |
| **Cells** | 466 single cells |
| **Platforms** | Illumina MiSeq (GPL15520), Illumina NextSeq 500 (GPL18573) |

---

## Cell Types Identified

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0 | Neurons | SOX11, MAP1B, DCX, STMN2 |
| 1, 2 | Fetal Neurons | MAP1B, DCX, BCL11A |
| 3 | Astrocytes | AQP4, GJA1, AGT, AGXT2L1 |
| 4 | OPC | PDGFRA |
| 5 | Oligodendrocytes | PLP1, MBP, ERMN, CNP |
| 6 | Microglia | CD74, CCL4, SPP1 |

---

## Analysis Pipeline

```
00_download_data.R        → Download GSE67835 from GEO (466 CSV files + metadata)
01_load_data_QC_filtering.R  → Build count matrix, QC metrics, filter low-quality cells
02_normalization.R        → LogNormalize, find variable features, scale, PCA
03_UMAP_clustering.R      → FindNeighbors, FindClusters (res=0.5), RunUMAP
04_cell_annotation.R      → Annotate clusters with known brain cell type markers
05_DESeq2.R               → Differential expression: Neurons vs Astrocytes
06_figures.R              → Violin plots and heatmap of top DEGs
```

---

## Key Results

### QC Summary (after filtering)
| Metric | Min | Median | Max |
|--------|-----|--------|-----|
| Genes per cell | 696 | 4,024 | 9,988 |
| Counts per cell | 431,090 | 2,224,995 | 6,298,797 |
| Cells after QC | — | **465** | — |

> Note: High counts are expected for SMART-seq full-length sequencing technology.
> Mitochondrial % = 0 across all cells (pre-filtered by authors before GEO deposit).

### Differential Expression: Neurons vs Astrocytes

| Direction | Genes |
|-----------|-------|
| Up in Neurons | SOX11, SOX4, TUBA1A, DCX, STMN2 |
| Up in Astrocytes | AGXT2L1, GJA1, AGT, VIP, TMEM144 |

Top findings are consistent with known brain cell biology:
- **SOX11** — master neuronal transcription factor, exclusively in Neurons
- **GJA1** (Connexin-43) — hallmark astrocyte gap junction protein
- **AGT** (Angiotensinogen) — textbook astrocyte secreted factor
- **DCX** (Doublecortin) — neuronal migration marker

---

## References

Darmanis S, Sloan SA, Zhang Y, Enge M, Caneda C, Shuer LM, Hayden Gephart MG, Barres BA, Quake SR.
*A survey of human brain transcriptome diversity at the single cell level.*
Proc Natl Acad Sci U S A. 2015 Jun 9;112(23):7285-90.
PMID: [26060301](https://pubmed.ncbi.nlm.nih.gov/26060301/)

---

## Author

**DivyaKaza** | dkdivya1138@gmail.com
