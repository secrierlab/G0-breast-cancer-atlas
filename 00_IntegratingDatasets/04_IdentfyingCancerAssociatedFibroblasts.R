library(monocle3)
library(Seurat)
library(tidyverse)
library(infercnv)
library(tibble)

projectDir <- "~/BreastCancerG0arrest"
setwd("00_IntegratingDatasets")

integrated <- readRDS("data/integrated.rds"))
DefaultAssay(integrated) <- "RNA"
integrated <- DietSeurat(integrated, assays = "RNA")
Idents(integrated) <- integrated@meta.data$celltype

FeaturePlot(integrated %>% subset(subset = seurat_clusters %in% c(5, 8, 47)), features = c("FAP","ACTA2","COL11A1","COL1A2"), pt.size = 0.1, raster = F, cols = c("#F7F7F7", "#CC3311")) &
  NoAxes() & coord_fixed()

integrated$celltype[integrated$seurat_clusters == 5] <- "CAF"

celltype_palette <- c(
  "#8B0000",  # B cell
  "#A52A2A",  # CD4 T cell
  "#1f77b4",  # CAF
  "#9467bd",  # cycling lactocyte
  "#9e9ac8",  # cycling mammary luminal progenitor
  "#FF69B4",  # EC
  "#6baed6",  # Fibroblast
  "#2ca02c",  # Macrophage
  "#bcbddc",  # Mammary basal
  "#FF6347",  # mast cell
  "#B22222",  # NK/CD8 T cell
  "#dadaeb",  # pip mammary luminal cell
  "#aec7e8",  # pericyte
  "#DC143C",  # plasma
  "#FF4500",  # Treg
  "#d594d7",  # SAA2 mammary luminal progenitor
  "#c3a2f2",  # SCGB3A1 mammary luminal progenitor
  "#cbb7e4",  # Secretoglobin mammary luminal cell
  "#66c2a5",  # TAM
  "#b2e2a8",  # cDC
  "#FF7F7F"   # pDC
)
# Figure 1e:
pdf("figures/umap_celltype.pdf", width = 10, height = 6)
DimPlot(integrated, label = T, repel = T, raster = F, group.by = "celltype", pt.size = 0.1, cols = celltype_palette) + ggtitle("Cell types") + 
  NoAxes() + coord_fixed()
dev.off()

setwd(projectDir)
sessionInfo()
