library(Seurat)
library(tidyverse)
library(pals)
library(ggplot2)
projectDir <- "~/BreastCancerG0arrest"
setwd("00_IntegratingDatasets")

integrated <- readRDS("data/SCTintegrated_leiden.rds")

# Annotation----
library(HGNChelper)

## Import functions & cell type marker database----
source("3a1.gene_sets_prepare.R")
source("3a2.sctype_score_.R")

dataBase <- readRDS("database/db.rds")
geneSetList <-  gene_sets_prepare(dataBase, "breast")

# if SCTransform-normalised, use "SCT" assay:
estimateMax <- sctype_score(scRNAseqData = integrated[["SCT"]]@scale.data, scaled = T, gs = geneSetList$gs_positive, gs2 = geneSetList$gs_negative)

clusterResults <- do.call("rbind", lapply(unique(integrated@meta.data$seurat_clusters), function(cl) {
  estimateMax.cl <- sort(rowSums(estimateMax[, rownames(integrated@meta.data[integrated@meta.data$seurat_clusters == cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(estimateMax.cl), scores = estimateMax.cl, ncells = sum(integrated@meta.data$seurat_clusters == cl)), 10)}))

celltypeScores <-  clusterResults %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

celltypeScores$type[as.numeric(as.character(celltypeScores$scores)) < celltypeScores$ncells / 4] = "Unknown"
print(celltypeScores[, 1:4], n = length(celltypeScores$cluster))

integrated$celltype = ""

for (j in unique(celltypeScores$cluster)) {
  cl_type = celltypeScores[celltypeScores$cluster == j, ]
  integrated$celltype[integrated$seurat_clusters == j] = as.character(cl_type$type[1])
}

Idents(integrated) <- paste0(integrated$seurat_clusters, "-", integrated$celltype)

## Rename celltypes----
integrated <- RenameIdents(integrated,
                           "1-Alveolar" = "SAA2+ mammary luminal progenitor",
                           "2-Alveolar" = "SAA2+ mammary luminal progenitor",
                           "3-Hormone sensing" = "PIP+ mammary luminal cell",
                           "4-Macrophages" = "Macrophage",
                           "5-Fibroblasts" = "Fibroblast", 
                           "6-Hormone sensing" = "PIP+ mammary luminal cell", 
                           "7-Hormone sensing" = "PIP+ mammary luminal cell", 
                           "8-Fibroblasts" = "Fibroblast", 
                           "9-Hormone sensing" = "PIP+ mammary luminal cell",
                           "10-Gamma delta T cells" = "Cycling lactocyte",
                           "11-Hormone sensing" = "PIP+ mammary luminal cell", 
                           "12-Hormone sensing" = "PIP+ mammary luminal cell", 
                           "13-Hormone sensing" = "Secretoglobin mammary luminal cell", 
                           "14-Dendritic cells" = "cDC",
                           "15-Hormone sensing" = "PIP+ mammary luminal cell", 
                           "16-Endothelial cells" = "Endothelial",
                           "17-Mammary epithelial cells" = "Cycling mammary luminal progenitor", 
                           "18-Macrophages" = "Macrophage",
                           "19-Alveolar" = "PIP+ mammary luminal cell", 
                           "20-NK cells" = "NK cell",
                           "21-T memory cells" = "CD4+ T cell", 
                           "22-Hormone sensing" = "PIP+ mammary luminal cell", 
                           "23-T memory cells" = "Regulatory T cell",
                           "24-Alveolar" = "SAA2+ mammary luminal progenitor", 
                           "25-Mammary epithelial cells" = "SCGB3A1 mammary luminal progenitor", 
                           "26-Alveolar" = "Secretoglobin mammary luminal cell", 
                           "27-Plasma cells" = "Plasma cell", 
                           "28-Alveolar" = "SAA2+ mammary luminal progenitor", 
                           "29-Alveolar" = "SAA2+ mammary luminal progenitor",
                           "30-Hormone sensing" = "PIP+ mammary luminal cell",
                           "31-Basal cells" = "Mammary basal cell",
                           "32-Basal cells" = "Mammary basal cell",
                           "33-Hormone sensing" = "PIP+ mammary luminal cell", 
                           "34-B cells naive" = "B cell", 
                           "35-Pericytes" = "Pericyte",
                           "36-Unknown" = "Secretoglobin mammary luminal cell",
                           "37-Hormone sensing" = "PIP+ mammary luminal cell",
                           "38-Adipocyte progenitor cells" = "Fibroblast",
                           "39-Hormone sensing" = "PIP+ mammary luminal cell",
                           "40-Mammary epithelial cells" = "Cycling mammary luminal progenitor",
                           "41-Hormone sensing" = "Secretoglobin mammary luminal cell",
                           "42-Plasmacytoid dendritic cells" = "pDC",
                           "43-Macrophages" = "Macrophage",
                           "44-Endothelial cells" = "Endothelial",
                           "45-Alveolar" = "SAA2+ mammary luminal progenitor",
                           "46-Mast cells" = "Mast cell",
                           "47-Fibroblasts" = "Fibroblast")

integrated$celltype <- Idents(integrated)
Idents(integrated) <- integrated$seurat_clusters

integrated[["S.Score"]] <- NULL
integrated[["G2M.Score"]] <- NULL

saveRDS(integrated, file = "data/integrated.rds")
# integrated <- readRDS(file = "data/integrated.rds")

# Visualisation----
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

subtype_palette <- c("#003f5c",  # Dark Blue
                     "#2f4b7c",  # Indigo
                     "#665191",  # Purple
                     "#a05195",  # Magenta
                     "#d45087",  # Pinkish Red
                     "#f95d6a",  # Coral
                     "#ff7c43")   # Orange)

study_palette <- c("#1A4383","#D5AC4C","#AE5B67")

# By BCa subtype (Figure 1b):
pdf("figures/integrated_by_type.pdf", width = 10, height = 6)
DimPlot(integrated, raster = F, group.by = "type", pt.size = 0.1, cols = subtype_palette) + NoAxes() + coord_fixed()
dev.off()

# By Seurat clusters (Supplementary Figure 1a):
pdf("figures/integrated_by_clusters.pdf", width = 4, height = 4)
DimPlot(integrated, label = T, raster = F, group.by = "seurat_clusters", pt.size = 0.1) + ggtitle("") + NoAxes() + coord_fixed() + NoLegend()
dev.off()

# By cell class (Supplementary Figure 1b):
pdf("figures/integrated_by_cellclass.pdf", width = 10, height = 6)
DimPlot(integrated, label = F, raster = F, group.by = "cellclass", pt.size = 0.1, cols = c("#4477AA","#228833","#CCBB44","#EE6677")) + NoAxes() + coord_fixed()+ NoLegend()+ggtitle("")
dev.off()

# By study (Supplementary Figure 1c):
pdf("figures/integrated_by_study.pdf", width = 10, height = 6)
DimPlot(integrated, label = T, repel = T, raster = F, group.by = "orig.ident", pt.size = 0.1, cols = study_palette) + NoAxes() + coord_fixed()
dev.off()

setwd(projectDir)

# Session----
sessionInfo()
