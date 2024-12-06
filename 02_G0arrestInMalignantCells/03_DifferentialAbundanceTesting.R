library(Seurat)
library(tidyverse)
library(ggpubr)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)

project_dir <- "~/BreastCancerG0arrest/"
setwd(paste0(project_dir, "/02_G0arrestInMalignantCells"))

# Read in & prepare the data----
integrated <- readRDS("data/integrated_with_quiescence.rds")
integrated$celltype[integrated$QuiescenceStatus == "Quiescent"] <- "G0 arrested"
integrated$celltype[integrated$QuiescenceStatus == "Proliferating"] <- "Fast-cycling"
integrated$celltype[integrated$QuiescenceStatus == "Slow-cycling"] <- "Slow-cycling"
integrated <- subset(integrated, subset = celltype %in% c("G0 arrested", "Fast-cycling", "Slow-cycling"))
integrated$celltype <- factor(integrated$celltype, levels = c("Fast-cycling", "Slow-cycling", "G0 arrested"))
Idents(integrated) <- integrated$celltype

# Differential abundance testing----
malignant <- readRDS("data/epi_integrated.rds")
## create milo object:
malignant <- subset(malignant, subset = malignancy %in% "malignant")
malignant <- subset(malignant, subset = celltype %in% c("G0 arrested", "Fast-cycling"))
Idents(malignant) <- malignant$celltype
DefaultAssay(malignant) <- "RNA"
malignant <- DietSeurat(malignant, assays = "RNA", dimreducs = c("pca", "umap"))
malignant <- JoinLayers(malignant)

sce <- as.SingleCellExperiment(malignant)

# add pre-calculated dimensional reduction to SingleCellExperiment object:
reducedDim(sce, "pca", withDimnames = T) <- malignant[['pca']]@cell.embeddings
reducedDim(sce, "umap", withDimnames = T) <- malignant[['umap']]@cell.embeddings

## convert to milo object:
milo <- Milo(sce)
milo <- milo[rowSums(logcounts(milo)) != 0, ]
## build graph:
d = 29
k = 15
milo <- buildGraph(milo, k = k, d = d, reduced.dim = "pca", transposed = T)
## make neighbourhoods:
milo <- makeNhoods(milo, prop = 0.2, k = k, d = d, refined = T, reduced_dims = "pca")
plotNhoodSizeHist(milo)
# rename patients to make them unique:
colData(milo)$patient <- paste(colData(milo)$orig.ident, colData(milo)$celltype, colData(milo)$type, sep = "_")
## count cells:
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample = "patient")
head(nhoodCounts(milo))

design <- data.frame(colData(milo))[,c("patient", "celltype", "type")]
design <- distinct(design)
rownames(design) <- design$patient
## Reorder row names to match columns of nhoodCounts(milo)
design <- design[colnames(nhoodCounts(milo)), , drop = F]

## calculate neighbour distances:
milo <- calcNhoodDistance(milo, d = d, reduced.dim = "pca")
## differential abundance testing:
differential_abundance_results <- testNhoods(milo, design = ~celltype, design.df = design, reduced.dim = "pca", robust = T)
head(differential_abundance_results)

differential_abundance_results %>% arrange(SpatialFDR) %>% head()

## visualize results p-values:
ggplot(differential_abundance_results, aes(PValue)) + geom_histogram(bins = 50)
## visualise differential expression:
ggplot(differential_abundance_results, aes(logFC, -log10(SpatialFDR))) + geom_point() + geom_hline(yintercept = 1)
## Build neighbourhood graph:
milo <- buildNhoodGraph(milo)

setwd("figures")

## annotate neighbourhoods:
differential_abundance_results <- annotateNhoods(milo, differential_abundance_results, coldata_col = "celltype")
head(differential_abundance_results)
## visualise fractions:
ggplot(differential_abundance_results, aes(celltype_fraction)) + geom_histogram(bins=50)
## set cut-off:
differential_abundance_results$celltype <- ifelse(differential_abundance_results$celltype_fraction < 0.7, "Mixed", differential_abundance_results$celltype)

# (Extended Data Fig. 2b)
pdf("differential_abundance_results.pdf", height = 8, width = 8)
plotDAbeeswarm(differential_abundance_results, group.by = "celltype") +
  scale_colour_gradient2(high = "#762A83", low = "#1B7837", mid = "#f7f7f7")
dev.off()

## plot neighbourhood graph (Extended Data Fig. 2c):
pdf("neighbourhoods_umap.pdf", width = 6, height = 6)
plotNhoodGraphDA(milo, differential_abundance_results, layout = "umap", alpha = 0.1) + 
  scale_fill_gradient2(high = "#762A83", low = "#1B7837", mid = "#f7f7f7") + 
  coord_fixed()
dev.off()

setwd(project_dir)
