library(Seurat)
library(dplyr)
library(ggplot2)

projectDir <- "~/BreastCancerG0arrest"
setwd("00_IntegratingDatasets")

# SCT Integration----
## Read Seurat Object----
seurat <- readRDS("data/seurat.merged.filtered.rds")

seurat.list <- SplitObject(seurat, split.by = "orig.ident")

## Estimate cell cycle scores to regress out----
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

## SCTransform each object----
seurat.list <- lapply(seurat.list, FUN = function(x) {
  x <-  CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  x$CC.Difference <- x$S.Score - x$G2M.Score
  x <- SCTransform(x, vars.to.regress = c("percent.mt", "CC.Difference"))
})

features <- SelectIntegrationFeatures(seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(seurat.list, anchor.features = features)
anchors <- FindIntegrationAnchors(seurat.list, normalization.method = "SCT", anchor.features = features)

# saveRDS(anchors, file = "./anchors.rds")
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

integrated$orig.ident[integrated$orig.ident == "Pal2021Group1"] <- "Pal2021"
integrated$orig.ident[integrated$orig.ident == "Pal2021Group2"] <- "Pal2021"
integrated$orig.ident[integrated$orig.ident == "Pal2021Group3"] <- "Pal2021"
integrated$orig.ident[integrated$orig.ident == "Pal2021Group4"] <- "Pal2021"
integrated$orig.ident[integrated$orig.ident == "Pal2021Group5"] <- "Pal2021"

integrated <- RunPCA(integrated)

# ElbowPlot(integrated, ndims = 50)
# integrated <- readRDS("data/SCTintegrated.rds")
pct <- integrated[["pca"]]@stdev / sum(integrated[["pca"]]@stdev) * 100
cumulative <- cumsum(pct)
co1 <- which(cumulative > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
nPCs <- min(co1, co2)
# nPCs <- 29

integrated <- FindNeighbors(integrated, dims = 1:nPCs, k.param = 15)

## Determine the resolution----
library(clustree)
seuratObj <- integrated
resolutionRange <- seq(from = 1.2, to = 2, by = 0.1)
seuratObj <- FindClusters(object = seuratObj, resolution = resolutionRange)
clustree(seuratObj, prefix = "integrated_snn_res.", highlight_core = T)
rm(seuratObj); gc(full = T)

## Find clusters----
# library(reticulate)
library(leidenAlg)
options(scipen = 1000000000)

integrated <- FindClusters(integrated, resolution = 2, method = "igraph", algorithm = "leiden", verbose = T)
integrated <- RunUMAP(integrated, dims = 1:nPCs, metric = "correlation", umap.method = "umap-learn")

DimPlot(integrated, label = T, repel = T, raster = F, group.by = "seurat_clusters", pt.size = 0.1) + coord_fixed()
DimPlot(integrated, label = T, repel = T, raster = F, group.by = "seurat_clusters", split.by = "seurat_clusters", ncol = 6, pt.size = 0.1) + coord_fixed()

saveRDS(integrated, file = "data/SCTintegrated_leiden.rds")

setwd(projectDir)
