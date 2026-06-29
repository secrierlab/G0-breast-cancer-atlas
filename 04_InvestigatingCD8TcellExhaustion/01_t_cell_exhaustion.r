library(Seurat)
library(tidyverse)
library(leidenAlg)
library(igraph)
library(ProjecTILs)
library(SignatuR)
library(clustree)
library(UCell)
project_dir <- "~/BreastCancerG0arrest"
# dir.create("figures")
setwd("04_InvestigatingCD8TcellExhaustion")

# Identify CD8+ T cell sub populations------------------------------------------
## Read in the data:
integrated <- readRDS(paste0(project_dir, "/02_G0arrestInMalignantCells/data/integrated_lite.rds"))

## Subset NK/CD8+ T cell population:
Idents(integrated) <- integrated$celltype
NK.CD8Tcells <- subset(integrated, subset = celltype == "NK/CD8+ T cell" & disease == "cancer")
NK.CD8Tcells[["SCT"]] <- NULL
NK.CD8Tcells[["integrated"]] <- NULL
# integration:
seurat.list <- SplitObject(NK.CD8Tcells, split.by = "orig.ident")
seurat.list$Gao2021 <- NULL

rm(integrated)

# SCTransform each object:
seurat.list <- lapply(seurat.list, FUN = function(x) {
  x <-  CellCycleScoring(x, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
  x$CC.Difference <- x$S.Score - x$G2M.Score
  x <- SCTransform(x, vars.to.regress = c("percent.mt", "CC.Difference"))
})

features <- SelectIntegrationFeatures(seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(seurat.list, anchor.features = features)
anchors <- FindIntegrationAnchors(seurat.list, anchor.features = features, normalization.method = "SCT")
NK.CD8Tcells <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Run PCA:
NK.CD8Tcells <- RunPCA(NK.CD8Tcells)
# Estimate PCs to include
pct <- NK.CD8Tcells[["pca"]]@stdev / sum(NK.CD8Tcells[["pca"]]@stdev) * 100
cumulative <- cumsum(pct)
co1 <- which(cumulative > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
nPCs <- min(co1, co2)
DimHeatmap(NK.CD8Tcells, dims = 1:15, balanced = T)
ElbowPlot(NK.CD8Tcells)
# Find Neighbours:
NK.CD8Tcells <- FindNeighbors(NK.CD8Tcells)
# Find Clusters:
## Determine the resolution:

seuratObj <- NK.CD8Tcells
resolutionRange <- seq(from = 0.3, to = 1, by = 0.1)
seuratObj <- FindClusters(object = seuratObj, resolution = resolutionRange, method = "igraph", algorithm = "leiden", verbose = T)
clustree(seuratObj, prefix = "integrated_snn_res.", highlight_core = T)
rm(seuratObj); gc(full = T)

NK.CD8Tcells <- FindClusters(NK.CD8Tcells, resolution = 0.6, method = "igraph", algorithm = "leiden", verbose = T)
# Run UMAP:
NK.CD8Tcells <- RunUMAP(NK.CD8Tcells, dims = 1:nPCs, metric = "correlation", umap.method = "umap-learn")

setwd("figures")
# Extended Data Fig. 7a:
pdf("CD8Tcells_clusters.pdf", width = 6, height = 6)
DimPlot(NK.CD8Tcells, group.by = "seurat_clusters", label = T) & NoAxes() & coord_fixed()
dev.off()
setwd("..")

# ProjectTILs reference mapping-------------------------------------------------
curl::curl_download("https://figshare.com/ndownloader/files/41414556", 
                    destfile = "data/CD8T_human_ref_v1.rds")
TICAtlas <- load.reference.map("data/CD8T_human_ref_v1.rds")

DefaultAssay(NK.CD8Tcells) <- "RNA"
NK.CD8Tcells <- NK.CD8Tcells %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c("percent.mt", "CC.Difference"))
# Run ProjecTILs:
NK.CD8Tcells <- Run.ProjecTILs(
  NK.CD8Tcells,
  ref = TICAtlas,
  filter.cells = F,
  split.by = "orig.ident",
  skip.normalize = T,
  reduction = "pca"
)

# Extended Data Fig. 7b (umap):
pdf("NK.CD8Tcells_projected_onto_reference_data.pdf", width = 6, height = 6)
plot.projection(TICAtlas, NK.CD8Tcells, linesize = 0.5, pointsize = 0.5) + 
  NoAxes() + theme(aspect.ratio = 1)
dev.off()

# Extended Data Fig. 7b (bars):
pdf("NK.CD8Tcells_projected_percentage.pdf", width = 4, height = 2)
plot.statepred.composition(TICAtlas, NK.CD8Tcells, metric = "Percent")
dev.off()

# Compare gene expressions:
genes4radar = c("FOXP3", "CD4", "CD8A", "TCF7", "CCR7", "SELL", "GZMB", "GZMK", "PDCD1", "HAVCR2", "TOX", "MKI67")

# Extended Data Fig. 7c:
pdf("NK.CD8Tcells_marker_genes.pdf", width = 12, height = 8)
plot.states.radar(TICAtlas, query = NK.CD8Tcells, genes4radar = genes4radar, min.cells = 20)
dev.off()

# Rename clusters:
NK.CD8Tcells$functional.cluster <- as.character(NK.CD8Tcells$functional.cluster)
NK.CD8Tcells$functional.cluster[NK.CD8Tcells$functional.cluster == "CD8.NaiveLike"] <- "Naive-like CD8+ T cell"
NK.CD8Tcells$functional.cluster[NK.CD8Tcells$functional.cluster == "CD8.CM"] <- "Naive CD8+ T cell"
NK.CD8Tcells$functional.cluster[NK.CD8Tcells$functional.cluster == "CD8.EM"] <- "Effector memory CD8+ T cell"
NK.CD8Tcells$functional.cluster[NK.CD8Tcells$functional.cluster == "CD8.TEMRA"] <- "Short-lived effector memory CD8+ T cell"
NK.CD8Tcells$functional.cluster[NK.CD8Tcells$functional.cluster == "CD8.TPEX"] <- "Precursor exhausted CD8+ T cell"
NK.CD8Tcells$functional.cluster[NK.CD8Tcells$functional.cluster == "CD8.TEX"] <- "Exhausted CD8+ T cell"
NK.CD8Tcells$functional.cluster[NK.CD8Tcells$functional.cluster == "CD8.MAIT"] <- "Mucosal-associated invariant CD8+ T cell"
NK.CD8Tcells$functional.cluster[is.na(NK.CD8Tcells$functional.cluster)] <- "Unknown"
NK.CD8Tcells$celltype <- NK.CD8Tcells$functional.cluster

# Remove redundant data:
NK.CD8Tcells$cellCycle.G1S <- NULL
NK.CD8Tcells$cellCycle.G2M <- NULL
NK.CD8Tcells$Tcell.cytotoxicity <- NULL
NK.CD8Tcells$Tcell.exhaustion <- NULL
NK.CD8Tcells$M2_macrophage <- NULL
NK.CD8Tcells$M1_macrophage <- NULL
NK.CD8Tcells$Tcell.stemness <- NULL
NK.CD8Tcells$IFN <- NULL
NK.CD8Tcells$S.Score <- NULL
NK.CD8Tcells$G2M.Score <- NULL
NK.CD8Tcells$functional.cluster <- NULL

# Stacked violin plot (Extended Data Fig. 7d):
pdf("immune_checkpoint_violin_stacked.pdf", width = 6, height = 4)
VlnPlot(NK.CD8Tcells, features = c("PDCD1", "CTLA4", "HAVCR2", "TIGIT", "LAG3"), group.by = "celltype", same.y.lims = T, pt.size = 0.1, fill.by = "ident", stack = T, flip = T) + theme(aspect.ratio = 1)
dev.off()

saveRDS(NK.CD8Tcells, "data/NK.CD8Tcells.rds")

setwd(project_dir)

sessionInfo()