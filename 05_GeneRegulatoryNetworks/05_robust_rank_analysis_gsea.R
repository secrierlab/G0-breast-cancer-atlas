library(dplyr)
library(Seurat)
library(irGSEA)
library(doMC)
library(readxl)

project_dir <- "~/BreastCancerG0arrest"
setwd(paste0(project_dir, "/05_GeneRegulatoryNetworks"))

# Load the data:
seurat_obj <- readRDS("../02_G0arrestInMalignantCells/data/malignant_integrated.rds")
# subset invasive types:
seurat_obj <- subset(
  seurat_obj,
  subset = (disease %in% "cancer") &
    (seurat_obj$QuiescenceStatus %in% c("Quiescent", "Proliferating", "Slow-cycling")) &
    !(seurat_obj$patient %in% c("Patient.12", "Patient.13", "Patient.14", "Patient.15")) & 
    !(seurat_obj$type %in% c("DCIS", "neoplasm")))

counts <- GetAssayData(seurat_obj, layer = "counts")
meta_data <- seurat_obj@meta.data
umap_loadings <- seurat_obj@reductions$umap

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta_data)
seurat_obj@reductions$umap <- umap_loadings
remove(counts, meta_data)

Idents(seurat_obj) <- seurat_obj$QuiescenceStatus

seurat_obj <- seurat_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# factor the QuiescenceStatus for the plot
seurat_obj$QuiescenceStatus <- factor(
  seurat_obj$QuiescenceStatus,
  levels = c("Proliferating", "Slow-cycling", "Quiescent")
)

# read XLSX file
er_proteostasis_df <- read_xlsx("database/ER_Proteostasis_250624.xlsx", sheet = "Sheet1")
er_proteostasis <- split(er_proteostasis_df$`Gene Symbol`, er_proteostasis_df$Group)

proteostasis_network_df <- read_xlsx("database/Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06.xlsx", sheet = "dense")
proteostasis_network_gene_set <- split(proteostasis_network_df$`Gene Symbol`, proteostasis_network_df$Branch)

# combine gene sets
geneset_list <- c(er_proteostasis, proteostasis_network_gene_set)

seurat_obj <- irGSEA.score(
  object = seurat_obj,
  assay = "RNA",
  slot = "data",
  seeds = 123,
  ncores = 2,
  min.cells = 3,
  min.feature = 0,
  custom = T,
  geneset = geneset_list,
  msigdb = F,
  geneid = "symbol",
  aucell.MaxRank = 2000,
  ucell.MaxRank = 2000,
  kcdf = 'Gaussian'
)

result_dge <- irGSEA.integrate(
  object = seurat_obj,
  group.by = "QuiescenceStatus",
  metadata = NULL,
  col.name = NULL
)

write.csv(result_dge$RRA, "irGSEA_AUCell_upr.csv")
# remove non significant pathways
result_dge <- result_dge$RRA[result_dge$RRA$pvalue < 0.05, ]

setwd("figures")

pdf("irGSEA_heatmap_upr.pdf")
irGSEA.heatmap(object = result_dge,
               method = "RRA", 
               top = 50, 
               cluster.color = c("#1B7837","#BBBBBB","#762A83"),
               direction.color = c("#364B9A", "#A50026"),
               significance.color = c("grey", "gold3"),
               cluster.levels = c("Proliferating", "Slow-cycling", "Quiescent"),
               show.geneset = NULL)
dev.off()

pdf("irGSEA_bubble_plot_upr.pdf")
irGSEA.bubble(object = result_dge,
              method = "RRA",
              cluster.color = c("#1B7837","#BBBBBB","#762A83"),
              direction.color = c("#364B9A", "#A50026"),
              significance.color = c("#f7f7f7", "black"),
              cluster.levels = c("Proliferating", "Slow-cycling", "Quiescent"),
              top = 50)
dev.off()

pdf("irGSEA_barplot_upr.pdf", width = 4, height = 4)
irGSEA.barplot(object = result_dge, 
               method = "RRA", 
               significance.color = c("#364B9A", "grey", "#A50026"),
               color.cluster = c("#1B7837","#BBBBBB","#762A83"),
               cluster.levels = c("Proliferating", "Slow-cycling", "Quiescent"))
dev.off()

pathways <- result_dge$AUCell$Name

pdf("irGSEA_hubplot_upr_all.pdf")
for (pathway in pathways) {
  hub.result <- irGSEA.hub(
    object = seurat_obj,
    assay = "RNA",
    slot = "data",
    method = c("ssgsea"),
    show.geneset = pathway,
    ncores = 2,
    type = "rank",
    maxRank = 2000,
    top = 10,
    correlation.color = c("#364B9A", "#EAECCC", "#A50026"),
    method.color = NULL
  )
  
  print(hub.result$hub_result)
  print(hub.result$hub_plot)
}
dev.off()

# Mitochondrial and ribosomal translation---------------------------------------
translation <- geneset_list$Translation

# grep genes starting with MRP and RPS or RPL
translation_geneset <- list()
translation_geneset$Mitochondrial <- translation[grepl("^MRP", translation)]
translation_geneset$Ribosomal <- translation[grepl("^RPS", translation) | grepl("^RPL", translation)]

# run irGSEA
seurat_obj <- irGSEA.score(
  object = seurat_obj,
  assay = "RNA",
  slot = "data",
  seeds = 123,
  ncores = 2,
  min.cells = 3,
  min.feature = 0,
  custom = T,
  geneset = translation_geneset,
  msigdb = F,
  geneid = "symbol",
  aucell.MaxRank = 2000,
  ucell.MaxRank = 2000,
  kcdf = 'Gaussian'
)

result_dge <- irGSEA.integrate(
  object = seurat_obj,
  group.by = "QuiescenceStatus",
  metadata = NULL,
  col.name = NULL
)

pdf("irGSEA_densityheatmap_mitochondrial_translation.pdf")
irGSEA.densityheatmap(object = seurat_obj, 
                      method = "AUCell",
                      group.by = "QuiescenceStatus",
                      show.geneset = "Mitochondrial",
                      heatmap_width = 6,
                      heatmap_height = 5)
dev.off()

pdf("irGSEA_densityheatmap_ribosomal_translation.pdf")
irGSEA.densityheatmap(object = seurat_obj, 
                      method = "AUCell",
                      group.by = "QuiescenceStatus",
                      show.geneset = "Ribosomal",
                      heatmap_width = 6,
                      heatmap_height = 5)
dev.off()

pathways <- result_dge$AUCell$Name

pdf("irGSEA_hubplot_translation.pdf")
for (pathway in pathways) {
  hub.result <- irGSEA.hub(
    object = seurat_obj,
    assay = "RNA",
    slot = "data",
    method = c("AUCell"),
    show.geneset = pathway,
    ncores = 2,
    type = "rank",
    maxRank = 2000,
    top = 20,
    correlation.color = c("#364B9A", "#EAECCC", "#A50026"),
    method.color = NULL
  )
  
  print(hub.result$hub_result)
  print(hub.result$hub_plot)
}
dev.off()

setwd("..")
save(hub.result, file = "hub_result.RData")

sessionInfo()