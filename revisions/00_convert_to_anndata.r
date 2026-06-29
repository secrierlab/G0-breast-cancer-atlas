library(Seurat)
library(SingleCellExperiment)

project_dir <- "~/BreastCancerG0arrest/"
working_dir <- paste0(project_dir, "00_revisions/")
figure_dir <- paste0(working_dir, "figures/")
data_dir <- paste0(working_dir, "data/")

setwd(working_dir)

seurat_obj <- readRDS(paste0(data_dir, "integrated_with_quiescence.rds"))
seurat_obj$QuiescenceType <- NULL
seurat_obj[["SCT"]] <- NULL
seurat_obj[["integrated"]] <- NULL
seurat_obj$quiescence_score <- seurat_obj$QuiescenceScore
seurat_obj$QuiescenceScore <- NULL

sce <- as.SingleCellExperiment(seurat_obj)

zellkonverter::writeH5AD(
    sce,
    paste0(data_dir, "integrated_revision.h5ad")
)
