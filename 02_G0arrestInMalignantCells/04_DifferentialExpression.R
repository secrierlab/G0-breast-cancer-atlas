library(Seurat)
library(tidyverse)
library(ggpubr)

project_dir <- "~/BreastCancerG0arrest/"
setwd(paste0(project_dir, "/02_G0arrestInMalignantCells"))

# Differential expression analysis----------------------------------------------
## Read in & prepare the data---------------------------------------------------
integrated <- readRDS("data/integrated_with_quiescence.rds")

## Filter out lowly expressed genes
integrated <- integrated[, colMeans(integrated@assays$RNA@counts) > 0.1]

## Perform differential expression analysis
Idents(integrated) <- integrated$celltype

integrated_markers <- FindAllMarkers(
  integrated,
  group.by = "celltype",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

integrated_markers <- integrated_markers %>% filter(p_val_adj < 0.1)

## Save the results
dir.create("results")
# Supplementary Table 1:
write.csv(integrated_markers, "results/integrated_markers.csv")

## Heat map of top markers------------------------------------------------------
markers <- c("TFF3", "ACKR1", "TOP2A", "C1QC", "ISG15", "ACTG2", "FDCSP", "MFAP4",
             "TNFRSF4", "IL7R", "RGS5", "SPP1", "CPVL", "COL10A1", "GNLY", "CLEC4C",
             "SIX3", "CALML3", "IGLC2", "MS4A1", "TPSD1")

integrated <- ScaleData(integrated, features = markers)
Idents(integrated) <- integrated$celltype

## Fig. 1f
pdf("heatmap_by_clusters.pdf", width = 8, height = 6)
DoHeatmap(
  integrated,
  features = markers,
  group.by = "ident",
  group.bar = FALSE,
  assay = "RNA",
  raster = TRUE
) + NoLegend() + scale_fill_gradientn(colors = c("white", "whitesmoke", "red"))
dev.off()
## Differential expression for cell cycle states---------------------------------
integrated$celltype[integrated$QuiescenceStatus == "Quiescent"] <- "G0 arrested"
integrated$celltype[integrated$QuiescenceStatus == "Proliferating"] <- "Fast-cycling"
integrated$celltype[integrated$QuiescenceStatus == "Slow-cycling"] <- "Slow-cycling"
integrated <- subset(integrated, subset = celltype %in% c("G0 arrested", "Fast-cycling", "Slow-cycling"))
integrated$celltype <- factor(integrated$celltype, levels = c("Fast-cycling", "Slow-cycling", "G0 arrested"))
Idents(integrated) <- integrated$celltype

cell_cycle_markers <- FindAllMarkers(
  integrated,
  group.by = "celltype",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Filter out p.adj > 0.1
cell_cycle_markers <- cell_cycle_markers %>% filter(p_val_adj < 0.1)

## Save the results
# Supplementary Table 2:
write.csv(cell_cycle_markers, "results/cell_cycle_markers.csv")
