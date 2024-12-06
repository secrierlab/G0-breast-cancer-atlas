library(Seurat)
library(tidyverse)
library(ggpubr)

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

## function to create violin plots with statistics:
violin_plot_with_stats <- function(gene, test_sign){
  plot <- function(signature){
    VlnPlot(integrated, features = signature, assay = "RNA", pt.size = 0, raster = F,
            idents = c("G0 arrested", "Fast-cycling", "Slow-cycling"),
            cols = c("Fast-cycling" = "#1B7837", "G0 arrested" = "#762A83", "Slow-cycling" = "#f7f7f7"),
            y.max = round(max(integrated@assays$RNA@data[integrated@assays$RNA@data %>% rownames(.) == gene])) + 1
    ) + stat_compare_means(comparisons = test_sign, method = "t.test", hide.ns = T) + NoLegend() + theme(aspect.ratio = 1) + xlab("") + geom_boxplot()
  }
  purrr::map(gene, plot)
  file_name <- paste0("figures/", gene, "_violin.pdf")
  ggsave(file_name, width = 4, height = 4)
}

comparisons <- list(c("G0 arrested", "Fast-cycling"), c("G0 arrested", "Slow-cycling"), c("Fast-cycling", "Slow-cycling"))

## genes to plot (Figure 2d):
violin_plot_with_stats(gene = "MKI67", test_sign = comparisons)
violin_plot_with_stats(gene = "TOP2A", test_sign = comparisons)
violin_plot_with_stats(gene = "MALAT1", test_sign = comparisons)
violin_plot_with_stats(gene = "ATF3", test_sign = comparisons)

setwd(project_dir)
