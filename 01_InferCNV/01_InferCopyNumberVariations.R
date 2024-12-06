library(Seurat)
library(tidyverse)
library(ggplot2)
library(infercnv)
library(tibble)
library(tsvio)

options("preferRaster" = T) # if using leiden clustering (preferred)

project_dir <- "~/BreastCancerG0arrest"
setwd(paste0(project_dir, "01_InferCNV"))

# inferCNV for epithelia--------------------------------------------------------
## Read in data-----------------------------------------------------------------
integrated <- readRDS(paste0(project_dir, "/00_IntegratingDatasets/data/integrated.rds"))
umap.coordinates <- integrated@reductions$umap@cell.embeddings # for later
Idents(integrated) <- integrated$celltype

dir.create("figures/", recursive = F)
out_dir = "infercnv_results/"

## Subset reference (healthy) cells---------------------------------------------
healthy <- subset(integrated, subset = type == "normal" & cellclass == "EPI")
set.seed(1234) # for reproducibility 
healthy <- healthy[ , sample(colnames(healthy), size = 1000, replace = F)]
healthy$celltype <- paste0(healthy$celltype, "_normal")
idents <- unique(healthy$celltype)

epithelial <- subset(integrated, subset = cellclass == "EPI" & type != "normal")
epithelial.idents <- unique(epithelial$celltype)

epithelial.merged <- merge(x = healthy, y = epithelial)
saveRDS(epithelial.merged, file = paste0(out_dir, "epithelial.merged.rds"))
# epithelial.merged <- readRDS(paste0(out_dir, "epithelial.merged.rds"))
counts_matrix = GetAssayData(epithelial.merged, slot = "counts")
annotations_file <- cbind(colnames(counts_matrix), epithelial.merged$celltype)
write.table(annotations_file, "infercnv_epi/annotations_file_epi.txt", sep = "\t", row.names = F)

## Create inferCNV object----
infercnv_obj <- CreateInfercnvObject(counts_matrix, annotations_file = paste0(out_dir, "annotations_file_epi.txt"), gene_order_file = "database/hg38_gencode_v27.txt", ref_group_names = idents)

## Run inferCNV----------------------------------------------------------------
infercnv_obj_default <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, # recommended for 10X
  out_dir = out_dir,
  cluster_by_groups = F,
  analysis_mode = "subclusters",
  tumor_subcluster_partition_method = "leiden",
  tumor_subcluster_pval = 0.05,
  scale_data = T,
  leiden_resolution = 0.000005,
  plot_steps = F,
  denoise = T,
  k_obs_groups = 5,
  HMM = T,
  HMM_type = "i6",
  num_threads = 10,
  no_prelim_plot = T,
  output_format = "pdf",
  useRaster = F,
  plot_chr_scale = T,
  no_plot = T
)
# as the infercnv::run() generate extensive output, these files will not be
# shared in the repository, but using integrated_with_quiescence.rds in 
# 02_G0ArrestInMalignantCells/data, one could regenerate the outputs.

## Re-plot if necessary (Extended Data Fig. 1d)---------------------------------
infercnv::plot_cnv(
  infercnv_obj,
  out_dir = out_dir,
  cluster_by_groups = F,
  output_filename = "infercnv",
  x.range = "auto",
  chr_lengths = T,
  x.center = 1,
  title = "Copy number variations in epithelia",
  output_format = "pdf",
  custom_color_pal = color.palette(c("#2166ac", "#f7f7f7", "#B2182B"), c(2, 2)),
  color_safe_pal = T
)

setwd(project_dir)
