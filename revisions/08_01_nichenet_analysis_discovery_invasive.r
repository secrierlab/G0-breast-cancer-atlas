library(nichenetr)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

# Set working directories
project_dir <- "~/BreastCancerG0arrest/"
working_dir <- paste0(project_dir, "00_revisions/")
figure_dir <- paste0(working_dir, "figures/")

# Define plot layout
layout <- "
#A
BC
"

setwd(working_dir)

# Process L-R databases for human
lr_network <- readRDS("data/lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS("data/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("data/weighted_networks_nsga2r_final.rds")

lr_network <- lr_network %>% distinct(from, to)
weighted_networks_lr <- weighted_networks$lr_sig %>%
  inner_join(lr_network, by = c("from", "to"))

## Load Seurat object----------------------------------------------------------
path <- "data/adata_full_harmony_integration_scenic.h5ad"
sce <- zellkonverter::readH5AD(path, reader = "R")
seurat_obj_full <- as.Seurat(sce, counts = "counts", data = NULL)
Idents(seurat_obj_full) <- seurat_obj_full$celltype

setwd(figure_dir)

# Shorten cell type names:
seurat_obj_full <- RenameIdents(seurat_obj_full, c(
  "PIP+ mammary luminal cell" = "Epithelial",
  "Endothelial" = "EC",
  "Cycling lactocyte" = "Epithelial",
  "Macrophage" = "Macrophage",
  "Secretoglobin mammary luminal cell" = "Epithelial",
  "Mammary basal cell" = "Epithelial",
  "SAA2+ mammary luminal progenitor" = "Epithelial",
  "Fibroblast" = "Fibroblast",
  "Cancer-associated fibroblast" = "CAF",
  "Tumour-associated macrophage" = "TAM",
  "Regulatory T cell" = "T_reg",
  "CD4+ T cell" = "CD4_T_cell",
  "Pericyte" = "Pericyte",
  "cDC" = "cDC",
  "NK/CD8+ T cell" = "CD8_T_cell",
  "pDC" = "pDC",
  "SCGB3A1+ mammary luminal progenitor" = "Epithelial",
  "Cycling mammary luminal progenitor" = "Epithelial",
  "Plasma cell" = "Plasma_cell",
  "B cell" = "B_cell",
  "Mast cell" = "Mast_cell"
))

seurat_obj_full$celltype <- Idents(seurat_obj_full)
# subset invasive breast carcinoma:
normal_patients <- c("Patient.12", "Patient.13", "Patient.14", "Patient.15")
subtypes <- c("ER", "PR", "HER", "TNBC")

seurat_obj <- subset(
  seurat_obj_full,
  subset = (type %in% subtypes) &
    !(seurat_obj_full$patient %in% normal_patients)
)

seurat_obj$celltype <- as.character(seurat_obj$celltype)
seurat_obj$celltype[seurat_obj$cycling_state_oren == "G0-arrested"] <- "G0_arrested"
seurat_obj$celltype[seurat_obj$cycling_state_oren == "Slow-cycling"] <- "Slow_cycling"
seurat_obj$celltype[seurat_obj$cycling_state_oren == "Fast-cycling"] <- "Fast_cycling"
Idents(seurat_obj) <- seurat_obj$celltype
seurat_obj$celltype_original <- seurat_obj$celltype

seurat_obj <- NormalizeData(seurat_obj)

pct <- 0.2

## Fast-cycling--------------------------------------------------------
receiver <- "Fast_cycling"
receiver_genes <- get_expressed_genes(receiver, seurat_obj, pct = pct)

all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, receiver_genes)

potential_ligands <- lr_network %>%
  filter(to %in% expressed_receptors) %>%
  pull(from) %>%
  unique()

sender_celltypes <- setdiff(
  unique(seurat_obj$celltype),
  c(receiver, "Slow_cycling", "G0_arrested")
)

# Get the expressed genes of every sender cell type separately
sender_gene_list <- sender_celltypes %>%
  unique() %>%
  lapply(get_expressed_genes, seurat_obj, pct)

sender_genes <- sender_gene_list %>%
  unlist() %>%
  unique()

potential_ligands_focused <- intersect(potential_ligands, sender_genes)

length(sender_genes)
length(potential_ligands)
length(potential_ligands_focused)

seurat_obj_receiver <- subset(
  seurat_obj,
  subset = celltype %in% c(receiver, "G0_arrested", "Slow_cycling")
)
seurat_obj_receiver <- NormalizeData(seurat_obj_receiver)

receiver_deg_table <- FindMarkers(
  object = seurat_obj_receiver,
  ident.1 = receiver,
  ident.2 = NULL,
  group.by = "celltype",
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.5,
  min.pct = pct
) %>% rownames_to_column("gene")

geneset_oi <- receiver_deg_table %>%
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.5) %>%
  pull(gene)

geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_genes <- receiver_genes %>% .[. %in% rownames(ligand_target_matrix)]

length(background_genes)
length(geneset_oi)

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands_focused
)

ligand_activities <- ligand_activities %>%
  arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected)))

p_hist_lig_activity <- ggplot(ligand_activities, aes(x = aupr_corrected)) +
  geom_histogram(color = "black", fill = "darkorange") +
  geom_vline(
    aes(
      xintercept = min(
        ligand_activities %>%
          top_n(25, aupr_corrected) %>%
          pull(aupr_corrected)
      ),
      color = "red",
      linetype = "dashed",
      linewidth = 1
    )
  ) +
  labs(x = "ligand activity (PCC)", y = "# ligands") +
  theme_classic() +
  theme(aspect.ratio = 1)

p_hist_lig_activity

# Focused approach on sender cell types:
ligand_activities <- ligand_activities %>%
  filter(test_ligand %in% potential_ligands_focused)

best_upstream_ligands <- ligand_activities %>%
  top_n(20, aupr_corrected) %>%
  arrange(-aupr_corrected) %>%
  pull(test_ligand) %>%
  unique()

ligand_aupr_matrix <- ligand_activities %>%
  filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>%
  select(aupr_corrected) %>%
  arrange(aupr_corrected)

vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1)

p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Prioritised ligands",
  "Ligand activity",
  legend_title = "AUPR\nligand activity",
  color = "brown"
) + theme(
  axis.text.x.top = element_blank(),
  aspect.ratio = nrow(vis_ligand_aupr) / ncol(vis_ligand_aupr),
  legend.position = "bottom",
  axis.title = element_blank(),
  axis.ticks.x = element_blank()
)

p_ligand_aupr

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = 200
  ) %>%
  bind_rows() %>%
  drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0
)

order_ligands <- intersect(
  best_upstream_ligands,
  colnames(active_ligand_target_links)
) %>% rev()

order_targets <- active_ligand_target_links_df$target %>%
  unique() %>%
  intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])
vis_ligand_aupr_to_plot <- vis_ligand_aupr[order_ligands, , drop = FALSE]
p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr_to_plot,
  "Prioritised ligands",
  "Ligand activity",
  legend_title = "AUPR\nligand activity",
  color = "brown"
) + theme(
  axis.text.x.top = element_blank(),
  aspect.ratio = nrow(vis_ligand_aupr_to_plot) / ncol(vis_ligand_aupr_to_plot),
  legend.position = "bottom",
  axis.title = element_blank(),
  axis.ticks.x = element_blank()
)

p_ligand_aupr

p_ligand_target <- make_heatmap_ggplot(
  vis_ligand_target,
  "Prioritised ligands",
  "Predicted target genes",
  color = "#2FBF71",
  legend_title = "Regulatory potential"
) + scale_fill_gradient2(low = "whitesmoke", high = "#2FBF71") +
  scale_x_discrete(position = "bottom") +
  theme(
    aspect.ratio = nrow(vis_ligand_target) / ncol(vis_ligand_target),
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

p_ligand_target

lr_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig
)

vis_lr_network <- prepare_ligand_receptor_visualization(
  lr_links_df,
  best_upstream_ligands,
  order_hclust = "both"
)

p_lr <- make_heatmap_ggplot(
  t(vis_lr_network),
  y_name = "Ligands",
  x_name = "Receptors",
  color = "mediumvioletred",
  legend_title = "Prior interaction potential"
) + theme(
  aspect.ratio = ncol(vis_lr_network) / nrow(vis_lr_network)
)

# Violin plot expression of target genes in receiver cells
exprs <- FetchData(
  seurat_obj_receiver %>%
    subset(celltype == receiver) %>%
    NormalizeData(),
  vars = order_targets
)

# Convert to long format
expr_long <- exprs %>%
  pivot_longer(everything(), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = order_targets))

# Plot all violins on one shared axis
p_target_expression <- ggplot(
  expr_long, aes(x = gene, y = expression, fill = gene)
) +
  geom_violin(scale = "width", trim = TRUE) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    aspect.ratio = 1 / length(unique(expr_long$gene))
  ) +
  ylab("Expression") +
  xlab("")

pdf(paste0("nichenet_", receiver, ".pdf"), width = 8, height = 8)
wrap_plots(
  A = p_target_expression,
  B = p_ligand_aupr,
  C = p_ligand_target,
  design = layout
) +
  plot_layout(
    widths = c(1, length(order_targets)),
    heights = c(1, length(order_ligands))
  )
dev.off()


## G0 arrested cells--------------------------------------------------------
receiver <- "G0_arrested"
receiver_genes <- get_expressed_genes(receiver, seurat_obj, pct = pct)

all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, receiver_genes)

potential_ligands <- lr_network %>%
  filter(to %in% expressed_receptors) %>%
  pull(from) %>%
  unique()

sender_celltypes <- setdiff(
  unique(seurat_obj$celltype),
  c(receiver, "Slow_cycling", "Fast_cycling")
)

# Get the expressed genes of every sender cell type separately
sender_gene_list <- sender_celltypes %>%
  unique() %>%
  lapply(get_expressed_genes, seurat_obj, pct)

sender_genes <- sender_gene_list %>%
  unlist() %>%
  unique()

potential_ligands_focused <- intersect(potential_ligands, sender_genes)

length(sender_genes)
length(potential_ligands)
length(potential_ligands_focused)

seurat_obj_receiver <- subset(
  seurat_obj,
  subset = celltype %in% c(receiver, "Fast_cycling", "Slow_cycling")
)
seurat_obj_receiver <- NormalizeData(seurat_obj_receiver)

receiver_deg_table <- FindMarkers(
  object = seurat_obj_receiver,
  ident.1 = receiver,
  ident.2 = NULL,
  group.by = "celltype",
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.5,
  min.pct = pct
) %>% rownames_to_column("gene")

geneset_oi <- receiver_deg_table %>%
  filter(p_val_adj <= 0.0001 & abs(avg_log2FC) >= 1) %>%
  pull(gene)

geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_genes <- receiver_genes %>% .[. %in% rownames(ligand_target_matrix)]

length(background_genes)
length(geneset_oi)

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands_focused
)

ligand_activities <- ligand_activities %>%
  arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected)))

p_hist_lig_activity <- ggplot(ligand_activities, aes(x = aupr_corrected)) +
  geom_histogram(color = "black", fill = "darkorange") +
  geom_vline(
    aes(
      xintercept = min(
        ligand_activities %>%
          top_n(30, aupr_corrected) %>%
          pull(aupr_corrected)
      ),
      color = "red",
      linetype = "dashed",
      linewidth = 1
    )
  ) +
  labs(x = "ligand activity (PCC)", y = "# ligands") +
  theme_classic() +
  theme(aspect.ratio = 1)

p_hist_lig_activity

# Focused approach on sender cell types:
ligand_activities <- ligand_activities %>%
  filter(test_ligand %in% potential_ligands_focused)

best_upstream_ligands <- ligand_activities %>%
  top_n(30, aupr_corrected) %>%
  arrange(-aupr_corrected) %>%
  pull(test_ligand) %>%
  unique()

ligand_aupr_matrix <- ligand_activities %>%
  filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>%
  select(aupr_corrected) %>%
  arrange(aupr_corrected)

vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1)

p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Prioritised ligands",
  "Ligand activity",
  legend_title = "AUPR\nligand activity",
  color = "brown"
) + theme(
  axis.text.x.top = element_blank(),
  aspect.ratio = nrow(vis_ligand_aupr) / ncol(vis_ligand_aupr),
  legend.position = "bottom",
  axis.title = element_blank(),
  axis.ticks.x = element_blank()
)

p_ligand_aupr

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = 50
  ) %>%
  bind_rows() %>%
  drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33
)

order_ligands <- intersect(
  best_upstream_ligands,
  colnames(active_ligand_target_links)
) %>% rev()

order_targets <- active_ligand_target_links_df$target %>%
  unique() %>%
  intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])
vis_ligand_aupr_to_plot <- vis_ligand_aupr[order_ligands, , drop = FALSE]
p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr_to_plot,
  "Prioritised ligands",
  "Ligand activity",
  legend_title = "AUPR\nligand activity",
  color = "brown"
) + theme(
  axis.text.x.top = element_blank(),
  aspect.ratio = nrow(vis_ligand_aupr_to_plot) / ncol(vis_ligand_aupr_to_plot),
  legend.position = "bottom",
  axis.title = element_blank(),
  axis.ticks.x = element_blank()
)

p_ligand_aupr

p_ligand_target <- make_heatmap_ggplot(
  vis_ligand_target,
  "Prioritised ligands",
  "Predicted target genes",
  color = "#6A5ACD",
  legend_title = "Regulatory potential"
) + scale_fill_gradient2(low = "whitesmoke", high = "#6A5ACD") +
  scale_x_discrete(position = "bottom") +
  theme(
    aspect.ratio = nrow(vis_ligand_target) / ncol(vis_ligand_target),
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

p_ligand_target

lr_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig
)

vis_lr_network <- prepare_ligand_receptor_visualization(
  lr_links_df,
  best_upstream_ligands,
  order_hclust = "both"
)

p_lr <- make_heatmap_ggplot(
  t(vis_lr_network),
  y_name = "Ligands",
  x_name = "Receptors",
  color = "mediumvioletred",
  legend_title = "Prior interaction potential"
) + theme(
  aspect.ratio = ncol(vis_lr_network) / nrow(vis_lr_network)
)

# Violin plot expression of target genes in receiver cells
exprs <- FetchData(
  seurat_obj_receiver %>%
    subset(celltype == receiver) %>%
    NormalizeData(),
  vars = order_targets
)

# Convert to long format
expr_long <- exprs %>%
  pivot_longer(everything(), names_to = "gene", values_to = "expression") %>%
  mutate(gene = factor(gene, levels = order_targets))

# Plot all violins on one shared axis
p_target_expression <- ggplot(
  expr_long, aes(x = gene, y = expression, fill = gene)
) +
  geom_violin(scale = "width", trim = TRUE) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    aspect.ratio = 1 / length(unique(expr_long$gene))
  ) +
  ylab("Expression") +
  xlab("")

pdf(paste0("nichenet_", receiver, ".pdf"), width = 12, height = 8)
wrap_plots(
  A = p_target_expression,
  B = p_ligand_aupr,
  C = p_ligand_target,
  design = layout
) +
  plot_layout(
    widths = c(1, length(order_targets)),
    heights = c(1, length(order_ligands))
  )
dev.off()
