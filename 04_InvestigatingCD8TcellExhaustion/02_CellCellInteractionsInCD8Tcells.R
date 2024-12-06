library(Seurat)
library(tidyverse)
library(nichenetr)
project_dir <- "~/BreastCancerG0arrest"
setwd(paste0(project_dir, "/04_InvestigatingCD8TcellExhaustion"))

# Cell-cell interactions in CD8+ T cells and its subtypes-----------------------
subtypes <- c("ER", "HER", "PR", "TNBC")
## Load in the databases:
lr_network <- readRDS(paste0(project_dir, "/03_CellCellInteractions/database/lr_network_human_21122021.rds")) %>% distinct(from, to)
ligand_target_matrix <- readRDS(paste0(project_dir, "/03_CellCellInteractions/database/ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(paste0(project_dir, "/03_CellCellInteractions/database/weighted_networks_nsga2r_final.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

## Load the CD8T cell data:
NK.CD8Tcells <- readRDS("data/NK.CD8Tcells_with_annotations.rds")
DefaultAssay(NK.CD8Tcells) <- "RNA"
NK.CD8Tcells <- DietSeurat(NK.CD8Tcells, assays = "RNA")
NK.CD8Tcells <- alias_to_symbol_seurat(NK.CD8Tcells, "human")
NK.CD8Tcells$celltype[is.na(NK.CD8Tcells$celltype)] <- "Unknown"
NK.CD8Tcells$celltype[(NK.CD8Tcells$celltype %in% c("Naive CD8+ T cell"))] <- "Reference"
Idents(NK.CD8Tcells) <- NK.CD8Tcells$celltype

## Load the integrated data:
integrated <- readRDS(paste0(project_dir, "/02_G0arrestInMalignantCells/data/integrated_with_quiescence.rds"))
integrated <- DietSeurat(integrated, assays = "RNA")
integrated <- alias_to_symbol_seurat(integrated, organism = "human")
Idents(integrated) <- integrated$celltype

# Shorten cell type names:
integrated <- RenameIdents(integrated, c("PIP+ mammary luminal cell" = "PIP.luminal",
                                         "Endothelial" = "EC",
                                         "Cycling lactocyte" = "Cyc.Lac",
                                         "Macrophage" = "Mac",
                                         "Secretoglobin mammary luminal cell" = "Sec.luminal",
                                         "Mammary basal cell" = "Basal",
                                         "SAA2+ mammary luminal progenitor" = "SAA2.luminal",
                                         "Fibroblast" = "Fib",
                                         "Cancer-associated fibroblast" = "CAF",
                                         "Tumour-associated macrophage" = "TAM",
                                         "Regulatory T cell" = "Treg",
                                         "CD4+ T cell" = "CD4.T",
                                         "Pericyte" = "Peri",
                                         "cDC" = "cDC",
                                         "NK/CD8+ T cell" = "NK/CD8.T",
                                         "pDC" = "pDC",
                                         "SCGB3A1+ mammary luminal progenitor" = "SCGB3A1.Pro",
                                         "Cycling mammary luminal progenitor" = "Cyc.Pro",
                                         "Plasma cell" = "Plas.B",
                                         "B cell" = "Bcell",
                                         "Mast cell" = "Mast"))

integrated$celltype <- Idents(integrated)
integrated <- subset(integrated, subset = (disease %in% "cancer") & 
                       !(integrated$patient %in% c("Patient.12", "Patient.13", "Patient.14", "Patient.15")) & 
                       !(integrated$type %in% c("DCIS", "neoplasm")))

## G0 arrested cells----
integrated$celltype <- as.character(integrated$celltype)
integrated$celltype[integrated$QuiescenceStatus == "Quiescent"] <- "G0_arrested"
integrated$celltype[integrated$QuiescenceStatus == "Proliferating"] <- "Proliferating"
integrated$celltype[integrated$QuiescenceStatus == "Slow-cycling"] <- "Reference"
Idents(integrated) <- integrated$celltype
integrated$celltype_original <- integrated$celltype

setwd("figures")

# Exhausted CD8+ T cell-cell interactions---------------------------------------
## Define receiver cell population:
receiver = "Exhausted CD8+ T cell"
pct = 0.1
Idents(NK.CD8Tcells) <- NK.CD8Tcells$celltype
expressed_genes_receiver = get_expressed_genes(receiver, NK.CD8Tcells, pct = pct)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## Define a gene set of interest:
DE_table_receiver = FindMarkers(object = NK.CD8Tcells, ident.1 = receiver, ident.2 = "Reference", group.by = "celltype", min.pct = pct, only.pos = T, test.use = "bimod") %>% rownames_to_column("gene")
geneset_of_interest = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.5) %>% pull(gene)
geneset_of_interest = geneset_of_interest %>% .[. %in% rownames(ligand_target_matrix)]

## Define sender cell populations:
sender_celltypes <- c("Mac", "Treg", "CD4.T", "cDC", "pDC", "TAM", "CAF", "Fib", "EC", "Peri", "Bcell", "Mast", "Plas.B", "G0_arrested", "Proliferating")
list_expressed_genes_sender = sender_celltypes %>% lapply(get_expressed_genes, integrated, pct)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

## Define a set of potential ligands:
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands, expressed_genes_sender)
expressed_receptors = intersect(receptors, expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# Perform NicheNet ligand activity analysis:
ligand_activities = predict_ligand_activities(geneset = geneset_of_interest, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

## Infer receptors and top-predicted target genes of ligands:
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_of_interest, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritised ligands", "Predicted target genes", color = "#F5D39F", legend_position = "top", x_axis_position = "top", legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "#F7F7F7",  high = "#F5D39F")

# Receptor of top-ranked ligands:
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight", fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
# receptors:
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
# ligands:
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "#EE3377", x_axis_position = "top", legend_title = "Prior interaction potential")

# Expression of target genes per tumour type:
Tex <- subset(NK.CD8Tcells, subset = celltype %in% receiver & type %in% subtypes)
order_targets_adopted <- str_replace_all(order_targets, "\\.", "-")
Idents(Tex) <- Tex$type
rotated_dotplot_target_gene_expression <- DotPlot(Tex, features = order_targets_adopted, scale = F) + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12), axis.text.x = element_text(angle = 90, hjust = 1)) + scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B")

## Summary visualisation of the NicheNet analysis----
# ligand activity heatmap
ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()
vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritised ligands", "Ligand activity", color = "brown",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)") + theme(legend.text = element_text(size = 9))
order_ligands_adapted <- str_replace_all(order_ligands, "\\.", "-")
rotated_dotplot = DotPlot(integrated %>% subset(celltype %in% sender_celltypes), features = order_ligands_adapted) + scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12))

rotated_dotplot_target_gene_expression <- VlnPlot(Tex, group.by = "celltype", features = order_targets_adopted, stack = T, pt.size = 0, cols = rep("#F5D39F", length(order_targets_adopted))) & theme_void() & NoLegend() & coord_flip()

pdf("TME_to_Texhausted_summary_ligand_receptor_interactions.pdf", width = 35, height = 15)
figures_without_legend = cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 8), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  NULL,
  NULL,
  rotated_dotplot_target_gene_expression + theme(legend.position = "none", axis.ticks = element_blank()) + xlab(""),
  align = "h",
  nrow = 2,
  rel_widths = c(ncol(vis_ligand_aupr) + 6, length(sender_celltypes) + 8, ncol(vis_ligand_target)),
  rel_heights = c(nrow(vis_ligand_aupr), nrow(rotated_dotplot_target_gene_expression)))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot_target_gene_expression)),
  ncol = 4,
  align = "h")

cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10, 2), nrow = 2, align = "hv")
dev.off()

setwd(project_dir)
sessionInfo()