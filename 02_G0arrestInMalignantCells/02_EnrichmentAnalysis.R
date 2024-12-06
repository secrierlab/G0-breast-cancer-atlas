library(Seurat)
library(ggplot2)
library(ggpubr)
library(scDECAF)

# set working dir:
project_dir <- "~/BreastCancerG0arrest"
setwd(paste0(project_dir, "/02_G0arrestInMalignantCells/"))

# Enrichment Analysis on G0 arrested, Fast cycling and Slow cycling cells-------
## read in + prepare the data:
integrated <- readRDS("data/integrated_with_quiescence.rds")
integrated$celltype <- as.character(integrated$celltype)
integrated$celltype[integrated$QuiescenceStatus == "Quiescent"] <- "G0 arrested"
integrated$celltype[integrated$QuiescenceStatus == "Proliferating"] <- "Fast-cycling"
integrated$celltype[integrated$QuiescenceStatus == "Slow-cycling"] <- "Slow-cycling"
integrated <- subset(integrated, subset = (disease %in% "cancer") &
                       !(integrated$patient %in% c("Patient.12", "Patient.13", "Patient.14", "Patient.15")))
integrated <- subset(integrated, subset = celltype %in% c("G0 arrested", "Fast-cycling", "Slow-cycling"))
integrated$celltype <- factor(integrated$celltype, levels = c("Fast-cycling", "Slow-cycling", "G0 arrested"))
Idents(integrated) <- integrated$celltype

## load pathways:
hallmark_pathways <- msigdbr::msigdbr("human", category = "H")
hallmark_pathways_list <- split(hallmark_pathways$gene_symbol, hallmark_pathways$gs_name)

gene_ontology_terms <- msigdbr::msigdbr("human", category = "C5")
gene_ontology_terms_list <- split(gene_ontology_terms$gene_symbol, gene_ontology_terms$gs_name)

## expression matrix:
expression_matrix <- as.matrix(GetAssayData(integrated, assay = "SCT", layer = "data"))
## cell embeddings:
cell_embedding <- integrated@reductions$umap@cell.embeddings
## cell states:
outcomes <- integrated$QuiescenceStatus
## model matrix:
model_matrix <- as.matrix(model.matrix(~ 0 + outcomes))
encoding = ifelse(outcomes == "Quiescent", 1, 0)
cell_embedding_with_quiescence = cbind(cell_embedding, encoding)

## Hallmark pathways------------------------------------------------------------
selected_genesets <- pruneGenesets(data = expression_matrix, genesetlist = hallmark_pathways_list, hvg = rownames(expression_matrix), embedding = cell_embedding_with_quiescence, min_gs_size = 10, lambda = exp(-4.5))
target_hallmark_genesets <- genesets2ids(expression_matrix[match(rownames(expression_matrix), rownames(expression_matrix)), ], 
                                         hallmark_pathways_list[selected_genesets])
## run scDECAF for hallmark pathways
hallmark_annotation_results <- scDECAF(
  data = expression_matrix,
  gs = target_hallmark_genesets,
  standardize = FALSE,
  hvg = rownames(expression_matrix),
  k = 20,
  embedding = cell_embedding_with_quiescence,
  n_components = ncol(target_hallmark_genesets) - 1,
  max_iter = 2,
  thresh = 0.5
)

hallmark_scores = data.frame(attributes(hallmark_annotation_results)$raw_scores)

# add metadata back to seurat object
integrated <- AddMetaData(integrated, metadata = hallmark_scores)
setwd("figures")

comparisons <- list(c("G0 arrested", "Fast-cycling"), c("G0 arrested", "Slow-cycling"), c("Fast-cycling", "Slow-cycling"))

# Hallmark EMT (Figure 4d):
pdf("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf", width = 4, height = 4)
VlnPlot(integrated, features = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", cols = c("#1B7837","#F7F7F7","#762A83"), pt.size=0, y.max = 1.5) + 
  NoLegend() + geom_boxplot(varwidth = T, notch = T) + theme(aspect.ratio = 1) +
  stat_kruskal_test(label.x.npc = "centre", p.adjust.method = "BH") + xlab("") +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons)
dev.off()

# Hallmark p53 (Figure 3f):
pdf("HALLMARK_P53_PATHWAY.pdf", width = 4, height = 4)
VlnPlot(integrated, features = "HALLMARK_P53_PATHWAY", cols = c("#1B7837","#F7F7F7","#762A83"), pt.size=0, y.max = 1.5) + 
  NoLegend() + geom_boxplot(varwidth = T, notch = T) +
  stat_kruskal_test(label.x.npc = "centre", p.adjust.method = "BH") + theme(aspect.ratio = 1) + xlab("") +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons)
dev.off()

setwd("..")

# Cell cycle pathways-----------------------------------------------------------
# filter cell cycle gene sets:
cell_cycle_genesets <- gene_ontology_terms_list[grepl("CELL_CYCLE", names(gene_ontology_terms_list))]
cell_cycle_genesets <- c(cell_cycle_genesets, "negative regulation of cyclin-dependent protein serine/threonine kinase activity")
# prune gene sets:
cell_cycle_selected_gene_sets <- pruneGenesets(data = expression_matrix, genesetlist = cell_cycle_genesets, hvg = rownames(expression_matrix), embedding = cell_embedding_with_quiescence, min_gs_size = 10, lambda = exp(-4.5))
target_cell_cycle_genesets <- genesets2ids(expression_matrix[match(rownames(expression_matrix), rownames(expression_matrix)),], cell_cycle_genesets[cell_cycle_selected_gene_sets])
# run scDECAF:
cell_cycle_annotation_results <- scDECAF(data = expression_matrix, gs = target_cell_cycle_genesets, standardize = F, hvg = rownames(expression_matrix), k = 20, embedding = cell_embedding_with_quiescence, n_components = ncol(target_cell_cycle_genesets) - 1, max_iter = 2, thresh = 0.5)

cell_cycle_scores <- data.frame(attributes(cell_cycle_annotation_results)$raw_scores)

# add metadata back to seurat object
integrated <- AddMetaData(integrated, metadata = cell_cycle_scores)
# plot results
setwd("figures")

# Figure 5d:
pdf("GOBP_CELL_CYCLE_DNA_REPLICATION.pdf", width = 4, height = 4)
VlnPlot(integrated, features = "GOBP_CELL_CYCLE_DNA_REPLICATION", cols = c("#1B7837","#F7F7F7","#762A83"), pt.size=0, y.max = 1.5) + NoLegend() + geom_boxplot(varwidth = T, notch = T) + 
  stat_kruskal_test(label.x.npc = "centre", p.adjust.method = "BH") + theme(aspect.ratio = 1) + xlab("") +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons)
dev.off()

# Figure 3f:
pdf("GOBP_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR_RESULTING_IN_CELL_CYCLE_ARREST.pdf", width = 4, height = 4)
VlnPlot(integrated, features = "GOBP_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR_RESULTING_IN_CELL_CYCLE_ARREST", cols = c("#1B7837","#F7F7F7","#762A83"), pt.size=0, y.max = 1.5) + NoLegend() + geom_boxplot(varwidth = T, notch = T) + 
  stat_kruskal_test(label.x.npc = "centre", p.adjust.method = "BH") + theme(aspect.ratio = 1) + xlab("") +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons)
dev.off()

# Figure 3f:
pdf("GOBP_REGULATION_OF_TRANSCRIPTION_INVOLVED_IN_G1_S_TRANSITION_OF_MITOTIC_CELL_CYCLE.pdf", width = 4, height = 4)
VlnPlot(integrated, features = "GOBP_REGULATION_OF_TRANSCRIPTION_INVOLVED_IN_G1_S_TRANSITION_OF_MITOTIC_CELL_CYCLE", cols = c("#1B7837","#F7F7F7","#762A83"), pt.size=0, y.max = 1.5) + NoLegend() + geom_boxplot(varwidth = T, notch = T) + 
  stat_kruskal_test(label.x.npc = "centre", p.adjust.method = "BH") + theme(aspect.ratio = 1) + xlab("") +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons)
dev.off()

# Figure 3f:
pdf("GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION.pdf", width = 4, height = 4)
VlnPlot(integrated, features = "GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION", cols = c("#1B7837","#F7F7F7","#762A83"), pt.size=0, y.max = 1.5) + NoLegend() + geom_boxplot(varwidth = T, notch = T) + 
  stat_kruskal_test(label.x.npc = "centre", p.adjust.method = "BH") + theme(aspect.ratio = 1) + xlab("") +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons)
dev.off()

setwd("..")

# Unfolded protein response pathways--------------------------------------------
# filter UPR gene sets:
upr_genesets <- gene_ontology_terms_list[grepl("UNFOLDED_PROTEIN", names(gene_ontology_terms_list))]
# prune gene sets:
upr_selected_gene_sets <- pruneGenesets(data = expression_matrix, genesetlist = upr_genesets, hvg = rownames(expression_matrix), embedding = cell_embedding_with_quiescence, min_gs_size = 10, lambda = exp(-4.5))
target_upr_genesets <- genesets2ids(expression_matrix[match(rownames(expression_matrix), rownames(expression_matrix)),], upr_genesets[upr_selected_gene_sets])
# run scDECAF:
upr_annotation_results <- scDECAF(data = expression_matrix, gs = target_upr_genesets, standardize = F, hvg = rownames(expression_matrix), k = 20, embedding = cell_embedding_with_quiescence, n_components = ncol(target_upr_genesets) - 1, max_iter = 2, thresh = 0.5)

upr_scores <- data.frame(attributes(upr_annotation_results)$raw_scores)

# add metadata back to seurat object
integrated <- AddMetaData(integrated, metadata = upr_scores)
# plot results
setwd("figures")

# Figure 5d:
pdf("GOBP_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE.pdf", width = 4, height = 4)
VlnPlot(integrated, features = "GOBP_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE", cols = c("#1B7837","#F7F7F7","#762A83"), pt.size=0, y.max = 1.5) + 
  NoLegend() + geom_boxplot(varwidth = T, notch = T) + 
  stat_kruskal_test(label.x.npc = "centre", p.adjust.method = "BH") + theme(aspect.ratio = 1) + xlab("") +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons)
dev.off()

# ssGSEA------------------------------------------------------------------------
library(escape)
library(dittoSeq)
library(fastcluster)

setwd("..")
integrated <- subset(integrated, idents = c("G0 arrested", "Fast-cycling"))

integrated <- integrated[, sample(colnames(integrated), size = 4000, replace = F)]
gene.sets <- getGeneSets(library = "H")
enrichment_scores <- escape.matrix(integrated,
                                   gene.sets = gene.sets, 
                                   groups = 1000, 
                                   min.size = 10)

integrated <- AddMetaData(integrated, enrichment_scores)

integrated@meta.data$active.idents <- integrated@active.ident

setwd("figures")

pathways <- c(
  "HALLMARK-E2F-TARGETS",
  "HALLMARK-G2M-CHECKPOINT",
  "HALLMARK-MITOTIC-SPINDLE",
  "HALLMARK-MYC-TARGETS-V1",
  "HALLMARK-MYC-TARGETS-V2",
  "HALLMARK-KRAS-SIGNALING-DN",
  "HALLMARK-OXIDATIVE-PHOSPHORYLATION",
  "HALLMARK-MTORC1-SIGNALING",
  "HALLMARK-P53-PATHWAY",
  "HALLMARK-PI3K-AKT-MTOR-SIGNALING",
  "HALLMARK-PROTEIN-SECRETION",
  "HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY",
  "HALLMARK-DNA-REPAIR",
  "HALLMARK-UNFOLDED-PROTEIN-RESPONSE")

# Figure 2e:
pdf("hallmarks_heatmap.pdf", width = 6, height = 3)
dittoHeatmap(
  integrated,
  genes = NULL,
  metas = pathways,
  annot.colors = c(
    "#3375B5",
    "#5FB8E9",
    "#439788",
    "#DF7E44",
    "#BC4024",
    "#DB4677",
    "#1B7837",
    "#762A83"
  ),
  annot.by = c("type", "celltype"),
  cluster_cols = T,
  cluster_rows = T,
  fontsize = 6,
  scale = "row",
  heatmap.colors = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(25))
dev.off()

setwd(project_dir)
