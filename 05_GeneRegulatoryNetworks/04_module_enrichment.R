library(dplyr)
library(enrichR)
library(ggplot2)
library(patchwork)

project_dir <- "~/BreastCancerG0arrest/"
setwd(paste0(project_dir, "/07_GeneRegulatoryNetworks"))

# define a function for enrichment analysis
perform_enrichment <- function(module, db) {
  enriched <- enrichr(module, db)
  enriched <- lapply(names(enriched), function(name) {
    transform(enriched[[name]], GO = name)
  }) %>% 
    lapply(., function(x) x %>% slice(1:10)) %>% 
    bind_rows() %>% arrange(by = .$GO)
  enriched$GO <- as.factor(enriched$GO)
  enriched$FDR <- -log10(enriched$Adjusted.P.value)
  return(enriched)
}

# define the databases to use for enrichment analysis:
database <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023", "GO_Cellular_Component_2023")
fdr_cutoff = 7.5

# Fast-cycling GRN enrichment---------------------------------------------------
# Fast cycling GRN modules:
gene_modules <- list(
  M1 = c("JUN","ATF3","EGR1","KLF4","CEBPB","TFDP1","PRRX2","IRF6","FOXC2","RELB","RUNX3","DLX5","JUNB","CEBPD","SOX9","FOS","ELF3","MAFF","FOSL1"),
  M2 = c("FOXM1","MYBL1","MYBL2","E2F1","PITX1","ID4","ETV7"),
  M3 = c("DLX1","ALX4","MYC","DLX2","HMGA1","TFAP2A","WT1","HDAC2","ATF5","TEAD2","HEY1","STAT3","BCLAF1","BPTF","TBL1XR1"),
  M4 = c("IRF1","LTF","PDLIM5","STAT1","PBX1"),
  M5 = c("SPDEF","NR2F1","GATA2","PGR","TFAP2B","ARX","AR","FOSB","JUND","GATA3","ZNF217","NFIA","FOXO3","SP5","NKX2-8","ZNF91","ISL2","HOXB2","CHURC1","PAX9","IRF7","SOX2","ETV1","GBX2","LHX1","ZBTB18","ESR1","GATA5","HIF1A","FOXA1","CREB3L2","HOX5B","IRX5")
)

# perform enrichment analysis for all modules
enrichment_results <- lapply(names(gene_modules), function(module_name) {
  module <- gene_modules[[module_name]]
  enriched <- perform_enrichment(module, database)
  return(list(enriched = enriched, module_name = module_name))
})

# remove (GO: and some number  patterns from the GO terms
for (i in 1:length(enrichment_results)) {
  enrichment_results[[i]]$enriched$Term <- gsub("\\(GO:.*?\\)", "", enrichment_results[[i]]$enriched$Term)
  enrichment_results[[i]]$enriched$Term <- stringr::str_to_sentence(enrichment_results[[i]]$enriched$Term)
}

# generate a heat map for the enrichment results
enrichment_results_df <- lapply(enrichment_results, function(result) {
  result$enriched$module_name <- result$module_name
  return(result$enriched)}) %>% bind_rows() %>% filter(FDR > fdr_cutoff)

enrichment_results_df$GO <- factor(enrichment_results_df$GO, levels = unique(enrichment_results_df$GO))
enrichment_results_df$module_name <- factor(enrichment_results_df$module_name, levels = unique(enrichment_results_df$module_name))

# heat map:
heatmap <- ggplot(enrichment_results_df, aes(x = module_name, y = Term, fill = FDR)) +
  geom_tile() + 
  scale_fill_gradient2(low = "#2166AC", mid = "black", midpoint = fdr_cutoff,  high = "#B2182B") +
  theme_classic() + coord_fixed() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), text=element_text(size=6), legend.position = "bottom") +
  labs(x = "Module enrichment", y = "", fill = "FDR") + coord_fixed()

# add bars to the right side of the heat map that defines GO
bars <- ggplot(enrichment_results_df, aes(x = 1, y = Term, fill = GO)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#00a087", "#3c5488", "#f39b7f")) +
  theme_void() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(x = "", y = "") + coord_fixed(ratio = 4)

setwd("figures")
# plot the heat map (Extended Data Figure 5c):
pdf("gene_regulatory_network_modules_enrichment_analysis_Fast_cycling.pdf", width = 8, height = 6)
heatmap + bars + plot_layout(ncol = 2)
dev.off()

setwd(project_dir)