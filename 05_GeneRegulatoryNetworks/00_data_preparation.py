import warnings
warnings.filterwarnings("ignore")

# setting up R dependencies
import anndata2ri
import rpy2
from rpy2.robjects import r
import random

anndata2ri.activate()

%load_ext rpy2.ipython

# Run this part of the code in R
# load the packages:
# library(Seurat)
# library(zellkonverter)
# library(SingleCellExperiment)
# set the working directory:
# projectDir <- "~/BreastCancerG0arrest"
# setwd(paste0(projectDir, "/05_GeneRegulatoryNetworks"))
# load the Seurat object:
# integrated <- readRDS(paste0(projectDir, "/02_G0arrestInMalignantCells/data/integrated_with_quiescence.rds"))
# extract pre-calculated UMAP and PCs:
# integrated$UMAP_1 <- integrated@reductions$umap@cell.embeddings[,"UMAP_1"]
# integrated$UMAP_2 <- integrated@reductions$umap@cell.embeddings[,"UMAP_2"]
# integrated$PC_1 <- integrated@reductions$pca@cell.embeddings[,"PC_1"]
# integrated$PC_2 <- integrated@reductions$pca@cell.embeddings[,"PC_2"]
# reduce the object size:
# integrated <- DietSeurat(integrated, assays = "RNA", layers = c("counts", "data"))
# rename cell types:
# integrated$celltype <- as.character(integrated$celltype)
# integrated$celltype[integrated$QuiescenceStatus == "Quiescent"] <- "G0_arrested"
# integrated$celltype[integrated$QuiescenceStatus == "Proliferating"] <- "Fast_cycling"
# integrated$celltype[integrated$QuiescenceStatus == "Slow-cycling"] <- "Slow_cycling"
# Idents(integrated) <- integrated$celltype
# subset celltypes of interest:
# integrated <- subset(integrated, idents = c("G0_arrested", "Fast_cycling", "Slow_cycling"))
# save the anndata object:
#  sce <- as.SingleCellExperiment(integrated)
# writeH5AD(sce, file = "resources/adata.h5ad")

# download the required databases
import requests
from multiprocessing.pool import ThreadPool

def url_response(url):
    path, url = url
    r = requests.get(url, stream = True)
    with open(path, "wb") as f:
        for ch in r:
            f.write(ch)

#Databases ranking the whole genome:            
urls = ["auxiliaries/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
         "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "auxiliaries/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather",
          "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather",
        "auxiliaries/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
         "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        "auxiliaries/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather",
         "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"]

for x in urls:
    url_response(x)

ThreadPool(5).imap_unordered(url_response, urls)

# Motif to TF annotations:
url = ["auxiliaries/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"]

for x in url:
    url_response(x)

ThreadPool(2).imap_unordered(url_response, url)