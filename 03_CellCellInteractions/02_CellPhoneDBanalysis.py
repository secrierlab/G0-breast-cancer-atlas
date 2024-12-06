# Prepare data for CellPhoneDB--------------------------------------------------
import warnings
warnings.filterwarnings("ignore")

# setting up R dependencies
import anndata2ri
import rpy2
from rpy2.robjects import r
import random

anndata2ri.activate()

%load_ext rpy2.ipython

# In R, do the following to prepare the data for CellPhoneDB v5:
# load the packages
# library(Seurat)
# library(zellkonverter)
# library(SingleCellExperiment)
# set the working directory
# projectDir <- "~/BreastCancerG0arrest"
# setwd(paste0(projectDir, "/03_CellCellInteractions"))
# load the Seurat object
# integrated <- readRDS(paste0(projectDir, "/02_G0arrestInMalignantCells/data/integrated_with_quiescence.rds"))
# reduce the object size
# integrated <- DietSeurat(integrated, assays = "RNA", layers = c("counts", "data"))
# Shorten cell type names:
# Idents(integrated) <- integrated$celltype
# integrated <- RenameIdents(integrated, c("PIP+ mammary luminal cell" = "PIP.luminal",
#                                          "Endothelial" = "EC",
#                                          "Cycling lactocyte" = "Cyc.Lac",
#                                          "Macrophage" = "Mac",
#                                          "Secretoglobin mammary luminal cell" = "Sec.luminal",
#                                          "Mammary basal cell" = "Basal",
#                                          "SAA2+ mammary luminal progenitor" = "SAA2.luminal",
#                                          "Fibroblast" = "Fib",
#                                          "Cancer-associated fibroblast" = "CAF",
#                                          "Tumour-associated macrophage" = "TAM",
#                                          "Regulatory T cell" = "Treg",
#                                          "CD4+ T cell" = "CD4.T",
#                                          "Pericyte" = "Peri",
#                                          "cDC" = "cDC",
#                                          "NK/CD8+ T cell" = "CD8.T",
#                                          "pDC" = "pDC",
#                                          "SCGB3A1+ mammary luminal progenitor" = "SCGB3A1.Pro",
#                                          "Cycling mammary luminal progenitor" = "Cyc.Pro",
#                                          "Plasma cell" = "Plas.B",
#                                          "B cell" = "Bcell",
#                                          "Mast cell" = "Mast"))

# rename cell types
# integrated$celltype <- Idents(integrated)
# integrated$celltype <- as.character(integrated$celltype)
# integrated$celltype[integrated$QuiescenceStatus == "Quiescent"] <- "G0_arrested"
# integrated$celltype[integrated$QuiescenceStatus == "Proliferating"] <- "Fast_cycling"
# integrated$celltype[integrated$QuiescenceStatus == "Slow-cycling"] <- "Slow_cycling"
# Idents(integrated) <- integrated$celltype

# subset the object to only include cancer patients
# integrated <- subset(integrated, subset = (disease %in% "cancer") & !(integrated$patient %in% c("Patient.12", "Patient.13", "Patient.14", "Patient.15")))
# subset the object to only include invasive breast cancer
# integrated <- subset(integrated, subset = type %in% c("ER", "PR", "HER", "TNBC"))

# save the metadata
# metadata <- integrated@meta.data
# metadata$barcode_sample <- rownames(metadata)
# metadata <- metadata[, c("barcode_sample", "celltype")]
# colnames(metadata) <- c("barcode_sample", "cell_type")
# write.table(metadata, file = "data/tumours_metadata.tsv", sep = "\t", quote = F, row.names = F)
# prepare the data for CellPhoneDB
# integrated$cell_type <- Idents(integrated)
# find DEGs for each cell type (for only method 3 of cellphonedb)
# DEGs <- FindAllMarkers(integrated, only.pos = T, min.pct = 0.1)
# DEGs_reordered <- DEGs[, c(6, 7, 5, 1, 2, 3, 4)]
# colnames(DEGs_reordered) <- c("cell_type", "gene", "p_val_adj", "p_val", "avg_log2FC", "pct.1", "pct.2")
# write.table(DEGs_reordered, file = "data/tumours_DEGs.tsv", row.names = F, sep = "\t", quote = F)
# save the anndata object
# sce <- as.SingleCellExperiment(integrated)
# writeH5AD(sce, file = "data/tumours.h5ad")

import pandas as pd
import os

import anndata as ad
import ktplotspy as kpy
import matplotlib.pyplot as plt
%matplotlib inline
plt.rcParams["pdf.fonttype"] = "truetype"

from IPython.display import HTML, display
from cellphonedb.utils import db_releases_utils

display(HTML(db_releases_utils.get_remote_database_versions_html()["db_releases_html_table"]))

# -- Version of the databse
cpdb_version = "v5.0.0"

# -- Path where the input files to generate the database are located
cpdb_target_dir = os.path.join("03_CellCellInteractions", cpdb_version)

# Download the CellPhoneDB database---------------------------------------------
from cellphonedb.utils import db_utils

db_utils.download_database(cpdb_target_dir, cpdb_version)

# Define inputs and path
pd.set_option("display.max_columns", 100)
os.chdir("03_CellCellInteractions")
out_path = "03_CellCellInteractions/output"
cpdb_file_path = "v5.0.0/cellphonedb.zip"
meta_file_path = "data/tumours_metadata.tsv"
counts_file_path = "data/tumours.h5ad"
degs_file_path = "data/tumours_DEGs.tsv"

# Run CellPhoneDB method 3------------------------------------------------------
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

cpdb_results_m3 = cpdb_degs_analysis_method.call(
         cpdb_file_path = cpdb_file_path,
         meta_file_path = meta_file_path,
         counts_file_path = counts_file_path,
         degs_file_path = degs_file_path,
         counts_data = "hgnc_symbol",
         threshold = 0.1,
         #threads = 4,
         score_interactions = True,
         output_suffix="",
         output_path = out_path)

# Generate plots----------------------------------------------------------------
import ktplotspy as kpy
import matplotlib
%matplotlib inline
plt.rcParams["pdf.fonttype"] = "truetype"
matplotlib.rcParams.update({"font.size": 5})

from IPython.display import HTML, display

adata = ad.read_h5ad("data/tumours.h5ad")

means = pd.read_csv("output/degs_analysis_means_.txt", sep="\t")
pvals = pd.read_csv("output/degs_analysis_relevant_interactions_.txt", sep="\t")
decon = pd.read_csv("output/degs_analysis_deconvoluted_.txt", sep="\t")
scores = pd.read_csv("output/degs_analysis_interaction_scores_.txt", sep="\t")

# Extended Data Fig 6e
p = kpy.plot_cpdb_heatmap(pvals = pvals, degs_analysis=True, figsize=(5, 5), title="Sum of significant interactions", symmetrical=False)
p.savefig("figures/degs_analysis_heatmaps_threshold.pdf")

# Figure 5c and Extended Data Fig. 6f
p=kpy.plot_cpdb(
    adata = adata,
    cell_type1 = "G0_arrested|Fast_cycling",
    cell_type2 = "CD4.T|CD8.T|Treg|Bcell|Plas.B|pDC|cDC|Mac|TAM|Mast|CAF|Fib|Peri|EC|G0_arrested|Fast_cycling",
    means = means,
    pvals = pvals,
    celltype_key = "cell_type",
    highlight_size = 1,
    degs_analysis=True,
    interaction_scores=scores,
    cmap_name="PuOr",
    scale_alpha_by_interaction_scores=True,
    keep_significant_only=True,
    highlight_col = "#000000",
    figsize = (18, 25),
    min_interaction_score=30,
    max_size=5
    )
p.save("figures/TME_to_all_dotplot.pdf")

# Figure 5b
p=kpy.plot_cpdb(
    adata = adata,
    cell_type1 = "G0_arrested|Fast_cycling",
    cell_type2 = "G0_arrested|Fast_cycling",
    means = means,
    pvals = pvals,
    celltype_key = "cell_type",
    highlight_size = 1,
    degs_analysis=True,
    interaction_scores=scores,
    cmap_name="PuOr",
    scale_alpha_by_interaction_scores=True,
    keep_significant_only=True,
    figsize = (10, 15),
    min_interaction_score=20,
)
p.save("figures/G0_arrested_Fast_cycling_dotplot.pdf")
