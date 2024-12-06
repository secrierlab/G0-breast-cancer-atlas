# Prepare the data--------------------------------------------------------------
import warnings
warnings.filterwarnings("ignore")

# setting up R dependencies
import anndata2ri
import rpy2
from rpy2.robjects import r
import random

anndata2ri.activate()

%load_ext rpy2.ipython

# In R, do the following:
# load the packages
# library(Seurat)
# library(zellkonverter)
# library(SingleCellExperiment)
# set the working directory
# projectDir <- "~/BreastCancerG0arrest"
# setwd(paste0(projectDir, "/05_CellPhoneDB"))
# CD8Tcells <- readRDS(paste0(projectDir, "/04_InvestigatingCD8TcellExhaustion/data/NK.CD8Tcells_with_annotations.rds"))
# patients_with_TEX <- c("Patient.6", "Patient.9", "Patient.10", "Patient.26", "Patient.30", "Patient.31", "Patient.36", "Patient.39", "Patient.40", "Patient.41", "Patient.42", "Patient.45","Patient.46", "Patient.5")
# Tex <- subset(CD8Tcells, subset = patient %in% patients_with_TEX)
# CD8Tcells$celltype <- as.character(CD8Tcells$celltype)
# Tex <- colnames(CD8Tcells)[CD8Tcells$celltype == "Exhausted CD8+ T cell"]
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
# integrated$celltype[colnames(integrated) %in% Tex] <- "Exhausted_CD8_T_cell"
# Idents(integrated) <- integrated$celltype

# subset the object to only include cancer patients
# integrated <- subset(integrated, subset = (disease %in% "cancer") & !(integrated$patient %in% c("Patient.12", "Patient.13", "Patient.14", "Patient.15")))
# subset the object to only include invasive breast cancer
# integrated <- subset(integrated, subset = type %in% c("ER", "PR", "HER", "TNBC"))
# subset Tex+ patients
# integrated <- subset(integrated, subset = patient %in% patients_with_TEX)

# integrated <- subset(integrated, subset = celltype != "CD8.T")

# save the metadata
# Idents(integrated) <- integrated$celltype
# metadata <- integrated@meta.data
# metadata$barcode_sample <- rownames(metadata)
# metadata <- metadata[, c("barcode_sample", "celltype")]
# colnames(metadata) <- c("barcode_sample", "cell_type")
# write.table(metadata, file = "data/cd8_metadata.tsv", sep = "\t", quote = F, row.names = F)

# prepare the data for CellPhoneDB
# integrated$cell_type <- Idents(integrated)

# find DEGs for each cell type (for only method 3 of cellphonedb)
# DEGs <- FindAllMarkers(integrated, only.pos = T, min.pct = 0.1)
# DEGs_reordered <- DEGs[, c(6, 7, 5, 1, 2, 3, 4)]
# colnames(DEGs_reordered) <- c("cell_type", "gene", "p_val_adj", "p_val", "avg_log2FC", "pct.1", "pct.2")
# write.table(DEGs_reordered, file = "data/cd8_DEGs.tsv", row.names = F, sep = "\t", quote = F)

# save the anndata object
# sce <- as.SingleCellExperiment(integrated)
# writeH5AD(sce, file = "data/cd8.h5ad")

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
cpdb_target_dir = os.path.join("/Users/cenkcelik/BreastCancerG0arrest/05_CellPhoneDB", cpdb_version)

# Define paths------------------------------------------------------------------
pd.set_option("display.max_columns", 100)
os.chdir("/Users/cenkcelik/BreastCancerG0arrest/05_CellPhoneDB")
out_path = "/Users/cenkcelik/BreastCancerG0arrest/05_CellPhoneDB/output"
cpdb_file_path = "v5.0.0/cellphonedb.zip"
meta_file_path = "data/cd8_metadata.tsv"
counts_file_path = "data/cd8.h5ad"
degs_file_path = "data/cd8_DEGs.tsv"

# Run CellPhoneDB method 3------------------------------------------------------
from cellphonedb.src.core.methods import cpdb_degs_analysis_method

cpdb_results_m3 = cpdb_degs_analysis_method.call(
         cpdb_file_path = cpdb_file_path,
         meta_file_path = meta_file_path,
         counts_file_path = counts_file_path,
         degs_file_path = degs_file_path,
         counts_data = "hgnc_symbol",
         threshold = 0.1,
         threads = 8,
         score_interactions = True,
         output_suffix="cd8",
         output_path = out_path)

import ktplotspy as kpy
import matplotlib
%matplotlib inline
plt.rcParams["pdf.fonttype"] = "truetype"
matplotlib.rcParams.update({"font.size": 5})

adata = ad.read_h5ad("data/cd8.h5ad")
means = pd.read_csv("output/degs_analysis_means_cd8.txt", sep="\t")
pvals = pd.read_csv("output/degs_analysis_relevant_interactions_cd8.txt", sep="\t")
decon = pd.read_csv("output/degs_analysis_deconvoluted_cd8.txt", sep="\t")
scores = pd.read_csv("output/degs_analysis_interaction_scores_cd8.txt", sep="\t")

# Extended Data Fig 7e
p=kpy.plot_cpdb(
    adata = adata,
    cell_type1 = "CD4.T||Treg|Bcell|Plas.B|pDC|cDC|Mac|TAM|Mast|CAF|Fib|Peri|EC|G0_arrested|Fast_cycling|Slow_cycling|Exhausted_CD8_T_cell",
    cell_type2 = "Exhausted_CD8_T_cell",
    means = means,
    pvals = pvals,
    celltype_key = "cell_type",
    highlight_size = 1,
    degs_analysis=True,
    interaction_scores=scores,
    cmap_name="PuOr",
    scale_alpha_by_interaction_scores=True,
    keep_significant_only=True,
    figsize = (15, 25),
    max_size=5,
    min_interaction_score=30,
)
p.save("figures/TME_to_Tex_dotplot.pdf")
