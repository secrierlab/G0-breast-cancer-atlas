import warnings
warnings.filterwarnings("ignore")

import os
import sys
import scanpy as sc
# Figure settings
sc.settings.set_figure_params(dpi = 80, facecolor="white", frameon=False, scanpy=True, vector_friendly=True)
sc.settings._vector_friendly = False

import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import zscore
from matplotlib import pyplot as plt

# Define the path to the spottedpy package
sys.path.append(r"/Users/cenkcelik/BreastCancerG0arrest/07_SpatialAnalysis/SpottedPy/spottedPy")
# Import the package
import spottedpy as sp

# Set the working directory
os.chdir(r"/Users/cenkcelik/BreastCancerG0arrest/07_SpatialAnalysis")

# Load the anndata object
adata = sc.read_h5ad("newly_deconvolved_data_with_hotspots3107.h5ad")

# Subset adata for cells where tumour_cells == 1
tumour_adata = adata[adata.obs["tumour_cells"] == 1].copy()

# Quiescence signature genes from Wiecek et al. 2023 Genome Biology
upregulated_genes = ["CFLAR", "CALCOCO1", "YPEL3", "CST3", "SERINC1", "CLIP4", "PCYOX1", "TMEM59","RGS2", "YPEL5", "CD63", "KIAA1109", "CDH13", "GSN", "MR1", "CYB5R1", "AZGP1","ZFYVE1", "DMXL1", "EPS8L2", "PTTG1IP", "MIR22HG", "PSAP", "GOLGA8B", "NEAT1", "TXNIP", "MTRNR2L12"]

downregulated_genes = ["NCAPD2", "PTBP1", "MPHOSPH9", "NUCKS1", "TCOF1", "SART3", "SNRPA", "KIF22", "HSP90AA1", "WBP11", "CAD", "SF3B2", "KHSRP", "WDR76", "NUP188", "HSP90AB1", "HNRNPM", "SMARCB1", "PNN", "RBBP7", "NPRL3", "USP10", "SGTA", "MRPL4",
                       "PSMD3", "KPNB1", "CBX1", "LRRC59", "TMEM97", "NSD2", "PRPF19", "PTGES3", "CPSF6", "SRSF3", "TCERG1", "SMC4", "EIF4G1", "ZNF142", "MSH6", "MRPL37", "SFPQ", "STMN1", "ARID1A", "PROSER1", "DDX39A", "EXOSC9", "USP22", "DEK",
                       "DUT", "ILF3", "DNMT1", "NASP", "HMGB1P5", "SRRM1", "GNL2", "RNF138", "SRSF1", "TRA2B", "SMPD4", "ANP32B", "HMGA1", "MDC1", "HADH", "ARHGDIA", "PRCC", "HDGF", "SF3B4", "UBAP2L", "ILF2", "PARP1", "LBR", "CNOT9", "PPRC1", "SSRP1", 
                       "CCT5", "DLAT", "HNRNPU", "LARP1", "SCAF4", "RRP1B", "RRP1", "CHCHD4", "GMPS", "RFC4", "SLBP", "PSIP1", "HNRNPK", "SKA3", "DIS3L", "USP39", "GPS1", "PA2G4", "HCFC1", "SLC19A1", "ETV4", "RAD23A", "DCTPP1", "RCC1", "EWSR1", "ALYREF", 
                       "PTMA", "HMGB1", "POM121", "MCMBP", "TEAD4", "CHAMP1", "TOP1", "PRRC2A", "RBM14", "HMGB1P6", "POM121C", "UHRF1"]

# Filter the genes for those that are in tumour_adata.var_names
upregulated_genes = [gene for gene in upregulated_genes if gene in tumour_adata.var_names]
downregulated_genes = [gene for gene in downregulated_genes if gene in tumour_adata.var_names]

# Precompute the z-scores for the tumour cells subset
gene_expr = pd.DataFrame(tumour_adata.X.toarray(), index=tumour_adata.obs_names, columns=tumour_adata.var_names).T
gene_expr_score = zscore(gene_expr, axis=1)

# Calculate the scores for tumour cells
pos_score = gene_expr_score.loc[upregulated_genes].sum(axis=0) / np.sqrt(len(upregulated_genes))
neg_score = gene_expr_score.loc[downregulated_genes].sum(axis=0) / np.sqrt(len(downregulated_genes))

# Compute the z-score differences for tumour cells
zscore_differences = pos_score - neg_score

# Assign the quiescence scores back to the original adata.obs
adata.obs.loc[adata.obs["tumour_cells"] == 1, "quiescence_score"] = zscore_differences

q25, q75 = adata.obs["quiescence_score"].quantile([0.25, 0.75])

# Create "G0_arrested" column: assign values where "quiescence" > Q75, otherwise 0
adata.obs["G0_arrested"] = adata.obs["quiescence_score"].where(adata.obs["quiescence_score"] > q75, 0)

# Create "Fast_cycling" column: assign negative values where "quiescence" < Q25, otherwise 0
adata.obs["Fast_cycling"] = adata.obs["quiescence_score"].where(adata.obs["quiescence_score"] < q25, 0) * -1

# Create "Intermediate" column: assign values where "quiescence" is between Q25 and Q75, otherwise 0
adata.obs["Slow_cycling"] = adata.obs["quiescence_score"].where((adata.obs["quiescence_score"] >= q25) & (adata.obs["quiescence_score"] <= q75), 0)

# Cell type labels
cell_type_labels = [
    "B cells Memory",
    "B cells Naive",
    "CAFs MSC iCAF-like s1",
    "CAFs MSC iCAF-like s2",
    "CAFs Transitioning s3",
    "CAFs myCAF like s4",
    "CAFs myCAF like s5",
    "Cycling PVL",
    "Cycling_Myeloid",
    "Endothelial ACKR1",
    "Endothelial CXCL12",
    "Endothelial Lymphatic LYVE1",
    "Endothelial RGS5",
    "Myeloid_c0_DC_LAMP3",
    "Myeloid_c10_Macrophage_1_EGR1",
    "Myeloid_c11_cDC2_CD1C",
    "Myeloid_c12_Monocyte_1_IL1B",
    "Myeloid_c1_LAM1_FABP5",
    "Myeloid_c2_LAM2_APOE",
    "Myeloid_c3_cDC1_CLEC9A",
    "Myeloid_c4_DCs_pDC_IRF7",
    "Myeloid_c5_Macrophage_3_SIGLEC1",
    "Myeloid_c7_Monocyte_3_FCGR3A",
    "Myeloid_c8_Monocyte_2_S100A9",
    "Myeloid_c9_Macrophage_2_CXCL10",
    "PVL Differentiated s3",
    "PVL Immature s1",
    "PVL_Immature s2",
    "Plasmablasts",
    "T_cells_c0_CD4+_CCR7",
    "T_cells_c10_NKT_cells_FCGR3A",
    "T_cells_c11_MKI67",
    "T_cells_c1_CD4+_IL7R",
    "T_cells_c2_CD4+_T-regs_FOXP3",
    "T_cells_c3_CD4+_Tfh_CXCL13",
    "T_cells_c4_CD8+_ZFP36",
    "T_cells_c5_CD8+_GZMK",
    "T_cells_c6_IFIT1",
    "T_cells_c7_CD8+_IFNG",
    "T_cells_c8_CD8+_LAG3",
    "T_cells_c9_NK_cells_AREG"
]

quiescence_labels = ["G0_arrested", "Fast_cycling", "Slow_cycling"]

labels = cell_type_labels + ["Slow_cycling"]

# Define a function to scale scores
def scale_scores(adata, score_names, batch_column="batch"):
    scaler = MinMaxScaler()
    
    for batch in adata.obs[batch_column].unique():
        adata_subset = adata[adata.obs[batch_column] == batch]
        
        # Scale scores for the given score names
        adata_subset.obs[score_names] = scaler.fit_transform(adata_subset.obs[score_names])
        
        # Update the main adata object with scaled scores
        adata.obs.loc[adata.obs[batch_column] == batch, score_names] = adata_subset.obs[score_names].values

# Apply scaling to TME labels
scale_scores(adata, labels)

## Hotspot analysis--------------------------------------------------------------------------------
def calculate_hotspots(adata, labels, filter_columns=None, filter_value=None):
    """
    Calculate hotspots for a list of labels.
    
    Parameters:
        adata: AnnData object
        labels: List of column names for which to calculate hotspots
        filter_columns: Column name to filter on (optional)
        filter_value: Value to filter for in filter_columns (optional)
    """
    for label in tqdm(labels):
        adata = sp.create_hotspots(
            adata,
            column_name=label,
            neighbours_parameters=12,
            relative_to_batch=True,
            number_hotspots=True,
            filter_columns=filter_columns,
            filter_value=filter_value
        )
    return adata

# Calculate hotspots for quiescence_labels with tumour cells only
adata = calculate_hotspots(adata, quiescence_labels, filter_columns="tumour_cells", filter_value=1)

# Calculate hotspots for cell_type_labels
adata = calculate_hotspots(adata, cell_type_labels)

# Generate cell_type_labels_hotspots by appending "_hot" to each label in cell_type_labels
cell_type_labels_hotspots = [label + "_hot" for label in cell_type_labels]

## G0 vs Fast cycling--------------------------------------------------------------------------------
# Define primary and comparison variables
primary_variables = ["G0_arrested_hot", "Fast_cycling_hot"]
comparison_variables = cell_type_labels_hotspots

# Calculate distances
distances = sp.calculateDistances(
    adata,
    primary_variables=primary_variables,
    comparison_variables=comparison_variables,
    empty_hotspot_default_to_max_distance=False,
    split_by_slide_in_batch=True,
    hotspot_number=True
)

sp.plot_custom_scatter(data=distances,
                       primary_vars=["G0_arrested_hot", "Fast_cycling_hot"],
                       comparison_vars=cell_type_labels_hotspots,
                       bubble_size=(300, 300),
                       fig_size=(20, 6),
                       file_save="figures/cell_type_main_custom_scatter_",
                       compare_distribution_metric="median")

# Define a list of tuples with column names and save paths
plots = [
    ("G0_arrested_hot", "G0_arrested_hotspots.pdf"),
    ("Fast_cycling_hot", "Fast_cycling_hotspots.pdf"),
    ("Myeloid_c2_LAM2_APOE_hot", "Myeloid_c2_LAM2_APOE_hotspots.pdf"),
    ("CAFs myCAF like s5_hot", "CAFs_myCAF_like_s5_hotspots.pdf"),
    ("PVL_Immature s2_hot", "PVL_Immature_s2_hotspots.pdf"),
    ("tumour_cells", "tumour_cell_spots.pdf"),
]

# Loop through the list and plot each hotspot
for column_name, save_path in plots:
    sp.plot_hotspots(
        adata,
        column_name=column_name,
        save_path=save_path
    )

# differences per batch (i.e. slide)
# returns distances for further analysis if needed
data=sp.plot_bubble_chart_by_batch(distances, primary_variable_value="G0_arrested_hot",
                               comparison_variable_values=cell_type_labels_hotspots,
                               reference_variable="Fast_cycling_hot",
                               save_path="figures/cell_type_bubble_chart_G0_vs_Fast_by_batch.pdf",
                               pval_cutoff=0.05,
                               fig_size=(16,10))