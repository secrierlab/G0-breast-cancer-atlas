import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams["pdf.fonttype"] = "truetype"

import warnings
warnings.simplefilter("ignore")

import scanpy as sc
import drug2cell as d2c
import blitzgsea as blitz
import numpy as np
import pandas as pd
import anndata as ad

import scipy.sparse as sp
from tqdm import tqdm

sc.set_figure_params(dpi=80, scanpy=True)
sc.settings._vector_friendly = False
sc.logging.print_header()

from IPython.display import HTML, display

# read the data
adata = sc.read_h5ad("/Users/cenkcelik/BreastCancerG0arrest/09_SpatialAnalysis/newly_deconvolved_data_with_hotspots1510.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()

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
    "T_cells_c9_NK_cells_AREG",
    "G0_arrested",
    "Fast_cycling"
]

quiescence_labels = ["G0_arrested", "Fast_cycling"]

cell_type_labels_hotspots = [label + '_hot' for label in cell_type_labels]

adata.obs['cell_type'] = adata.obs[cell_type_labels_hotspots].idxmax(axis=1)
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

# if tumour_cells == 1 and (cell_type != G0_arrested or Fast_cycling), then it is a Slow_cycling cell
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
adata.obs["cell_type"] = adata.obs["cell_type"].cat.add_categories('Slow_cycling_hot')
adata.obs.loc[(adata.obs["tumour_cells"] == 1) & (~adata.obs["cell_type"].isin(["G0_arrested_hot", "Fast_cycling_hot"])), "cell_type"] = "Slow_cycling_hot"
# Subset tumour_cells
adata = adata[adata.obs["tumour_cells"] == 1].copy()

def select_slide(adata, s, batch_key="sample"):
    r"""This function selects the data for one slide from the spatial anndata object.

    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param batch_key: column in adata.obs listing experiment name for each location
    """

    slide = adata[adata.obs[batch_key].isin([s]), :].copy()
    s_keys = list(slide.uns["spatial"].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]

    slide.uns["spatial"] = {s_spatial: slide.uns["spatial"][s_spatial]}

    return slide

# Placeholder lists for AnnData objects and their spatial data
score_results = []
drug2cell_results = []

# Iterate over each slide batch with a progress bar
for slide in tqdm(adata.obs["batch"].unique(), desc="Processing slides", unit="slide"):
    # Subset the AnnData object for this slide
    adata_slide = adata[adata.obs["batch"] == slide]
    print(f"Processing slide {slide} (patient {int(slide) + 1}).")

    # Apply the d2c.score() function
    d2c.score(adata_slide, layer="log1p")

    # Run differential drug scores
    sc.tl.rank_genes_groups(adata_slide.uns["drug2cell"], method="t-test_overestim_var", groupby="cell_type")
    sc.tl.filter_rank_genes_groups(adata_slide.uns["drug2cell"], min_in_group_fraction=0.5, max_out_group_fraction=0.5, min_fold_change=2)
    sc.pl.rank_genes_groups_dotplot(adata_slide.uns["drug2cell"], n_genes=5, swap_axes=True, dendrogram=False, cmap="RdYlBu_r")
    # Add drug scores to adata.obs
    drug_data = adata_slide.uns["drug2cell"].copy()
    drug_data.obs[drug_data.var_names] = drug_data.X.tolist()  # Store drug scores
    drug2cell_results.append(drug_data)

    # Append the processed slide to score_results
    score_results.append(adata_slide)

# Merge all the scored AnnData objects
merged_adata = score_results[0].concatenate(*score_results[1:], batch_key="batch")
merged_drug2cell = drug2cell_results[0].concatenate(*drug2cell_results[1:], batch_key="batch")

# Add merged drug scores and spatial info to merged AnnData objects
merged_adata.uns["drug2cell"] = merged_drug2cell
merged_adata.uns["drug2cell"].uns["spatial"] = adata.uns["spatial"]
merged_adata.uns["spatial"] = adata.uns["spatial"]

# Display final message
print("Finished processing and merging AnnData object.")

print("Adding drug scores to adata.obs...")
merged_adata.obs[merged_adata.uns["drug2cell"].var_names] = merged_adata.uns["drug2cell"].X.tolist()
print("Finished adding drug scores to adata.obs.")
adata = merged_adata
del merged_adata

adata.write_h5ad("adata_with_drug_scores.h5ad")