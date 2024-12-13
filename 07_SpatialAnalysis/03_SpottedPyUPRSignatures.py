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

# Define the path to the spottedpy package
# sys.path.append(r"/Users/cenkcelik/BreastCancerG0arrest/10_SpatialAnalysis/SpottedPy/spottedPy")
# Import the package
import spottedpy as sp

# Set the working directory
os.chdir(r"/Users/cenkcelik/BreastCancerG0arrest/10_SpatialAnalysis")

# Load the anndata object
adata = sc.read_h5ad("newly_deconvolved_data_with_hotspots1410.h5ad")

# read in the UPR signature xlsx file in sheet "dense"
upr_1 = pd.read_excel("/Users/cenkcelik/BreastCancerG0arrest/05_GeneRegulatoryNetworks/database/Human-Proteostasis-Network-2.0-2024-0415-1-e1e497a1e8cd5a06.xlsx", sheet_name="dense")
upr_2 = pd.read_excel("/Users/cenkcelik/BreastCancerG0arrest/05_GeneRegulatoryNetworks/database/ER_Proteostasis_250624.xlsx", sheet_name="Sheet1")

# convert to a dict with Branch as keys and all gene symbols as values
upr_1_dict = {branch: upr_1[upr_1["Branch"] == branch]["Gene Symbol"].values for branch in upr_1["Branch"].unique()}
upr_2_dict = {branch: upr_2[upr_2["Group"] == branch]["Gene Symbol"].values for branch in upr_2["Group"].unique()}

# combine the two dicts
upr_dict = {**upr_1_dict, **upr_2_dict}

# Signature labels
signature_labels = upr_dict.keys()
quiescence_labels = ["G0_arrested","Fast_cycling"]

# Check if genes are in the adata object, remove those that are not
for branch, genes in upr_dict.items():
    upr_dict[branch] = [gene for gene in genes if gene in adata.var_names]

signatures = []

for key, value in upr_dict.items():
    # Replace spaces with underscores to create valid Python variable names
    key = key.replace(" ", "_")
    # Add the variable (not its name as a string) to the tuple
    signatures.append((value, key))

# Process each batch and calculate scores
for batch in adata.obs["batch"].unique():
    batch_subset = adata[adata.obs["batch"] == batch]
    
    for (list, branch) in signatures:
        # Calculate scores for the hallmark genes
        sc.tl.score_genes(batch_subset, list, score_name=branch)
        
        # Update the scores in the main dataframe
        adata.obs.loc[adata.obs["batch"] == batch, branch] = batch_subset.obs[branch].values

# Define a function to scale scores
def scale_scores(adata, score_names, batch_column="batch"):
    scaler = MinMaxScaler()
    
    for batch in adata.obs[batch_column].unique():
        adata_subset = adata[adata.obs[batch_column] == batch]
        
        # Scale scores for the given score names
        adata_subset.obs[score_names] = scaler.fit_transform(adata_subset.obs[score_names])
        
        # Update the main adata object with scaled scores
        adata.obs.loc[adata.obs[batch_column] == batch, score_names] = adata_subset.obs[score_names].values

# Apply scaling to signatures
labels = [t[1] for t in signatures]
scale_scores(adata, labels)

# Hotspot analysis--------------------------------------------------------------
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

adata = calculate_hotspots(adata, labels, filter_columns="tumour_cells", filter_value=1)

# Generate signature_labels_hotspots by appending "_hot" to each label in signature_labels
signature_labels_hotspots = [label + "_hot" for label in labels]

# Define primary and comparison variables
primary_variables = ["G0_arrested_hot", "Fast_cycling_hot"]
comparison_variables = signature_labels_hotspots

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
                       comparison_vars=signature_labels_hotspots,
                       bubble_size=(300, 300),
                       fig_size=(15, 6),
                       save_path="figures/signature_upr_custom_scatter.pdf",
                       compare_distribution_metric="median")

# Generate a list of tuples with column names and save paths using drug_hotspots
plots = [(label, f"figures/{label}.pdf") for label in signature_labels_hotspots]

for slide in adata.obs["batch"].unique():
    # create folder for each slide
    os.makedirs(f"figures/show_{slide}_figures", exist_ok=True)

# Loop through the list and plot each hotspot
for column_name, save_path in plots:
    sp.plot_hotspots(adata, column_name=column_name, save_path=save_path)

# differences per batch (i.e. slide)
data=sp.plot_bubble_chart_by_batch(distances,
                                   primary_variable_value="G0_arrested_hot",
                                   comparison_variable_values=signature_labels_hotspots,
                                   reference_variable="Fast_cycling_hot",
                                   save_path="figures/upr_bubble_chart_by_batch.pdf",
                                   pval_cutoff=0.05,
                                   fig_size=(12,12))