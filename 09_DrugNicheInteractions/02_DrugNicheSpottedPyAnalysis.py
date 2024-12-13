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

# Import the package
import spottedpy as sp

# Set the working directory
os.chdir(r"/Users/cenkcelik/BreastCancerG0arrest/09_Drug2Cell")

# Load the anndata object
adata = sc.read_h5ad("adata_with_drug_scores.h5ad")

# Convert the 'drug_names' column into a list
drug_names = [
    "CHEMBL2095222|OCRIPLASMIN",
    "CHEMBL2105662|LOMITAPIDE MESYLATE",
    "CHEMBL1200475|DAUNORUBICIN CITRATE",
    "CHEMBL1703|METFORMIN HYDROCHLORIDE",
    "CHEMBL4298211|PEGCETACOPLAN",
    "CHEMBL1201505|FIBRINOLYSIN, HUMAN",
    "CHEMBL1201733|PAZOPANIB HYDROCHLORIDE",
    "CHEMBL1086440|TRICLABENDAZOLE",
    "CHEMBL1200645|ETOPOSIDE PHOSPHATE"
]

# quiescence labels
quiescence_labels = ["G0_arrested","Fast_cycling"]

# Define a function to scale scores
def scale_scores(adata, drug_names, batch_column="batch"):
    scaler = MinMaxScaler()
    
    for batch in adata.obs[batch_column].unique():
        adata_subset = adata[adata.obs[batch_column] == batch]
        
        # Scale scores for the given score names
        adata_subset.obs[drug_names] = scaler.fit_transform(adata_subset.obs[drug_names])
        
        # Update the main adata object with scaled scores
        adata.obs.loc[adata.obs[batch_column] == batch, drug_names] = adata_subset.obs[drug_names].values

# Apply scaling to signatures
scale_scores(adata, drug_names)

## Hotspot analysis--------------------------------------------------------------
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
            neighbours_parameters=10,
            relative_to_batch=True,
            number_hotspots=True,
            filter_columns=filter_columns,
            filter_value=filter_value
        )
    return adata

adata = calculate_hotspots(adata, drug_names)

# Generate signature_labels_hotspots by appending "_hot" to each label in signature_labels
drug_hotspots = [label + "_hot" for label in drug_names]

# Define primary and comparison variables
primary_variables = ["G0_arrested_hot", "Fast_cycling_hot"]
comparison_variables = drug_hotspots

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
                       comparison_vars=drug_hotspots,
                       bubble_size=(300, 300),
                       fig_size=(6, 6),
                       save_path="figures/drug_main_custom_scatter.pdf",
                       compare_distribution_metric="median")

# save bar plot distances per drug 
for drug in drug_hotspots:
    sp.plot_bar_plot_distance(distances,
                              primary_variables=primary_variables,
                              comparison_variables=[drug],
                              fig_size=(3, 3),
                              save_path=f"figures/{drug}_custom_barplot.pdf")
    
# Generate a list of tuples with column names and save paths using drug_hotspots
plots = [(label, f"figures/{label}.pdf") for label in drug_hotspots]

for slide in adata.obs["batch"].unique():
    # create folder for each slide
    os.makedirs(f"figures/show_{slide}_figures", exist_ok=True)

# Loop through the list and plot each hotspot
for column_name, save_path in plots:
    sp.plot_hotspots(adata, column_name=column_name, save_path=save_path)

# differences per batch (i.e. slide)
data=sp.plot_bubble_chart_by_batch(distances,
                                   primary_variable_value="G0_arrested_hot",
                                   comparison_variable_values=drug_hotspots,
                                   reference_variable="Fast_cycling_hot",
                                   save_path="figures/drug_bubble_chart_by_batch.pdf",
                                   pval_cutoff=0.05,
                                   fig_size=(6,12))