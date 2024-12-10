import warnings
warnings.filterwarnings("ignore")

import os, pickle, requests, glob

import operator as op
from cytoolz import compose

from collections import OrderedDict

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import loompy as lp
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams["pdf.fonttype"] = "truetype"
sc.set_figure_params(scanpy=True, dpi=80, dpi_save=300, frameon=False, vector_friendly=False, fontsize=14, figsize=(6,6), color_map=None, format="pdf", facecolor=None, transparent=True, ipython_format="png2x")
sc.settings._vector_friendly = False
sc.settings.verbosity = 0
sc.logging.print_header()

from pyscenic.export import add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss

import scipy.cluster.hierarchy as spc
from IPython.display import HTML, display

# set project directory:
project_dir = "/Users/cenkcelik/BreastCancerG0arrest/05_GeneRegulatoryNetworks/"
os.chdir(project_dir)

# some helper functions
def derive_regulons(motifs, db_names=("hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings", 
                                      "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores", 
                                      "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings", 
                                      "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores")):
    
    motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method "weight>50.0%" (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains("weight>50.0%")), motifs.Context), dtype=np.bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=np.bool) & \
        np.fromiter(map(contains("activating"), motifs.Context), dtype=np.bool)]

    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter(lambda r: len(r) >= 10, df2regulons(motifs[(motifs["NES"] >= 3.0) 
                                                                      & ((motifs["Annotation"] == "gene is directly annotated")
                                                                        | (motifs["Annotation"].str.startswith("gene is orthologous to")
                                                                           & motifs["Annotation"].str.endswith("which is directly annotated for motif")))
                                                                     ])))
    
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))

def download_file(url, destination_directory):
    # Ensure the destination directory exists
    os.makedirs(destination_directory, exist_ok=True)

    # Get the file name from the URL
    file_name = os.path.join(destination_directory, url.split("/")[-1])

    # Download the file
    response = requests.get(url)
    with open(file_name, "wb") as file:
        file.write(response.content)

    print(f"File downloaded to {file_name}")

def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ["k"] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment="center", verticalalignment="center")
    return f

def savepdf(fname: str, fig, folder: str="figures") -> None:
    """
    Save figure as PDF image format.
    """
    fig.tight_layout()
    fig.savefig(os.path.join(folder, fname), format="pdf")

# Prep data for SCENIC--------------------------------------------
# load the data
adata = sc.read_h5ad("resources/adata.h5ad")
adata.var_names_make_unique()
# normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
# log transform the data
sc.pp.log1p(adata)
adata.raw = adata

# subset invasive breast carcinoma
adata = adata[adata.obs["type"].isin(["ER", "PR", "HER", "TNBC"])]
# Embed PCs precalculated in Seurat
adata.obsm["X_pca"] = np.vstack((adata.obs["PC_1"].to_numpy(), adata.obs["PC_2"].to_numpy())).T
# Embed UMAP precalculated in Seurat
adata.obsm["X_umap"] = np.vstack((adata.obs["UMAP_1"].to_numpy(), adata.obs["UMAP_2"].to_numpy())).T

# select highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

# subset G0 arrested cells
adata = adata[adata.obs.celltype == "Slow_cycling", :]
adata.shape

# Downloaded TFs from pySCENIC github repo: https://github.com/aertslab/pySCENIC/tree/master/resources
url = "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"

tfs_path = "resources/allTFs_hg38.txt"
if not os.path.exists(tfs_path):
    download_file(url, "resources")

# loom directory
loom_path = os.path.join("resources", "Slow_cycling_processed_input.loom")
loom_path_output = os.path.join("resources", "Slow_cycling_processed_output.loom")
tfs = [tf.strip() for tf in open(tfs_path)]

# as a general QC. We inspect that our object has transcription factors listed in our main annotations.
print(f"{np.sum(adata.var.index.isin(tfs))} out of {len(tfs)} TFs are found in the object")

# use highly variable genes
use_hvg = True
if use_hvg:
    mask = (adata.var["highly_variable"] == True) | adata.var.index.isin(tfs)
    adata = adata[:, mask]

# save the loom file
row_attributes = {"Gene": np.array(adata.var.index), }
col_attributes = {
    "CellID": np.array(adata.obs.index),
    "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
}

lp.create(loom_path, adata.X.transpose(), row_attributes, col_attributes)

# run pySCENIC---------------------------------------------------
num_workers = 8

outpath_adj = os.path.join(project_dir, "results/Slow_cycling_adjacencies.csv")
if not os.path.exists(outpath_adj):
    !pyscenic grn {loom_path} {tfs_path} -o $outpath_adj --num_workers {num_workers} --seed 42 --method grnboost2

# see the top TF-target gene associations
adjacencies = pd.read_csv(outpath_adj, index_col=False, sep=",")
print(f"Number of associations: {adjacencies.shape[0]}")
adjacencies.head()

# visualise the distribution of weights
plt.hist(np.log10(adjacencies["importance"]), bins=50)
plt.xlim([-10, 5])

# load TSS annotations
db_glob = "auxiliaries/*feather"
db_names = " ".join(glob.glob(db_glob))

# load the catalog of motif-to-TF associations
motif_path = "auxiliaries/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

# run the context-specific co-expression network inference--------
outpath_motifs = "results/Slow_cycling_motifs.csv"

if not os.path.exists(outpath_motifs):
    !pyscenic ctx $outpath_adj $db_names --annotations_fname "auxiliaries/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" --expression_mtx_fname {loom_path} --output {outpath_motifs} --mask_dropouts --num_workers {num_workers}

# explore the regulons
n_genes_detected_per_cell = np.sum(adata.X > 0, axis=1)
percentiles = pd.Series(n_genes_detected_per_cell.flatten().A.flatten()).quantile([0.01, 0.05, 0.10, 0.50, 1])
print(percentiles)

# plot the distribution of # of genes detected per cell to determine the threshold
fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=100)
sns.distplot(n_genes_detected_per_cell, norm_hist=False, kde=False, bins="fd")
for i, x in enumerate(percentiles):
    fig.gca().axvline(x=x, ymin=0, ymax=1, color="red")
    ax.text(
        x=x,
        y=ax.get_ylim()[1],
        s=f"{int(x)} ({percentiles.index.values[i]*100}%)",
        color="red",
        rotation=30,
        size="x-small",
        rotation_mode="anchor",
    )
ax.set_xlabel("# of genes")
ax.set_ylabel("# of cells")
fig.tight_layout()

# run AUCell-----------------------------------------------------
if not os.path.exists(loom_path_output):
    !pyscenic aucell $loom_path "results/Slow_cycling_motifs.csv" --output {loom_path_output} --num_workers {num_workers}

# collect SCENIC AUCell output
lf = lp.connect(loom_path_output, mode="r+", validate=False)
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

# add metadata to the adata
ad_auc_mtx = ad.AnnData(auc_mtx)
sc.pp.neighbors(ad_auc_mtx)
sc.tl.umap(ad_auc_mtx)
sc.tl.leiden(ad_auc_mtx)

# add to the adata
adata.obsm["X_umap_aucell"] = ad_auc_mtx.obsm["X_umap"]

# visualise the new embedding space
sc.pl.embedding(adata, basis="X_umap_aucell", color="type")

auc_mtx["type"] = adata.obs["type"]
mean_auc_by_type = auc_mtx.groupby("type").mean()
del auc_mtx["type"]

# top regulons
top_n = 100
top_tfs = mean_auc_by_type.max(axis=0).sort_values(ascending=False).head(top_n).index.tolist()

# subset the regulons
auc_mtx = auc_mtx[top_tfs]

# aucell correlation matrix
correlation_matrix = auc_mtx.corr(method="spearman")

# fix infinities
correlation_matrix = np.nan_to_num(correlation_matrix)

pdist = spc.distance.pdist(correlation_matrix)
linkage = spc.linkage(pdist, method="complete")
idx = spc.fcluster(linkage, 0.5 * pdist.max(), "distance")

spc.dendrogram(linkage, labels=auc_mtx.columns, orientation="top", leaf_rotation=90)

threshold = 0.2
labels = spc.fcluster(linkage, threshold, criterion="distance")

# Keep the indices to sort labels
labels_order = np.argsort(labels)

# Build a new dataframe with the sorted columns
for idx, i in enumerate(auc_mtx.columns[labels_order]):
    if idx == 0:
        clustered = pd.DataFrame(auc_mtx[i])
    else:
        df_to_append = pd.DataFrame(auc_mtx[i])
        clustered = pd.concat([clustered, df_to_append], axis=1)

correlations_clustered = clustered.corr(method="spearman")

# plot the heatmap (Extended Data Figure 5d):
plt.figure(figsize=(25, 25), dpi=150)
sns.set_theme(font_scale=1)
g=sns.heatmap(correlations_clustered, cmap="RdYlBu_r", square=True, annot=False, center=0, cbar_kws={"shrink": 0.2})
plt.savefig("figures/Slow_cycling_correlation_heatmap.pdf", dpi=72, bbox_inches='tight')
plt.show()

top_tfs = correlations_clustered.columns.tolist()

# average AUC per regulon
mean_auc_by_type_top_n = mean_auc_by_type[[c for c in mean_auc_by_type.columns if c in top_tfs]]
mean_auc_by_type_top_n = mean_auc_by_type_top_n[top_tfs]

# Regulon enrichment
# load the motifs
df_motifs = load_motifs("results/Slow_cycling_motifs.csv")
# derive regulons
regulons = derive_regulons(df_motifs)

with open("results/Slow_cycling_motifs.dat", "wb") as f: 
    pickle.dump(regulons, f)

# rename regulons
auc_mtx.columns = auc_mtx.columns.str.replace("\(\+\)", "")
# add metadata to AnnData
add_scenic_metadata(adata, auc_mtx, regulons)

TYPE_COLORS = ["#5FB8E9","#439788","#BC4024","#DB4677"]

# calculate regulon specificity scores
rss_score = regulon_specificity_scores(auc_mtx, adata.obs.type)

# Function to convert gene names
def convert_gene_name(gene_name):
    return f"Regulon({gene_name})"

# extract top regulons without (+)
tf_names = [gene.rstrip("(+)") for gene in top_tfs]
adata_top_tfs = adata[:, adata.var_names.isin(tf_names)]
# Create a new list with converted gene names
regulon_names = [convert_gene_name(gene) for gene in tf_names]

# write regulons MYC, FOS, MYBL2, FOXM1 AUCell scores to a csv file
outpath = os.path.join("results", "Slow_cycling_regulons_aucell.csv")
auc_mtx[["MYC", "FOS", "MYBL2", "FOXM1"]].to_csv(outpath)

# plot them on the new embedding space (Extended Data Figure 4e):
sc.pl.embedding(adata, basis="X_umap_aucell", 
                color=["Regulon(FOS)", "Regulon(FOXM1)", "Regulon(MYC)", "Regulon(MYBL2)", "Regulon(NR2F2)", "Regulon(SOX9)", "type"], 
                title=["Regulon(FOS)", "Regulon(FOXM1)", "Regulon(MYC)", "Regulon(MYBL2)", "Regulon(NR2F2)", "Regulon(SOX9)", "type"], 
                palette=TYPE_COLORS, cmap="Reds", 
                save="_Slow_cycling_rss_regulon_activity.pdf")