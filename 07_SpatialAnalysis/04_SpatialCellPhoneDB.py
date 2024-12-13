import warnings

import scanpy as sc
import liana as li
from liana.method import cellphonedb

from numba.core.errors import NumbaPerformanceWarning
warnings.simplefilter("ignore", NumbaPerformanceWarning)
warnings.filterwarnings('ignore')

adata = sc.read_h5ad("newly_deconvolved_data_with_hotspots1510.h5ad")
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
adata.obs["cell_type"].value_counts()

# Run the cellphonedb function
cellphonedb(
    adata,
    resource_name="cellphonedb",
    layer="log1p",
    min_cells=20,
    expr_prop=0.1,
    de_method="t-test_overestim_var",
    groupby="cell_type",
    seed=42,
    use_raw=False,
    return_all_lrs=False,
    verbose=True,
    n_jobs=6
)

lr_plot = li.pl.dotplot(
    adata=adata,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    source_labels=["Myeloid_c2_LAM2_APOE_hot"],
    target_labels=["G0_arrested_hot", "Fast_cycling_hot"],
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    filter_fun=lambda x: x['cellphone_pvals'] <= 0.001,
    # expression magnitude
    orderby="lr_means",
    orderby_ascending=False,
    top_n=20,
    figure_size=(4, 6),
    cmap="Reds",
    size_range=(1, 5)
)

lr_plot.save('figures/Myeloid_c2_LAM2_APOE_hot_lr_plot.pdf')

lr_plot = li.pl.dotplot(
    adata=adata,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    source_labels=["CAFs myCAF like s5_hot"],
    target_labels=["G0_arrested_hot", "Fast_cycling_hot"],
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    filter_fun=lambda x: x['cellphone_pvals'] <= 0.001,
    # expression magnitude
    orderby="lr_means",
    orderby_ascending=False,
    top_n=20,
    figure_size=(4, 6),
    cmap="Reds",
    size_range=(1, 5)
)

lr_plot.save('figures/CAFs myCAF like s5_hot_lr_plot.pdf')

lr_plot = li.pl.dotplot(
    adata=adata,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    source_labels=["PVL_Immature s2_hot"],
    target_labels=["G0_arrested_hot", "Fast_cycling_hot"],
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    filter_fun=lambda x: x['cellphone_pvals'] <= 0.001,
    # expression magnitude
    orderby="lr_means",
    orderby_ascending=False,
    top_n=20,
    figure_size=(4, 6),
    cmap="Reds",
    size_range=(1, 5)
)

lr_plot.save('figures/PVL_Immature s2_hot_lr_plot.pdf')

lr_plot = li.pl.dotplot(
    adata=adata,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    source_labels=["G0_arrested_hot", "Fast_cycling_hot"],
    target_labels=["CAFs myCAF like s5_hot"],
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    filter_fun=lambda x: x['cellphone_pvals'] <= 0.001,
    # expression magnitude
    orderby="lr_means",
    orderby_ascending=False,
    top_n=30,
    figure_size=(6, 8),
    cmap="Reds",
    size_range=(1, 5)
)
print(lr_plot)
lr_plot.save('figures/G0_and_fast_to_CAFs myCAF like s5_lr_plot.pdf')

lr_plot = li.pl.dotplot(
    adata=adata,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    source_labels=["G0_arrested_hot", "Fast_cycling_hot"],
    target_labels=["PVL_Immature s2_hot"],
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    filter_fun=lambda x: x['cellphone_pvals'] <= 0.001,
    # expression magnitude
    orderby="lr_means",
    orderby_ascending=False,
    top_n=30,
    figure_size=(6, 8),
    cmap="Reds",
    size_range=(1, 5)
)
print(lr_plot)
lr_plot.save('figures/G0_and_fast_to_PVL_Immature s2_hot_lr_plot.pdf')

lr_plot = li.pl.dotplot(
    adata=adata,
    colour="lr_means",
    size="cellphone_pvals",
    inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
    # We choose only the cell types which we wish to plot
    source_labels=["G0_arrested_hot", "Fast_cycling_hot"],
    target_labels=["Myeloid_c2_LAM2_APOE_hot"],
    # as this type of methods tends to result in large numbers
    # of predictions, we can also further order according to
    filter_fun=lambda x: x['cellphone_pvals'] <= 0.001,
    # expression magnitude
    orderby="lr_means",
    orderby_ascending=False,
    top_n=30,
    figure_size=(6, 8),
    cmap="Reds",
    size_range=(1, 5)
)
print(lr_plot)
lr_plot.save('figures/G0_and_fast_to_Myeloid_c2_LAM2_APOE_hot_lr_plot.pdf')