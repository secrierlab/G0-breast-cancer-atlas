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
sys.path.append(r"/Users/cenkcelik/BreastCancerG0arrest/10_SpatialAnalysis/SpottedPy/spottedPy")
# Import the package
import spottedpy as sp

# Set the working directory
os.chdir(r"/Users/cenkcelik/BreastCancerG0arrest/10_SpatialAnalysis")

# Load the anndata object
adata = sc.read_h5ad("newly_deconvolved_data_with_hotspots1410.h5ad")

# Signature labels
signature_labels = [
    "Hypoxia_hallmarks",
    "Angiogenesis_hallmarks",
    "EMT_hallmarks",
    "Hybrid_EMT",
    "Mesenchymal",
    "Epithelial",
    "UPR_hallmarks",
    "IRE1_pathway",
    "PERK_pathway",
    "ATF6_pathway",
    "Semaphorin_Plexin",
    "Necrosis",
    "Apoptosis",
    "Complement"
]

quiescence_labels = ["G0_arrested","Fast_cycling", "Slow_cycling"]

# List of hypoxia hallmark genes
Hypoxia_hallmarks = [
    "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2",
    "BCAN", "BCL2", "BGN", "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CCNG2", "CCRN4L", "CDKN1A", "CDKN1B", "CDKN1C",
    "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", "CTGF", "CXCR4", "CXCR7", "CYR61", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", "DUSP1",
    "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1L", "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2",
    "FOXO3", "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1",
    "HDLBP", "HEXA", "HK1", "HK2", "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", "ILVBL", "INHA", "IRS2",
    "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE", "LDHA", "LDHC", "LOX", "LXN", "MAFF",
    "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NR3C1", "P4HA1",
    "P4HA2", "PAM", "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1",
    "PKLR", "PKP1", "PLAC8", "PLAUR", "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", "PRKCA", "PRKCDBP", "PTRF",
    "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", "SCARB1", "SDC2", "SDC3", "SDC4", "SELENBP1", "SERPINE1", "SIAH2", "SLC25A1", "SLC2A1",
    "SLC2A3", "SLC2A5", "SLC37A4", "SLC6A6", "SRPX", "STBD1", "STC1", "STC2", "SULT2B1", "TES", "TGFB3", "TGFBI", "TGM2", "TIPARP", "TKTL1",
    "TMEM45A", "TNFAIP3", "TPBG", "TPD52", "TPI1", "TPST2", "UGP2", "VEGFA", "VHL", "VLDLR", "WISP2", "WSB1", "XPNPEP1", "ZFP36", "ZNF292"
]

# List of angiogenesis hallmark genes
Angiogenesis_hallmarks = [
    "APOH", "APP", "CCND2", "COL3A1", "COL5A2", "CXCL6", "FGFR1", "FSTL1", "ITGAV", "JAG1", "JAG2", "KCNJ8", "LPL", "LRPAP1", "LUM", "MSX1",
    "NRP1", "OLR1", "PDGFA", "PF4", "PGLYRP1", "POSTN", "PRG2", "PTK2", "S100A4", "SERPINA5", "SLCO2A1", "SPP1", "STC1", "THBD", "TIMP1",
    "TNFRSF21", "VAV2", "VCAN", "VEGFA", "VTN"
]

# List of epithelial-mesenchymal transition (EMT) hallmark genes
EMT_hallmarks = [
    "ABI3BP", "ACTA2", "ADAM12", "ANPEP", "APLP1", "AREG", "BASP1", "BDNF", "BGN", "BMP1", "CADM1", "CALD1", "CALU", "CAP2", "CAPG", "CD44", "CD59",
    "CDH11", "CDH2", "CDH6", "COL11A1", "COL12A1", "COL16A1", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL5A3",
    "COL6A2", "COL6A3", "COL7A1", "COL8A2", "COMP", "COPA", "CRLF1", "CTGF", "CTHRC1", "CXCL1", "CXCL12", "CXCL6", "CYR61", "DAB2", "DCN", "DKK1",
    "DPYSL3", "DST", "ECM1", "ECM2", "EDIL3", "EFEMP2", "ELN", "EMP3", "ENO2", "FAP", "FAS", "FBLN1", "FBLN2", "FBLN5", "FBN1", "FBN2",
    "FERMT2", "FGF2", "FLNA", "FMOD", "FN1", "FOXC2", "FSTL1", "FSTL3", "FUCA1", "FZD8", "GADD45A", "GADD45B", "GAS1", "GEM", "GJA1", "GLIPR1",
    "GLT25D1", "GPC1", "GPX7", "GREM1", "HTRA1", "ID2", "IGFBP2", "IGFBP3", "IGFBP4", "IL15", "IL32", "IL6", "IL8", "INHBA", "ITGA2", "ITGA5",
    "ITGAV", "ITGB1", "ITGB3", "ITGB5", "JUN", "LAMA1", "LAMA2", "LAMA3", "LAMC1", "LAMC2", "LEPRE1", "LGALS1", "LOX", "LOXL1", "LOXL2",
    "LRP1", "LRRC15", "LUM", "MAGEE1", "MATN2", "MATN3", "MCM7", "MEST", "MFAP5", "MGP", "MMP1", "MMP14", "MMP2", "MMP3", "MSX1", "MXRA5",
    "MYL9", "MYLK", "NID2", "NNMT", "NOTCH2", "NT5E", "NTM", "OXTR", "PCOLCE", "PCOLCE2", "PDGFRB", "PDLIM4", "PFN2", "PLAUR", "PLOD1",
    "PLOD2", "PLOD3", "PMEPA1", "PMP22", "POSTN", "PPIB", "PRRX1", "PRSS2", "PTHLH", "PTX3", "PVR", "QSOX1", "RGS4", "RHOB", "SAT1",
    "SCG2", "SDC1", "SDC4", "SERPINE1", "SERPINE2", "SERPINH1", "SFRP1", "SFRP4", "SGCB", "SGCD", "SGCG", "SLC6A8", "SLIT2", "SLIT3",
    "SNAI2", "SNTB1", "SPARC", "SPOCK1", "SPP1", "TAGLN", "TFPI2", "TGFB1", "TGFBI", "TGFBR3", "TGM2", "THBS1", "THBS2", "THY1","TIMP1",
    "TIMP3", "TNC", "TNFAIP3", "TNFRSF11B", "TNFRSF12A", "TPM1", "TPM2", "TPM4", "VCAM1", "VCAN", "VEGFA", "VEGFC", "VIM", "WIPF1", "WNT5A"
]

# pEMT, E and M lists from Malagoli Tagliazucchi, G. et al. (2023) Nat Commun https://doi.org/10.1038/s41467-023-36439-7
Hybrid_EMT = ["PDPN", "ITGA5", "ITGA6", "TGFBI", "LAMC2", "MMP10", "LAMA3", "CDH13", "SERPINE1", "P4HA2", "TNC", "MMP1"]

Mesenchymal = ["VIM", "FOXC2", "SNAI1", "SNAI2", "TWIST1", "FN1", "ITGB6", "MMP2", "MMP3", "MMP9", "SOX10", "GCS", "ZEB1", "ZEB2", "TWIST2"]

Epithelial = ["CDH1", "DSP", "OCLN", "CRB3"]

# List of endoplasmic reticulum (ER) stress response genes
UPR_hallmarks = [
    "ALDH18A1", "ARFGAP1", "ASNS", "ATF3", "ATF4", "ATF6", "ATP6V0D1", "BAG3", "BANF1", "CALR", "CCL2", "CEBPB", "CEBPG", "CHAC1", "CKS1B",
    "CNOT2", "CNOT4", "CNOT6", "CXXC1", "DCP1A", "DCP2", "DCTN1", "DDIT4", "DDX10", "DKC1", "DNAJA4", "DNAJB9", "DNAJC3", "EDC4", "EDEM1", 
    "EEF2", "EIF2AK3", "EIF2S1", "EIF4A1", "EIF4A2", "EIF4A3", "EIF4E", "EIF4EBP1", "EIF4G1", "ERN1", "ERO1L", "EXOC2", "EXOSC1", "EXOSC10",
    "EXOSC2", "EXOSC4", "EXOSC5", "EXOSC9", "FKBP14", "FUS", "GEMIN4", "GOSR2", "H2AFX", "HERPUD1", "HSP90B1", "HSPA5", "HSPA9", "HYOU1", 
    "IARS", "IFIT1", "IGFBP1", "IMP3", "KDELR3", "KHSRP", "KIF5B", "LSM1", "LSM4", "MTHFD2", "NFYA", "NFYB", "NHP2", "NOLC1", "NOP14", 
    "NOP56", "NPM1", "OBFC2A", "PAIP1", "PARN", "PDIA5", "PDIA6", "POP4", "PREB", "PSAT1", "RPS14", "RRP9", "SDAD1", "SEC11A", "SEC31A", 
    "SERP1", "SHC1", "SKIV2L2", "SLC1A4", "SLC30A5", "SLC7A5", "SPCS1", "SPCS3", "SRPR", "SRPRB", "SSR1", "STC2", "TARS", "TATDN2", "TSPYL2", 
    "TTC37", "TUBB2A", "VEGFA", "WFS1", "WIPI1", "XBP1", "XPOT", "YIF1A", "YWHAZ", "ZBTB17"
]

IRE1_pathway = [
    "AGR2", "BAK1", "BAX", "BCL2L11", "BFAR", "COPS5", "DAB2IP", "DDRGK1", "DNAJB9", "ERN1", "ERN2", "FICD", "HSPA5", "PARP16", "PTPN1",
    "TMEM33", "UFL1", "VAPB", "XBP1"
]

PERK_pathway = [
    "ABCA7", "AGR2", "ATF4", "BOK", "DDIT3", "EIF2AK3", "EIF2S1", "HSPA5", "NCK1", "NCK2", "NFE2L2", "PPP1R15A", "PPP1R15B", "PTPN1",
    "PTPN2", "QRICH1", "TMED2", "TMEM33"
]

ATF6_pathway = ["ATF6", "ATF6B", "DDIT3", "HSPA5", "MBTPS1", "MBTPS2", "WFS1"]

Semaphorin_Plexin = [
    "PLXNC1", "SEMA3A", "SEMA6C", "SEMA6B", "SEMA4F", "SEMA4D", "SEMA4B", "SEMA3C", "ECE1", "EDN1", "EDNRA", "ERBB2", "SEMA3D",
    "PLXND1", "FLNA", "SH3BP1", "PLXNB2", "GDNF", "KDR", "RHOA", "ARHGDIA", "MET", "NCAM1", "PLXNA1", "PLXNA2", "PLXNB1", "PLXNB3",
    "SEMA5B", "SEMA4C", "PLXNA3", "SEMA3G", "SEMA6A", "SEMA4G", "RAC1", "SEMA3F", "SEMA4A", "SEMA3B", "SEMA6D", "SEMA7A", "NRP2",
    "NRP1", "SEMA5A", "PLXNA4", "HAND2", "SEMA3E", "FARP2"
]

Necrosis = ["RIPK3", "MLKL", "FAS", "FASLG", "TLR3", "TNF", "RIPK1", "FADD"]

Apoptosis = [
    "AKT3", "IRAK3", "CHP1", "CHUK", "CSF2RB", "DFFA", "DFFB", "ENDOG", "AKT1", "AKT2", "ENDOD1", "PIK3R5", "APAF1", "BIRC2", "BIRC3", "XIAP", "FAS",
    "IKBKB", "IL1A", "IL1B", "IL1R1", "IL1RAP", "FASLG", "IL3", "IL3RA", "IRAK1", "IRAK2", "MYD88", "ATM", "NFKB1", "NFKBIA", "NGF", "NTRK1", "IRAK4",
    "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "CYCS", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", "PRKACA", "PRKACB", "PRKACG",
    "PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B", "PRKX", "BAD", "BAX", "BCL2", "RELA", "BCL2L1", "BID", "CHP2", "TNF", "TNFRSF1A", "TP53", "TRAF2",
    "CAPN1", "CAPN2", "CASP3", "CASP6", "CASP7", "CASP8", "CASP9", "CASP10", "PIK3R3", "IKBKG", "TRADD", "RIPK1", "TNFSF10", "FADD", "TNFRSF10D",
    "TNFRSF10C", "TNFRSF10B", "TNFRSF10A", "CFLAR", "MAP3K14", "AIFM1", "EXOG"
]

Complement = [
    "ACTN2", "ADAM9", "ADRA2B", "AKAP10", "ANG", "ANXA5", "APOA4", "APOBEC3F", "APOBEC3G", 
    "APOC1", "ATOX1", "BRPF3", "C1QA", "C1QC", "C1R", "C1S", "C2", "C3", "C4BPB", "C9", 
    "CA2", "CALM1", "CALM3", "CASP1", "CASP10", "CASP3", "CASP4", "CASP5", "CASP7", "CASP9", 
    "CBLB", "CCL5", "CD36", "CD40LG", "CD46", "CD55", "CD59", "CDA", "CDH13", "CDK5R1", 
    "CEBPB", "CFB", "CFH", "CLU", "COL4A2", "CP", "CPM", "CR1", "CR2", "CSRP1", "CTSB", 
    "CTSC", "CTSD", "CTSH", "CTSL", "CTSV", "CTSO", "CTSS", "CXCL1", "DGKG", "DGKH", 
    "DOCK10", "DOCK4", "DOCK9", "DPP4", "DUSP5", "DUSP6", "DYRK2", "EHD1", "ERAP2", "F10", 
    "F2", "F3", "F5", "F7", "F8", "FCER1G", "FCN1", "FDX1", "FN1", "FYN", "GATA3", "GCA", 
    "GMFB", "GNAI2", "GNAI3", "GNB2", "GNB4", "GNG2", "GNGT2", "GP1BA", "GP9", "GPD2", 
    "GRB2", "GZMA", "GZMB", "GZMK", "HNF4A", "HPCAL4", "HSPA1A", "HSPA5", "IL6", "IRF1", 
    "IRF2", "IRF7", "ITGAM", "ITIH1", "JAK2", "KCNIP2", "KCNIP3", "KIF2A", "KLK1", "KLKB1", 
    "KYNU", "L3MBTL4", "LAMP2", "LAP3", "LCK", "LCP2", "LGALS3", "LGMN", "LIPA", "LRP1", 
    "LTA4H", "LTF", "LYN", "MAFF", "ME1", "MMP12", "MMP13", "MMP14", "MMP15", "MMP8", 
    "MT3", "NOTCH4", "OLR1", "PCLO", "PCSK9", "PDGFB", "PDP1", "PFN1", "CPQ", "PHEX", 
    "PIK3CA", "PIK3CG", "PIK3R5", "PIM1", "PLA2G4A", "PLA2G7", "PLAT", "PLAUR", "PLEK", 
    "PLG", "PLSCR1", "PPP2CB", "PPP4C", "PRCP", "PRDM4", "PREP", "PRKCD", "PRSS3", 
    "PRSS36", "PSEN1", "PSMB9", "RABIF", "RAF1", "RASGRP1", "RCE1", "RHOG", "RNF4", 
    "S100A12", "S100A13", "S100A9", "SCG3", "MSRB1", "SERPINA1", "SERPINB2", "SERPINC1", 
    "SERPINE1", "SERPING1", "SH2B3", "SIRT6", "SPOCK2", "SRC", "STX4", "TFPI2", "TIMP1", 
    "TIMP2", "TMPRSS6", "TNFAIP3", "USP14", "USP15", "USP16", "USP8", "VCPIP1", "WAS", 
    "XPNPEP1", "ZEB1", "ZFPM2", "RBSN"
]


# filter upregulated genes for those that are in the adata.var_names
Hypoxia_hallmarks = [gene for gene in Hypoxia_hallmarks if gene in adata.var_names]
Angiogenesis_hallmarks = [gene for gene in Angiogenesis_hallmarks if gene in adata.var_names]
EMT_hallmarks = [gene for gene in EMT_hallmarks if gene in adata.var_names]
Hybrid_EMT = [gene for gene in Hybrid_EMT if gene in adata.var_names]
Mesenchymal = [gene for gene in Mesenchymal if gene in adata.var_names]
Epithelial = [gene for gene in Epithelial if gene in adata.var_names]
UPR_hallmarks = [gene for gene in UPR_hallmarks if gene in adata.var_names]
IRE1_pathway = [gene for gene in IRE1_pathway if gene in adata.var_names]
PERK_pathway = [gene for gene in PERK_pathway if gene in adata.var_names]
ATF6_pathway = [gene for gene in ATF6_pathway if gene in adata.var_names]
Semaphorin_Plexin = [gene for gene in Semaphorin_Plexin if gene in adata.var_names]
Necrosis = [gene for gene in Necrosis if gene in adata.var_names]
Apoptosis = [gene for gene in Apoptosis if gene in adata.var_names]
Complement = [gene for gene in Complement if gene in adata.var_names]

# Define a list of hallmarks and corresponding column names
hallmarks = [
    (Hypoxia_hallmarks, "Hypoxia_hallmarks"),
    (Angiogenesis_hallmarks, "Angiogenesis_hallmarks"),
    (EMT_hallmarks, "EMT_hallmarks"),
    (Hybrid_EMT, "Hybrid_EMT"),
    (Mesenchymal, "Mesenchymal"),
    (Epithelial, "Epithelial"),
    (UPR_hallmarks, "UPR_hallmarks"),
    (IRE1_pathway, "IRE1_pathway"),
    (PERK_pathway, "PERK_pathway"),
    (ATF6_pathway, "ATF6_pathway"),
    (Semaphorin_Plexin, "Semaphorin_Plexin"),
    (Necrosis, "Necrosis"),
    (Apoptosis, "Apoptosis"),
    (Complement, "Complement")
]

# Process each batch and calculate scores
for batch in adata.obs["batch"].unique():
    batch_subset = adata[adata.obs["batch"] == batch]
    
    for hallmark_genes, score_name in hallmarks:
        # Calculate scores for the hallmark genes
        sc.tl.score_genes(batch_subset, hallmark_genes, score_name=score_name)
        
        # Update the scores in the main dataframe
        adata.obs.loc[adata.obs["batch"] == batch, score_name] = batch_subset.obs[score_name].values

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
scale_scores(adata, signature_labels)

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
            neighbours_parameters=12,
            relative_to_batch=True,
            number_hotspots=True,
            filter_columns=filter_columns,
            filter_value=filter_value
        )
    return adata

tumour_signature_labels = ["Semaphorin_Plexin", "Hybrid_EMT", "Mesenchymal", "Epithelial", "UPR_hallmarks", "IRE1_pathway", "PERK_pathway", "ATF6_pathway", "Complement"]

adata = calculate_hotspots(adata, [label for label in signature_labels if label not in tumour_signature_labels])
adata = calculate_hotspots(adata, tumour_signature_labels, filter_columns="tumour_cells", filter_value=1)

# Generate signature_labels_hotspots by appending "_hot" to each label in signature_labels
signature_labels_hotspots = [label + "_hot" for label in signature_labels]

## G0 vs Fast cycling------------------------------------------------------------
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
                       fig_size=(8, 6),
                       file_save="figures/signature_main_custom_scatter_",
                       compare_distribution_metric="median")

# Define a list of tuples with column names and save paths
plots = [
    ("EMT_hallmarks_hot", "EMT_hallmarks_hotspots.pdf"),
    ("Hybrid_EMT_hot", "Hybrid_EMT_hotspots.pdf"),
    ("UPR_hallmarks_hot", "UPR_hallmarkss_hotspots.pdf"),
    ("Angiogenesis_hallmarks_hot", "Angiogenesis_hallmarks_hotspots.pdf"),
    ("Hypoxia_hallmarks_hot", "Hypoxia_hallmarks_hotspots.pdf"),
    ("Semaphorin_Plexin_hot", "Semaphorin_Plexin_hotspots.pdf"),
    ("Necrosis_hot", "Necrosis_hotspots.pdf"),
    ("Apoptosis_hot", "Apoptosis_hotspots.pdf"),
    ("IRE1_pathway_hot", "IRE1_pathway_hotspots.pdf"),
    ("PERK_pathway_hot", "PERK_pathway_hotspots.pdf"),
    ("ATF6_pathway_hot", "ATF6_pathway_hotspots.pdf")
]

# Loop through the list and plot each hotspot
for column_name, save_path in plots:
    sp.plot_hotspots(
        adata,
        column_name=column_name,
        save_path=save_path
    )

# differences per batch (i.e. slide)
data=sp.plot_bubble_chart_by_batch(distances,
                                   primary_variable_value="G0_arrested_hot",
                                   comparison_variable_values=signature_labels_hotspots,
                                   reference_variable="Fast_cycling_hot",
                                   save_path="figures/signature_bubble_chart_by_batch.pdf",
                                   pval_cutoff=0.05,
                                   fig_size=(8,7))