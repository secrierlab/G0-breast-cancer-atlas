remove(list = ls())

# load libraries
library(dplyr)
library(Seurat)
library(irGSEA)
library(doMC)
library(readxl)

project_dir <- "~/BreastCancerG0arrest"
setwd(paste0(project_dir, "/05_GeneRegulatoryNetworks"))

# Load the data:
seurat_obj <- readRDS("../02_G0arrestInMalignantCells/data/malignant_integrated.rds")
# subset invasive types:
seurat_obj <- subset(
  seurat_obj,
  subset = (disease %in% "cancer") &
    (seurat_obj$QuiescenceStatus %in% c("Quiescent", "Proliferating", "Slow-cycling")) &
    !(seurat_obj$patient %in% c("Patient.12", "Patient.13", "Patient.14", "Patient.15")) & 
    !(seurat_obj$type %in% c("DCIS", "neoplasm")))

counts <- GetAssayData(seurat_obj, layer = "counts")
meta_data <- seurat_obj@meta.data
umap_loadings <- seurat_obj@reductions$umap

seurat_obj <- CreateSeuratObject(counts = counts, meta.data = meta_data)
seurat_obj@reductions$umap <- umap_loadings
remove(counts, meta_data)

Idents(seurat_obj) <- seurat_obj$QuiescenceStatus

seurat_obj <- seurat_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# factor the QuiescenceStatus for the plot
seurat_obj$QuiescenceStatus <- factor(
  seurat_obj$QuiescenceStatus,
  levels = c("Proliferating", "Slow-cycling", "Quiescent")
)

## senescence signature:
senescence_geneset <- c("CDC27", "E2F2", "SCMH1", "MRE11", "MAP2K3", "MAPK9", "ANAPC4", "MAP2K4", "MAP4K4", "RPS6KA2", 
                        "UBE2D1", "EED", "MAP2K7", "TNRC6C", "MAPKAPK5", "ANAPC5", "TNRC6A", "TINF2", "AGO1", "CDC23", 
                        "CABIN1", "MAPK1", "HIRA", "TNRC6B", "E2F1", "RBBP7", "MAPK3", "ACD", "NBN", "CCNE1", "FZR1", 
                        "ERF", "CDK6", "H2AZ2", "EZH2", "MAPK8", "UBE2S", "MAP2K6", "NFKB1", "MAPK10", "ANAPC15", 
                        "CDKN1B", "PHC1", "ASF1A", "MAPK14", "E2F3", "LMNB1", "RAD50", "TFDP2", "MAPKAPK3", "IL1A", 
                        "RPS6KA1", "UBN1", "RNF2", "CDKN2C", "CDK2", "H1-3", "H1-1", "H2BC11", "CDKN1A", "ID1", "AGO3", 
                        "POT1", "CDKN2D", "CDC16", "H3-3B", "KDM6B", "TERF2", "CCNA1", "PHC2", "AGO4", "ETS1", "CDK4", 
                        "MDM2", "IL6", "TXN", "HMGA1", "RB1", "MINK1", "TP53", "ANAPC11", "CBX8", "CBX4", "RPS27A", 
                        "CCNA2", "H2BC1", "TERF1", "CDKN2B", "CDKN2A", "ATM", "HMGA2", "UBC", "VENTX", "ANAPC1", 
                        "TNIK", "MOV10", "ETS2", "H2BC5", "H4C8", "RBBP4", "MAPKAPK2", "H3-3A", "IGFBP7", "ANAPC10", 
                        "ANAPC16", "MAPK7", "TERF2IP", "H3-4", "BMI1", "H1-4", "STAT3", "CXCL8", "UBE2E1", "UBB", 
                        "FOS", "IFNB1", "CEBPB", "KAT5", "RELA", "PHC3", "CBX2", "UBE2C", "CCNE2", "ANAPC2", "CDC26", 
                        "RPS6KA3", "JUN", "SUZ12", "H2AC6", "H2BC4", "EHMT1", "EP400", "H3C13", "CBX6", "H2AC20", 
                        "H1-5", "H2BC21", "H2BC13", "MAPK11", "SP1", "H1-2", "H2AX", "H1-0", "ANAPC7", "H2AC7", 
                        "H2BC26", "H4C3", "H3C12", "H4C11", "H3C4", "MAP3K5", "H4C16", "H2BC12", "TFDP1", "MDM4", 
                        "H3C14", "H3C15", "RING1", "EHMT2", "UBA52", "H2BC15", "H2BC12L", "H2AJ", "H4C15", "H4C14", 
                        "H4C12", "H2BC14", "H2BC8", "H3C 8", "H2AB1", "H2BC6", "H4C6", "H2BC17", "H3C6", "H4C13", 
                        "H3C11", "H2BC9", "H3C1", "H4C9", "H2AC14", "H2BC3", "H4C5", "H2AC8", "H4C4", "H2BC7", 
                        "H3C7", "H2AC4", "H2BC10", "H4C1", "H4C2", "H3C10", "RPS6KA1", "E2F2", "MIR24-2", "MIR24-1", 
                        "TINF2", "H3-4", "H2BC26", "H3C2", "H3C3", "H2AC18", "H2AC19")

# dormancy signature (0.1186/s13058-022-01503-5:
dormancy_geneset_Ren <- c("CFH", "GAS6", "MME", "OGN", "POSTN", "PDGFRB", "DHRS3", 
                      "MGP", "ALDH1A1", "ALDH3A1", "PRELP", "AAK1", "THBD", 
                      "BCL2L11", "GLIS1")
# dormancy signature (10.1186/s13058-021-01416-9):
dormancy_geneset_Ruth <- c("BHLHE41", "TGFBR3", "TGFBR2", "THBS1", "FOXM1", "PLAUR")

# dormancy signature (10.1186/s13058-014-0407-9):
dormancy_geneset_Cheng <- c(
  "SKAP2", "LIMS1", "PLXDC1", "MAFB", "SHOX2", "LOH3CR2A", "NNMT", "TNC", "MRC2", "LAMA4",
  "DYRK2", "SEMA5A", "SHOX2", "SNAI2", "SEC23A", "GNS", "HTRA1", "DDEF1IT1", "LAMB1", "COL5A1",
  "JAG1", "COL1A1", "EMILIN1", "PHTF2", "COL5A1", "CLEC11A", "JAG1", "ANGPTL2", "COL3A1", "FAM65A",
  "COL6A3", "THY1", "PKD2", "NID2", "VCAN", "TGFB1I1", "CALM1", "CXCR7", "IL13RA1", "COL1A2",
  "THY1", "PDE8A", "ATP2B1", "ITGA5", "RABGGTA", "LRRC32", "CDKN1A", "POSTN", "TWIST1", "ATP2B1",
  "M6PRBP1", "ADAM9", "DYRK2", "HSPG2", "TPP1", "AEBP1", "MRPL52"
)

# dormancy signature (10.1371/journal.pone.0035569):
dormancy_geneset_Kim <- c(
  "ACVR1", "ADAM10", "AMOT", "BHLHE41", "COL1A1", "COL4A5", "CTSD", "DDR1", "EPHA5", 
  "GATA6", "HIST1H2BK", "IGFBP5", "MMP2", "NR2F1", "P4HA1", "SOX9", "SREBF1", 
  "STAT3", "TGFB2", "THBS1", "TP53", "TPM1"
)

# dormancy signature (10.1038/s41556-020-0474-3):
dormancy_geneset_Montagner <- c(
  "S1PR3", "PENK", "CHST8", "CACNA2D2", "ENTPD2", "SFRP2", "RARRES1", "HMGCS2", "CPZ", "SRPX2",
  "CORO2B", "SLC17A7", "TSPAN11", "LRRC17", "FAM101A", "EMID1", "RFTN2", "QPCT", "COL14A1", "RAMP2",
  "DPT", "FGL2", "PTGFR", "TCEAL1", "HEYL", "IL16", "CTSF", "WIF1", "SLC34A2", "C1QTNF5", "ECM2",
  "POSTN", "PCDHB12", "ABCA8A", "USP2", "PID1", "SHISA2", "SOD3", "FLT3L", "RARRES2", "CERCAM",
  "DIO2", "COL3A1", "CILP", "CRYAA", "MFAP2", "TRF", "CXCL12", "CXXC5", "EFS"
)

# Create a list of gene sets
gene_sets <- list(
  "Dormancy_Ren" = dormancy_geneset_Ren,
  "Dormancy_Ruth" = dormancy_geneset_Ruth,
  "Dormancy_Cheng" = dormancy_geneset_Cheng,
  "Dormancy_Kim" = dormancy_geneset_Kim,
  "Dormancy_Montagner" = dormancy_geneset_Montagner
)

seurat_obj <- irGSEA.score(
  object = seurat_obj,
  assay = "RNA",
  slot = "data",
  seeds = 123,
  ncores = 2,
  min.cells = 3,
  method = "AUCell",
  min.feature = 0,
  custom = T,
  geneset = gene_sets,
  msigdb = F,
  geneid = "symbol",
  aucell.MaxRank = 2000,
  ucell.MaxRank = 2000,
  kcdf = 'Gaussian'
)

result_dge <- irGSEA.integrate(
  object = seurat_obj,
  method = "AUCell",
  group.by = "QuiescenceStatus",
  metadata = NULL,
  col.name = NULL
)

setwd("figures")

pdf("irGSEA_heatmap_dormancy.pdf")
irGSEA.bubble(object = result_dge,
              method = "AUCell",
              cluster.color = c("#1B7837","#BBBBBB","#762A83"),
              direction.color = c("#364B9A", "#A50026"),
              significance.color = c("#f7f7f7", "black"),
              cluster.levels = c("Proliferating", "Slow-cycling", "Quiescent"),
              top = 50)
dev.off()

## Senescence-------------------------------------------------------------------
seurat_obj <- irGSEA.score(
  object = seurat_obj,
  assay = "RNA",
  slot = "data",
  seeds = 123,
  ncores = 2,
  min.cells = 3,
  method = "AUCell",
  min.feature = 0,
  custom = T,
  geneset = list("Senescence" = senescence_geneset),
  msigdb = F,
  geneid = "symbol",
  aucell.MaxRank = 2000,
  ucell.MaxRank = 2000,
  kcdf = 'Gaussian'
)

result_dge <- irGSEA.integrate(
  object = seurat_obj,
  method = "AUCell",
  group.by = "QuiescenceStatus",
  metadata = NULL,
  col.name = NULL
)

pdf("irGSEA_heatmap_senescence.pdf")
irGSEA.densityheatmap(object = seurat_obj, 
                      method = "AUCell",
                      group.by = "QuiescenceStatus",
                      show.geneset = "Senescence",
                      heatmap_width = 6,
                      heatmap_height = 5)
dev.off()

setwd(project_dir)