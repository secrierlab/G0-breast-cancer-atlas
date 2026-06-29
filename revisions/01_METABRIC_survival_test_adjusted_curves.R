rm(list = ls())
library(dplyr)
library(data.table)
library(genefu)
library(survival)
library(survminer)
library(ggplot2)
library(GSVA)
library(maxstat)
library(DGEobj.utils)

projectDir <- "~/BreastCancerG0arrest/"
setwd(paste0(projectDir, "/06_Survival"))

# METABRIC cohort----
# download the METABRIC data
link = "https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz"
curl::curl_download(url = link, destfile = "data/brca_metabric.tar.gz")
# untar the file
system("tar -xvzf data/brca_metabric.tar.gz -C data/")
# read the expression data
expression_data_full <- fread("data/brca_metabric/data_mrna_illumina_microarray.txt") %>% 
  data.frame(check.names = FALSE)
# average replicate probes
expression_data_full <- expression_data_full %>% 
  group_by(Hugo_Symbol) %>% 
  summarise_all(mean) %>% 
  as.data.frame()
# remove the gene names
rownames(expression_data_full) <- expression_data_full$Hugo_Symbol
expression_data_full$Hugo_Symbol <- NULL
expression_data_full$Entrez_Gene_Id <- NULL
# replace the missing values with zero
expression_data_full[is.na(expression_data_full)] <- 0
# transpose the data
expression_matrix_full <- t(expression_data_full)
## log2 + 1 transformation
# expression_data_full <- log2(expression_matrix_full + 1) %>% as.data.frame()
expression_data_full <- convertCounts(
  expression_matrix_full, 
  unit = "TPM", 
  log = TRUE, 
  geneLength = rowMeans(expression_matrix_full)
) %>% as.data.frame()
# save the expression data
# save(expression_data_full, file = "data/expression_data_full.RData")
# load the expression data
# load("data/expression_data_full.RData")

expression_data = expression_data_full

# G0 scoring gene sets:
upregulated_common <- c(
  "HLA-B", "TACSTD2", "LCN2", "HLA-A", "CD24", "GSTK1", "PDLIM1", "SMIM22",
  "PPBP", "KLK5", "SLPI", "C15orf48", "CST6", "IL32", "ELF3", "ID3", "TFPI2",
  "SAA1", "JUN", "IFI27", "CYP1B1", "IFI6", "FOS", "TSPAN8", "CXCL16", "CYR61",
  "IGFBP3", "VSIG2", "AC006262.5", "PSCA", "XAGE2", "RARRES1", "FOSB", "NFKBIA",
  "NOV", "HPGD", "PMP22", "SAA2", "DUSP1", "BPGM", "HLA-DRB1", "HLA-DRB5", "UBD",
  "HLA-DPB1", "KRT16", "RP3-522D1.1", "PDZK1IP1", "OLR1", "BIRC3", "SEMA5A",
  "COL5A2", "NREP", "ATF3", "C6orf132", "ALDH1A1", "HIST1H2AC", "CLIC3",
  "RP11-95M15.1", "FAM134B", "ALPP", "RP11-54H7.4", "MT1E", "IFIT2", "NTN4",
  "LAMC2", "CRYAB", "CD82", "CTA-293F17.1", "MUC16", "LINC01133", "COL16A1",
  "SPNS2", "SORBS2", "MAF", "STXBP6", "DRAM1", "IFI44L", "SDCBP2", "CDH19",
  "TAC1", "IFIT1", "DDX60L", "UPK3A", "PHGR1", "C2orf54", "SPINK1", "IFI44",
  "LOXL3", "DAPP1", "LYPD6B", "GABBR2", "TTC9", "HLA-DQA1", "CTD-3252C9.4",
  "CXCL14", "CXCL3", "KCNS3", "MYLK", "MMP7", "SLC40A1", "HLA-DQA2", "NGFR",
  "TRIM31", "CXCL10", "ALOX5", "RP11-465N4.4", "DEGS2", "RARRES2", "TMEM47",
  "CCL2", "C1orf116", "CRYGS", "LCTL", "INHA", "MOXD1", "WNT4", "NEURL3",
  "NOXO1", "FRK", "GLIS3", "HIST1H2BD", "FAM3D", "TNNC1", "COL8A1", "PTPRZ1",
  "TNFSF10", "OLFML2A", "AMN", "HES2", "POF1B", "ARRB1", "HERC6", "UPK3B",
  "RGL1", "RGCC", "TCEA3", "LGALS2", "PRR15L", "ARHGAP24", "CHST15", "OSR2",
  "MUC20", "XAF1", "PVRL4", "ASCL2", "EPSTI1", "APOL1", "PAQR7", "CLEC7A",
  "RP11-736N17.10", "SCD5", "CTSS", "RND1", "ATP6V1B1", "PTN", "RHOU", "ANGPTL7",
  "SLITRK6"
)

downregulated_common <- c(
  "SLC45A2", "SYTL5", "GMPR", "TOX2", "TFF3", "PLEKHH2", "CPQ", "DTL", "RGS5",
  "MCM10", "TRPM1", "SCG2", "NTS", "CDCA2", "STMN3", "NCAPH", "HOXD-AS2",
  "DLX2", "RAB3A", "BUB1", "LINC00326", "TMEM98", "HAPLN1", "PGF", "KIF14",
  "ESCO2", "SPRY4-IT1", "SLCO4A1-AS1", "LRRC17", "RGS17", "TBX2", "ST3GAL5",
  "LONRF2", "STC1", "BUB1B", "NEK2", "SERPINF1", "CNIH3", "SNCA", "TROAP",
  "DLX1", "ENO2", "C6orf141", "BAALC", "HIST1H2BB", "RAB38", "LINC01419",
  "NCAPG", "CENPA", "TRIB2", "SH2B3", "HJURP", "PARVB", "DSCR8", "CCDC88A",
  "NUF2", "CEP55", "CAPN3", "ETV1", "ANLN", "PLK1", "MAGEB2", "TESC", "ST3GAL6",
  "NEFL", "DEPDC1", "APOC1", "IL13RA2", "INPP5F", "APOE", "TNC", "BIRC7", "F2R",
  "SLC5A3", "PRR11", "NES", "SPRY4", "ETV5", "DLGAP5", "GTSE1", "HMMR", "MLANA",
  "ASPM", "FOSL1", "IGFBP2", "FRMD4A", "SLC20A1", "ETV4", "PBK", "AURKA", "CDK1",
  "PCOLCE", "THBS2", "HMGA2", "PLN", "CTHRC1", "PMEL", "BIRC5", "CTAG2", "PRC1",
  "CDKN3", "BCYRN1", "TPX2", "PRAME", "DUSP6", "DUSP4", "GYPC", "IGFBP5", "CENPF",
  "PEG10", "TOP2A", "PTTG1", "VGF", "MT-ATP8", "LGALS1"
)

# Keep only the genes that are present in the expression data
upregulated_geneset <- upregulated_common[upregulated_common %in% colnames(expression_data)]
downregulated_geneset <- downregulated_common[downregulated_common %in% colnames(expression_data)]

# rownames to column
expression_data$PATIENT_ID <- rownames(expression_data)
# subset the signature genes
expression_data <- expression_data[, c("PATIENT_ID", upregulated_geneset, downregulated_geneset, "PLXNB1", "SEMA4D")]

# read the clinical data
clinical_data_patient <- fread("data/brca_metabric/data_clinical_patient.txt", skip = 4) %>% data.frame(check.names = F)
clinical_data_sample <- fread("data/brca_metabric/data_clinical_sample.txt", skip = 4) %>% data.frame(check.names = F)
# merge clinical data
clinical_data <- merge(clinical_data_patient, clinical_data_sample, by = "PATIENT_ID")
# relabel OS_STATUS
clinical_data$OS_STATUS <- ifelse(clinical_data$OS_STATUS == "1:DECEASED", 1, 0)
# relabel RFS_STATUS
clinical_data$RFS_STATUS <- ifelse(clinical_data$RFS_STATUS == "1:Recurred", 1, 0)
# convert months to years
clinical_data$OS_YEARS <- clinical_data$OS_MONTHS / 12
clinical_data$RFS_YEARS <- clinical_data$RFS_MONTHS / 12

# merge the expression and clinical data
metabric_data <- merge(expression_data, clinical_data, by = "PATIENT_ID")
# remove OS_YEARS > 25
metabric_data <- metabric_data %>% filter(OS_YEARS < 20)
metabric_data$TUMOR_STAGE <- factor(metabric_data$TUMOR_STAGE)
metabric_data$GRADE <- factor(metabric_data$GRADE)

## Calculate signature scores
G0_score_genesets <- list(upregulated_geneset, downregulated_geneset)
expression_matrix <- expression_data
expression_matrix$PATIENT_ID <- NULL
expression_matrix <- as.matrix(expression_matrix) %>% t()
expression_matrix_2 <- expression_matrix[G0_score_genesets %>% unlist(), metabric_data$PATIENT_ID]
expression_matrix_2 <- expression_matrix_2[!duplicated(rownames(expression_matrix_2)), ]

G0_score_parameters <- zscoreParam(
  exprData = expression_matrix_2, 
  geneSets=G0_score_genesets
)
G0_score <- GSVA::gsva(G0_score_parameters) %>% t() %>% data.frame()
metabric_data$G0_score = G0_score$X1 - G0_score$X2
hist(metabric_data$G0_score, breaks = 50)
# Cut off Q3
metabric_data$Quiescence = ifelse(metabric_data$G0_score > quantile(metabric_data$G0_score, 0.75), "G0_arrested",
                                  ifelse(metabric_data$G0_score < quantile(metabric_data$G0_score, 0.25), "Proliferating", "Slow_cycling"))

## Adjusted curves:
library(adjustedCurves)
library(riskRegression)
library(pammtools)

## Survival Plots----
### ER+ tumours----
# relapse free survival:
metabric_er_pos <- metabric_data %>% 
  filter(ER_STATUS == "Positive" & HER2_STATUS == "Negative") %>% 
  filter(Quiescence != "Slow_cycling")

metabric_er_pos$GRADE <- factor(metabric_er_pos$GRADE)
metabric_er_pos$Quiescence <- factor(metabric_er_pos$Quiescence)

outcome_model <- coxph(Surv(RFS_YEARS, RFS_STATUS) ~ GRADE + Quiescence, data = metabric_er_pos, x = TRUE)

adjusted_curve <- adjustedsurv(
  data = metabric_er_pos,
  variable = "Quiescence",
  ev_time = "RFS_YEARS",
  event = "RFS_STATUS",
  method = "direct",
  outcome_model = outcome_model,
  conf_int = TRUE,
  bootstrap = TRUE
)

adjusted_surv_quantile(adjusted_curve, p = 0.5, conf_int = TRUE, contrast = "ratio") -> median_survival

pdf("figures/METABRIC_G0_ER_signature_survival_RFS.pdf", width = 5, height = 5)
plot(adjusted_curve,
     conf_int = TRUE,
     median_surv_lines = TRUE,
     risk_table = TRUE,
     censoring_ind = "lines",
     risk_table_stratify = TRUE,
     risk_table_warn = FALSE,
     risk_table_digits = 0,
     pval = TRUE,
     title = "METABRIC ER+ G0 arrest signature",
     subtitle = paste0("p value = ", median_survival$p_value),
     legend.position = "top",
     custom_colors = c("#6A5ACD", "#2FBF71", "#8C8C8C"))
dev.off()

### TNBC tumours----
# relapse free survival:
metabric_tnbc <- metabric_data %>%
  filter(ER_STATUS == "Negative", PR_STATUS == "Negative", HER2_STATUS == "Negative") %>% 
  filter(Quiescence != "Slow_cycling")

metabric_tnbc$GRADE <- factor(metabric_tnbc$GRADE)
metabric_tnbc$Quiescence <- factor(metabric_tnbc$Quiescence)

outcome_model <- coxph(Surv(RFS_YEARS, RFS_STATUS) ~ GRADE + Quiescence, data = metabric_tnbc, x = TRUE)

adjusted_curve <- adjustedsurv(
  data = metabric_tnbc,
  variable = "Quiescence",
  ev_time = "RFS_YEARS",
  event = "RFS_STATUS",
  method = "direct",
  outcome_model = outcome_model,
  conf_int = TRUE
)

adjusted_surv_quantile(adjusted_curve, p = 0.5, conf_int = TRUE, contrast = "ratio") -> median_survival

pdf("figures/METABRIC_G0_TNBC_signature_survival_RFS.pdf", width = 5, height = 5)
plot(adjusted_curve,
     conf_int = TRUE,
     median_surv_lines = TRUE,
     risk_table = TRUE, 
     censoring_ind = "lines",
     risk_table_stratify = TRUE,
     risk_table_warn = FALSE,
     risk_table_digits = 0,
     legend.position = "top",
     title = "METABRIC TNBC G0 arrest signature",
     subtitle = paste0("p value = ", median_survival$p_value),
     custom_colors = c("#6A5ACD", "#2FBF71", "#8C8C8C"))
dev.off()

sessionInfo()