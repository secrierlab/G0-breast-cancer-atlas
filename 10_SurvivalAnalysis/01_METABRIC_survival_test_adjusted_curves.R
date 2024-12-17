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
setwd(paste0(projectDir, "/10_SurvivalAnalysis"))

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
load(paste0(projectDir, "/02_G0ArrestInMalignantCells/data/upregulated_common.RData"))
load(paste0(projectDir, "/02_G0ArrestInMalignantCells/data/downregulated_common.RData"))
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
G0_score_parameters <- zscoreParam(expression_matrix_2, G0_score_genesets)
G0_score <- gsva(G0_score_parameters) %>% t() %>% data.frame()
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

## Extended Data Fig. 2f:
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
     custom_colors = c("#762A83", "#1B7837", "#BBBBBB"))
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

## Extended Data Fig. 2g:
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
     custom_colors = c("#762A83", "#1B7837", "#BBBBBB"))
dev.off()

sessionInfo()