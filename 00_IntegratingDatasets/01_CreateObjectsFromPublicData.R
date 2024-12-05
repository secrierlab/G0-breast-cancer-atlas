library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggplot2)

projectDir <- "~/BreastCancerG0arrest"
setwd("00_IntegratingDatasets")

# Download data----
dir.create(paste0("rawdata"), recursive = F)

## Gao et al. 2021----
url <- "https://www.dropbox.com/sh/6k3zkgai4qlcvaq/AAAMsgKN88T6ZK43FX6frrHTa?dl=1"
destfile <- "rawdata/Data_Gao2021_Breast.zip"

## Pal et al. 2021----
url <- "https://www.dropbox.com/sh/0qskk6dlk80j690/AABfXa9meRfTt5OHLvPdV-9ua?dl=1"
destfile <- "rawdata/Data_Pal2021_Breast.zip"

## Qian et al. 2020----
url <- "https://www.dropbox.com/sh/nbx7v3om85wkfoq/AACpeZEZ4RNQwMW37Q7AHxExa?dl=1"
destfile <- "rawdata/Data_Qian2020_Breast.zip"

files <- list.files(path = "rawdata/", pattern = ".zip$")
outDir <- "rawdata"

for (i in files) {unzip(paste0("rawdata/", i), exdir = outDir)}

# Read 10X data----
dirs <- list.dirs(path = "rawdata", recursive = F, full.names = F)

for (x in dirs) {
  counts <- ReadMtx(mtx = paste0("rawdata/", x, "/matrix.mtx"), features = paste0("rawdata/", x, "/Genes.txt"), feature.column = 1, cells = paste0("rawdata/", x, "/barcodes.tsv"), skip.cell = 1)
  assign(x, CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = x))
}

Pal2021Group1@meta.data$orig.ident <- "Pal2021Group1"
Pal2021Group2@meta.data$orig.ident <- "Pal2021Group2"
Pal2021Group3@meta.data$orig.ident <- "Pal2021Group3"
Pal2021Group4@meta.data$orig.ident <- "Pal2021Group4"
Pal2021Group5@meta.data$orig.ident <- "Pal2021Group5"
Qian2020@meta.data$orig.ident <- "Qian2020"

rm(counts)

## Add metadata to the object----
gao.metadata <- read.csv(file = "rawdata/Gao2021/barcodes.tsv", sep = "\t")
colnames(gao.metadata) <- c("barcode","type","annotation","disease")
gao.metadata <- gao.metadata[ , 2:4]
gao.metadata$annotation <- sub("^$", "Unknown", gao.metadata$annotation)

Gao2021@meta.data$type <- gao.metadata$type
Gao2021@meta.data$patient <- gao.metadata$type
Gao2021@meta.data$type[Gao2021@meta.data$type == "DCIS1"] <- "DCIS"
Gao2021@meta.data$type[Gao2021@meta.data$type == "TNBC1"] <- "TNBC"
Gao2021@meta.data$type[Gao2021@meta.data$type == "TNBC2"] <- "TNBC"
Gao2021@meta.data$type[Gao2021@meta.data$type == "TNBC3"] <- "TNBC"
Gao2021@meta.data$patient[Gao2021@meta.data$patient == "DCIS1"] <- "Patient.1"
Gao2021@meta.data$patient[Gao2021@meta.data$patient == "TNBC1"] <- "Patient.2"
Gao2021@meta.data$patient[Gao2021@meta.data$patient == "TNBC2"] <- "Patient.3"
Gao2021@meta.data$patient[Gao2021@meta.data$patient == "TNBC3"] <- "Patient.4"
Gao2021@meta.data$annotation <- gao.metadata$annotation
Gao2021@meta.data$disease <- gao.metadata$disease

qian.metadata <- read.csv(file = "rawdata/Qian2020/barcodes.tsv", sep = "\t")
colnames(qian.metadata) <- c("barcode","sample","annotation", "subclone")
qian.metadata$annotation <- sub("^$", "Unknown", qian.metadata$annotation)
qian.metadata <- qian.metadata[ , 1:3]
qian.metadata$type <- "neoplasm"
Qian2020@meta.data$type <- qian.metadata$type
Qian2020@meta.data$patient <- qian.metadata$sample
Qian2020@meta.data$patient[Qian2020@meta.data$patient == 42] <- "Patient.5"
Qian2020@meta.data$patient[Qian2020@meta.data$patient == 43] <- "Patient.6"
Qian2020@meta.data$patient[Qian2020@meta.data$patient == 47] <- "Patient.7"
Qian2020@meta.data$patient[Qian2020@meta.data$patient == 49] <- "Patient.8"
Qian2020@meta.data$patient[Qian2020@meta.data$patient == 51] <- "Patient.9"
Qian2020@meta.data$patient[Qian2020@meta.data$patient == 53] <- "Patient.10"
Qian2020@meta.data$patient[Qian2020@meta.data$patient == 54] <- "Patient.11"
Qian2020@meta.data$annotation <- qian.metadata$annotation
Qian2020@meta.data$disease <- "cancer"

# pal
library(plyr)
pal.dirs <- dirs[2:6]

for (x in pal.dirs) {
  name <- paste0(x, "_metadata")
  assign(name, read.csv(file = paste0("rawdata/", x, "/barcodes.tsv"), sep = "\t"))
}

df <- Pal2021Group1@meta.data
df$Cell_names <- rownames(df)

df <- match_df(Pal2021Group1_metadata, df, on = "Cell_names")
Pal2021Group1@meta.data$type <- df$Sample
Pal2021Group1@meta.data$patient <- df$Sample
Pal2021Group1@meta.data$annotation <- df$Annotation
Pal2021Group1@meta.data$disease <- df$Is_Cancer

Pal2021Group1@meta.data$disease[Pal2021Group1@meta.data$disease == 0] <- "reference"
Pal2021Group1@meta.data$disease[Pal2021Group1@meta.data$disease == 1] <- "cancer"

Pal2021Group1@meta.data$type[Pal2021Group1@meta.data$type == "BRCA1_pre_neoplastic_0023"] <-  "normal"
Pal2021Group1@meta.data$type[Pal2021Group1@meta.data$type == "BRCA1_pre_neoplastic_0033"] <- "normal"
Pal2021Group1@meta.data$type[Pal2021Group1@meta.data$type == "BRCA1_pre_neoplastic_0090"] <- "normal"
Pal2021Group1@meta.data$type[Pal2021Group1@meta.data$type == "BRCA1_pre_neoplastic_0894"] <-"normal"
Pal2021Group1@meta.data$type[Pal2021Group1@meta.data$type == "ER_positive__LN_0056"] <- "ER"
Pal2021Group1@meta.data$type[Pal2021Group1@meta.data$type == "ER_positive_0001"] <- "ER"
Pal2021Group1@meta.data$type[Pal2021Group1@meta.data$type == "ER_positive_0025"] <- "ER"
Pal2021Group1@meta.data$type[Pal2021Group1@meta.data$type == "ER_positive_0029_7C"] <- "ER"
Pal2021Group1@meta.data$type[Pal2021Group1@meta.data$type == "ER_positive_0029_9C"] <- "ER"

Pal2021Group1@meta.data$patient[Pal2021Group1@meta.data$patient == "BRCA1_pre_neoplastic_0023"] <- "Patient.12"
Pal2021Group1@meta.data$patient[Pal2021Group1@meta.data$patient == "BRCA1_pre_neoplastic_0033"] <- "Patient.13"
Pal2021Group1@meta.data$patient[Pal2021Group1@meta.data$patient == "BRCA1_pre_neoplastic_0090"] <- "Patient.14"
Pal2021Group1@meta.data$patient[Pal2021Group1@meta.data$patient == "BRCA1_pre_neoplastic_0894"] <- "Patient.15"
Pal2021Group1@meta.data$patient[Pal2021Group1@meta.data$patient == "ER_positive__LN_0056"] <- "Patient.16"
Pal2021Group1@meta.data$patient[Pal2021Group1@meta.data$patient == "ER_positive_0001"] <- "Patient.17"
Pal2021Group1@meta.data$patient[Pal2021Group1@meta.data$patient == "ER_positive_0025"] <- "Patient.18"
Pal2021Group1@meta.data$patient[Pal2021Group1@meta.data$patient == "ER_positive_0029_7C"] <- "Patient.19"
Pal2021Group1@meta.data$patient[Pal2021Group1@meta.data$patient == "ER_positive_0029_9C"] <- "Patient.19"



df <- Pal2021Group2@meta.data
df$Cell_names <- rownames(df)

df <- match_df(Pal2021Group2_metadata, df, on = "Cell_names")
Pal2021Group2@meta.data$type <- df$Sample
Pal2021Group2@meta.data$patient <- df$Sample
Pal2021Group2@meta.data$annotation <- df$Annotation
Pal2021Group2@meta.data$disease <- df$Is_Cancer

Pal2021Group2@meta.data$disease[Pal2021Group2@meta.data$disease == 0] <- "reference"
Pal2021Group2@meta.data$disease[Pal2021Group2@meta.data$disease == 1] <- "cancer"

Pal2021Group2@meta.data$type[Pal2021Group2@meta.data$type == "ER_positive_0032"] <-  "ER"
Pal2021Group2@meta.data$type[Pal2021Group2@meta.data$type == "ER_positive_0040"] <- "ER"
Pal2021Group2@meta.data$type[Pal2021Group2@meta.data$type == "ER_positive_0042"] <- "ER"
Pal2021Group2@meta.data$type[Pal2021Group2@meta.data$type == "ER_positive_0043"] <-"ER"
Pal2021Group2@meta.data$type[Pal2021Group2@meta.data$type == "ER_positive_0056"] <- "ER"
Pal2021Group2@meta.data$type[Pal2021Group2@meta.data$type == "ER_positive_0064"] <- "ER"
Pal2021Group2@meta.data$type[Pal2021Group2@meta.data$type == "ER_positive_0068"] <- "ER"
Pal2021Group2@meta.data$type[Pal2021Group2@meta.data$type == "ER_positive_0114"] <- "ER"
Pal2021Group2@meta.data$type[Pal2021Group2@meta.data$type == "ER_positive_0125"] <- "ER"

Pal2021Group2@meta.data$patient[Pal2021Group2@meta.data$patient == "ER_positive_0032"] <- "Patient.20"
Pal2021Group2@meta.data$patient[Pal2021Group2@meta.data$patient == "ER_positive_0040"] <- "Patient.21"
Pal2021Group2@meta.data$patient[Pal2021Group2@meta.data$patient == "ER_positive_0042"] <- "Patient.22"
Pal2021Group2@meta.data$patient[Pal2021Group2@meta.data$patient == "ER_positive_0043"] <- "Patient.23"
Pal2021Group2@meta.data$patient[Pal2021Group2@meta.data$patient == "ER_positive_0056"] <- "Patient.16"
Pal2021Group2@meta.data$patient[Pal2021Group2@meta.data$patient == "ER_positive_0064"] <- "Patient.24"
Pal2021Group2@meta.data$patient[Pal2021Group2@meta.data$patient == "ER_positive_0068"] <- "Patient.25"
Pal2021Group2@meta.data$patient[Pal2021Group2@meta.data$patient == "ER_positive_0114"] <- "Patient.26"
Pal2021Group2@meta.data$patient[Pal2021Group2@meta.data$patient == "ER_positive_0125"] <- "Patient.27"



df <- Pal2021Group3@meta.data
df$Cell_names <- rownames(df)

df <- match_df(Pal2021Group3_metadata, df, on = "Cell_names")
Pal2021Group3@meta.data$type <- df$Sample
Pal2021Group3@meta.data$patient <- df$Sample
Pal2021Group3@meta.data$annotation <- df$Annotation
Pal2021Group3@meta.data$disease <- df$Is_Cancer

Pal2021Group3@meta.data$disease[Pal2021Group3@meta.data$disease == 0] <- "reference"
Pal2021Group3@meta.data$disease[Pal2021Group3@meta.data$disease == 1] <- "cancer"

Pal2021Group3@meta.data$type[Pal2021Group3@meta.data$type == "ER_positive_0151"] <-  "ER"
Pal2021Group3@meta.data$type[Pal2021Group3@meta.data$type == "ER_positive_0163"] <- "ER"
Pal2021Group3@meta.data$type[Pal2021Group3@meta.data$type == "ER_positive_0167"] <- "ER"
Pal2021Group3@meta.data$type[Pal2021Group3@meta.data$type == "ER_positive_0173"] <-"ER"
Pal2021Group3@meta.data$type[Pal2021Group3@meta.data$type == "ER_positive_0178"] <- "ER"
Pal2021Group3@meta.data$type[Pal2021Group3@meta.data$type == "ER_positive_0360"] <- "ER"
Pal2021Group3@meta.data$type[Pal2021Group3@meta.data$type == "ER_positive_LN_0040"] <- "ER"
Pal2021Group3@meta.data$type[Pal2021Group3@meta.data$type == "ER_positive_LN_0043"] <- "ER"
Pal2021Group3@meta.data$type[Pal2021Group3@meta.data$type == "ER_positive_LN_0064"] <- "ER"

Pal2021Group3@meta.data$patient[Pal2021Group3@meta.data$patient == "ER_positive_0151"] <- "Patient.28"
Pal2021Group3@meta.data$patient[Pal2021Group3@meta.data$patient == "ER_positive_0163"] <- "Patient.29"
Pal2021Group3@meta.data$patient[Pal2021Group3@meta.data$patient == "ER_positive_0167"] <- "Patient.30"
Pal2021Group3@meta.data$patient[Pal2021Group3@meta.data$patient == "ER_positive_0173"] <- "Patient.31"
Pal2021Group3@meta.data$patient[Pal2021Group3@meta.data$patient == "ER_positive_0178"] <- "Patient.32"
Pal2021Group3@meta.data$patient[Pal2021Group3@meta.data$patient == "ER_positive_0360"] <- "Patient.33"
Pal2021Group3@meta.data$patient[Pal2021Group3@meta.data$patient == "ER_positive_LN_0040"] <- "Patient.21"
Pal2021Group3@meta.data$patient[Pal2021Group3@meta.data$patient == "ER_positive_LN_0043"] <- "Patient.23"
Pal2021Group3@meta.data$patient[Pal2021Group3@meta.data$patient == "ER_positive_LN_0064"] <- "Patient.24"



df <- Pal2021Group4@meta.data
df$Cell_names <- rownames(df)

df <- match_df(Pal2021Group4_metadata, df, on = "Cell_names")
Pal2021Group4@meta.data$type <- df$Sample
Pal2021Group4@meta.data$patient <- df$Sample
Pal2021Group4@meta.data$annotation <- df$Annotation
Pal2021Group4@meta.data$disease <- df$Is_Cancer

Pal2021Group4@meta.data$disease[Pal2021Group4@meta.data$disease == 0] <- "reference"
Pal2021Group4@meta.data$disease[Pal2021Group4@meta.data$disease == 1] <- "cancer"

Pal2021Group4@meta.data$type[Pal2021Group4@meta.data$type == "ER_positive_LN_0068"] <-  "ER"
Pal2021Group4@meta.data$type[Pal2021Group4@meta.data$type == "ER_positive_LN_0167"] <- "ER"
Pal2021Group4@meta.data$type[Pal2021Group4@meta.data$type == "ER_positive_LN_0173"] <- "ER"
Pal2021Group4@meta.data$type[Pal2021Group4@meta.data$type == "HER2_0069"] <-"HER"
Pal2021Group4@meta.data$type[Pal2021Group4@meta.data$type == "HER2_0161"] <- "HER"
Pal2021Group4@meta.data$type[Pal2021Group4@meta.data$type == "HER2_0176"] <- "HER"
Pal2021Group4@meta.data$type[Pal2021Group4@meta.data$type == "HER2_0308"] <- "HER"
Pal2021Group4@meta.data$type[Pal2021Group4@meta.data$type == "HER2_0331"] <- "HER"
Pal2021Group4@meta.data$type[Pal2021Group4@meta.data$type == "HER2_0337"] <- "HER"

Pal2021Group4@meta.data$patient[Pal2021Group4@meta.data$patient == "ER_positive_LN_0068"] <- "Patient.25"
Pal2021Group4@meta.data$patient[Pal2021Group4@meta.data$patient == "ER_positive_LN_0167"] <- "Patient.30"
Pal2021Group4@meta.data$patient[Pal2021Group4@meta.data$patient == "ER_positive_LN_0173"] <- "Patient.31"
Pal2021Group4@meta.data$patient[Pal2021Group4@meta.data$patient == "HER2_0069"] <- "Patient.34"
Pal2021Group4@meta.data$patient[Pal2021Group4@meta.data$patient == "HER2_0161"] <- "Patient.35"
Pal2021Group4@meta.data$patient[Pal2021Group4@meta.data$patient == "HER2_0176"] <- "Patient.36"
Pal2021Group4@meta.data$patient[Pal2021Group4@meta.data$patient == "HER2_0308"] <- "Patient.37"
Pal2021Group4@meta.data$patient[Pal2021Group4@meta.data$patient == "HER2_0331"] <- "Patient.38"
Pal2021Group4@meta.data$patient[Pal2021Group4@meta.data$patient == "HER2_0337"] <- "Patient.39"



df <- Pal2021Group5@meta.data
df$Cell_names <- rownames(df)

df <- match_df(Pal2021Group5_metadata, df, on = "Cell_names")
Pal2021Group5@meta.data$type <- df$Sample
Pal2021Group5@meta.data$patient <- df$Sample
Pal2021Group5@meta.data$annotation <- df$Annotation
Pal2021Group5@meta.data$disease <- df$Is_Cancer

Pal2021Group5@meta.data$disease[Pal2021Group5@meta.data$disease == 0] <- "reference"
Pal2021Group5@meta.data$disease[Pal2021Group5@meta.data$disease == 1] <- "cancer"

Pal2021Group5@meta.data$type[Pal2021Group5@meta.data$type == "PR_positive_0319"] <-  "PR"
Pal2021Group5@meta.data$type[Pal2021Group5@meta.data$type == "Triple_negative_0106"] <- "TNBC"
Pal2021Group5@meta.data$type[Pal2021Group5@meta.data$type == "Triple_negative_0114"] <- "TNBC"
Pal2021Group5@meta.data$type[Pal2021Group5@meta.data$type == "Triple_negative_0126"] <-"TNBC"
Pal2021Group5@meta.data$type[Pal2021Group5@meta.data$type == "Triple_negative_0135"] <- "TNBC"
Pal2021Group5@meta.data$type[Pal2021Group5@meta.data$type == "Triple_negative_BRCA1_0131"] <- "TNBC"
Pal2021Group5@meta.data$type[Pal2021Group5@meta.data$type == "Triple_negative_BRCA1_0177"] <- "TNBC"
Pal2021Group5@meta.data$type[Pal2021Group5@meta.data$type == "Triple_negative_BRCA1_0554"] <- "TNBC"
Pal2021Group5@meta.data$type[Pal2021Group5@meta.data$type == "Triple_negative_BRCA1_4031"] <- "TNBC"

Pal2021Group5@meta.data$patient[Pal2021Group5@meta.data$patient == "PR_positive_0319"] <- "Patient.40"
Pal2021Group5@meta.data$patient[Pal2021Group5@meta.data$patient == "Triple_negative_0106"] <- "Patient.41"
Pal2021Group5@meta.data$patient[Pal2021Group5@meta.data$patient == "Triple_negative_0114"] <- "Patient.26"
Pal2021Group5@meta.data$patient[Pal2021Group5@meta.data$patient == "Triple_negative_0126"] <- "Patient.42"
Pal2021Group5@meta.data$patient[Pal2021Group5@meta.data$patient == "Triple_negative_0135"] <- "Patient.43"
Pal2021Group5@meta.data$patient[Pal2021Group5@meta.data$patient == "Triple_negative_BRCA1_0131"] <- "Patient.44"
Pal2021Group5@meta.data$patient[Pal2021Group5@meta.data$patient == "Triple_negative_BRCA1_0177"] <- "Patient.45"
Pal2021Group5@meta.data$patient[Pal2021Group5@meta.data$patient == "Triple_negative_BRCA1_0554"] <- "Patient.46"
Pal2021Group5@meta.data$patient[Pal2021Group5@meta.data$patient == "Triple_negative_BRCA1_4031"] <- "Patient.47"

# Merge datasets----
seurat.merged <- merge(x = Gao2021, y = c(Pal2021Group1, Pal2021Group2, Pal2021Group3, Pal2021Group4, Pal2021Group5, Qian2020), add.cell.ids = dirs, project = "breast")

seurat.merged[["percent.mt"]] <- PercentageFeatureSet(seurat.merged, pattern = "^MT-")
Idents(seurat.merged) <- seurat.merged@meta.data$orig.ident

pdf("figures/low_quality_cells.pdf")
FeatureScatter(seurat.merged, feature1 = "nFeature_RNA", feature2 = "percent.mt", raster = F) + theme(aspect.ratio = 1)
dev.off()

saveRDS(seurat.merged, file = "data/seurat.merged.rds")
# seurat.merged <- readRDS("data/seurat.merged.rds"))

# Filter-out low quality cells----
pdf("figures/QCmetrics_before_filtering.pdf")
VlnPlot(seurat.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F) + theme(legend.position = "none", axis.title.x = element_blank())
dev.off()

## Estimate cells-to-remove----
library(miQC)
library(flexmix)
library(splines)
library(scater)

x <- as.SingleCellExperiment(seurat.merged)

mt_genes <- grepl("^MT-",  rownames(x))
feature_ctrls <- list(mito = rownames(x)[mt_genes])

x <- addPerCellQC(x, subsets = feature_ctrls)
plotMetrics(x)
model <- mixtureModel(x, model_type = "spline")
plotModel(x, model)
metrics <- as.data.frame(colData(x))

intercept1 <- parameters(model, component = 1)[1]
intercept2 <- parameters(model, component = 2)[1]

if (intercept1 > intercept2) {
  compromised_dist <- 1
  intact_dist <- 2
} else {
  compromised_dist <- 2
  intact_dist <- 1
}

post <- posterior(model); posterior_cutoff = 0.75
prob_compromised <- post[, compromised_dist]
keep <- prob_compromised <= posterior_cutoff

metrics <- cbind(metrics, prob_compromised = prob_compromised, keep = keep)

predictions <- fitted(model)[, intact_dist]
metrics$intact_prediction <- predictions
metrics[metrics$subsets_mito_percent < metrics$intact_prediction, ]$keep <- TRUE

print(low_quality_percentage <- min(metrics[!metrics$keep, ]$subsets_mito_percent))

pdf("figures/low_quality_threshold.pdf")
print(plotFiltering(x, model, detected = "detected", subsets_mito_percent = "subset_mito_percent") + geom_hline(yintercept=low_quality_percentage, linetype="dashed", color = "red") + theme_minimal())
dev.off()

rm(x)

## Subset cells----
seurat.merged.filtered <- subset(seurat.merged, subset = nFeature_RNA < 7000 & nCount_RNA < 40000 & percent.mt < low_quality_percentage & annotation != "Filtered-out")
rm(seurat.merged)

pdf("figures/QCmetrics_after_filtering.pdf")
VlnPlot(seurat.merged.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F) + theme(legend.position = "none", axis.title.x = element_blank())
dev.off()

saveRDS(seurat.merged.filtered, file = "data/seurat.merged.filtered.rds")
# seurat.merged.filtered <- readRDS("data/seurat.merged.filtered.rds")

setwd(projectDir)
#----session info----
sessionInfo()
