library(Seurat)
# devtools::install_github("miccec/yaGST")
# devtools::install_github("AntonioDeFalco/SCEVAN")
library(SCEVAN)
library(dplyr)

project_dir <- "~/BreastCancerG0arrest"
setwd(paste0(project_dir, "/06_TumourClones"))

# Load the data
integrated <- readRDS(paste0(project_dir, "/02_G0arrestInMalignantCells/data/integrated_with_quiescence.rds"))
integrated <- DietSeurat(integrated, assays = "RNA")

# Subset the data to only include malignant cells
malignant <- subset(
  integrated,
  subset = type %in% c("TNBC", "ER", "PR", "HER", "DCIS", "neoplasm") &
    malignancy %in% "malignant"
)

rm(integrated)

# Run SCEVAN for ER G0 subtype-------------------------------------------------
er_g0 <- subset(
  malignant,
  subset = QuiescenceStatus %in% "Quiescent" & type %in% "ER"
)

dir.create("ER_G0_arrested")
setwd("ER_G0_arrested")

# Get the counts
count_mtx <- GetAssayData(er_g0, layer = "counts", assay = "RNA")

set.seed(42)
if (ncol(count_mtx) > 3000) {
  count_mtx <- count_mtx[, sample(ncol(count_mtx), 3000)]
}

gc(full = TRUE)
# Run SCEVAN
results <- SCEVAN::pipelineCNA(
  count_mtx,
  sample = "ER_G0_arrested",
  par_cores = 4,
  SUBCLONES = TRUE,
  plotTree = TRUE
)

setwd("..")

# Run SCEVAN for TNBC G0 subtype-----------------------------------------------
tnbc_g0 <- subset(
  malignant,
  subset = QuiescenceStatus %in% "Quiescent" & type %in% "TNBC"
)

dir.create("TNBC_G0_arrested")
setwd("TNBC_G0_arrested")

# Get the counts
count_mtx <- GetAssayData(tnbc_g0, layer = "counts", assay = "RNA")

gc(full = TRUE)
# Run SCEVAN
results <- SCEVAN::pipelineCNA(
  count_mtx,
  sample = "TNBC_G0_arrested",
  par_cores = 20,
  SUBCLONES = TRUE,
  plotTree = TRUE
)

setwd("..")

# Run SCEVAN for ER Fast cycling subtype---------------------------------------
er_fast <- subset(
  malignant,
  subset = QuiescenceStatus %in% "Proliferating" & type %in% "ER"
)

dir.create("ER_Fast_cycling")
setwd("ER_Fast_cycling")

# Get the counts
count_mtx <- GetAssayData(er_fast, layer = "counts", assay = "RNA")

gc(full = TRUE)
# Run SCEVAN
results <- SCEVAN::pipelineCNA(
  count_mtx,
  sample = "ER_Fast_cycling",
  par_cores = 20,
  SUBCLONES = TRUE,
  plotTree = TRUE
)

setwd("..")

# Run SCEVAN for TNBC Fast cycling subtype-------------------------------------
tnbc_fast <- subset(
  malignant,
  subset = QuiescenceStatus %in% "Proliferating" & type %in% "TNBC"
)

dir.create("TNBC_Fast_cycling")
setwd("TNBC_Fast_cycling")

# Get the counts
count_mtx <- GetAssayData(tnbc_fast, layer = "counts", assay = "RNA")

set.seed(42)
if (ncol(count_mtx) > 3000) {
  count_mtx <- count_mtx[, sample(ncol(count_mtx), 3000)]
}

gc(full = TRUE)
# Run SCEVAN
results <- SCEVAN::pipelineCNA(
  count_mtx,
  sample = "TNBC_Fast_cycling",
  par_cores = 20,
  SUBCLONES = TRUE,
  plotTree = TRUE
)

setwd("..")