# G0 breast cancer atlas

### Author: Cenk Celik, UCL Genetics Institute
This repository contains scripts for the evaluation of intrinsic and extrinsic regulation of G0 cell cycle arrest in breast cancer using single-cell and spatial transcriptomics data.

![Graphical abstract](img/graphical_abstract.jpg)

#### Integrating public data

Scripts in [00_IntegratingDatasets](00_IntegratingDatasets) describe how to [download and preprocess](00_IntegratingDatasets/01_CreateObjectsFromPublicData.R) publicly available data, integrate them using [`SCTransform` v2](00_IntegratingDatasets/02_IntegrateDatasetswithSCTransform.R), annotate cell types using marker gene expressions from the [`PanglaoDB`](https://panglaodb.se). The full cell type annotation database can be download from the [link](https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz). Using the full cell type databae, generate the input for scAnnotate [here](00_IntegratingDatasets/03_a_CreateAnnotationDatabase.R). Make sure to also download helper function [1](00_IntegratingDatasets/03_a1_gene_sets_prepare.r) and [2](00_IntegratingDatasets/03_a1_gene_sets_prepare.R) during cell type [annotation](00_IntegratingDatasets/03_b_AnnotateIntegratedDataset.R). Since, the `PanglaoDB` do not have cancer associated fibroblast markers, we checked for this separately in [04_IdentifyingCancerAssociatedFibroblasts](00_IntegratingDatasets/04_IdentifyingCancerAssociatedFibroblasts.R).

#### Inferring copy number alterations

We leveraged [`infercnv` package](https://github.com/broadinstitute/infercnv) to [infer](01_InferCNV/01_InferCopyNumberVariations.R) copy number variations in epithelial signle cells.

#### G0 arrest scoring in malignant cells

First, apply G0 arrest scoring using the gene sets [downregulated_common.RData](02_G0arrestInMalignantCells/data/downregulated_common.RData) and [upregulated_common.RData](02_G0arrestInMalignantCells/data/upregulated_common.RData).
> **Note:** if the data has ENSEMBL ID rather than HGNC symbols, visit the original study by [Wiecek *et al.* 2023](https://github.com/secrierlab/CancerG0Arrest) for the appropriate gene sets.
To evaluate pathway enrichment in malignant cells based on their cell cycle status, use [02_EnrichmentAnalysis](02_G0arrestInMalignantCells/02_EnrichmentAnalysis.R). One could use [03_DifferentialAbundanceTesting](02_G0arrestInMalignantCells/03_DifferentialAbundanceTesting.R), [04_DifferentialExpression](02_G0arrestInMalignantCells/04_DifferentialExpression.R) to compute differential abundance and gene expresssion per cell cycle state, and [05_GeneExpression](02_G0arrestInMalignantCells/05_GeneExpression.R) to plot genes of interest.
