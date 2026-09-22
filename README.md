# Predicting Hippocampal Regional Vulnerability in Alzheimer's Disease

## Overview

The hippocampus has distinct anterior and posterior functions. This project investigates which hippocampal region is more vulnerable to Alzheimer's disease (AD)-related functional loss by integrating spatial transcriptomics, cell-type annotation, plaque-associated latent factors, and cell-cell interaction analysis.

## Workflow

### 1. Spatial annotation and regional comparison

Spatial cells were annotated in two complementary ways: high-confidence reference label transfer from the Allen Mouse Brain Atlas hippocampal formation scRNA-seq reference, and manual Common Coordinate Framework (CCF) region mapping. The mapped regions were then used to compare cell-type distributions between AD and wild-type samples.

![Posterior and anterior hippocampal regional functions](assets/01-regional-functions.png)

![Reference label transfer and CCF mapping](assets/02-label-transfer-ccf.png)

Microglia were elevated across mapped AD regions, supporting a disease-associated relationship between microglia and regional vulnerability.

![Cell-type distributions and regional proportions](assets/03-regional-cell-types.png)

### 2. Plaque-associated latent-factor modeling

Plaque protein size and spatial coordinates were combined into a plaque score. A conditional variational autoencoder represented spatial gene-expression programs in a latent space, and ridge-regularized linear regression related plaque scores to latent dimensions.

![Conditional VAE and plaque-score modeling workflow](assets/04-cvae-plaque-model.png)

Microglia showed the strongest plaque association. Latent dimension z15, which was negatively associated with plaque score, was decoded to genes enriched for lipid metabolism and chronic inflammatory pathways.

![Latent-factor decoding and pathway enrichment](assets/05-latent-decoding.png)

### 3. Regional vulnerability and cell-cell interactions

Decoded z15 genes were prioritized by expression density. APOE and TREM2 were retained despite sparse expression because of their established AD relevance; both were more highly expressed in the posterior hippocampus, indicating increased regional vulnerability.

![APOE and TREM2 selection and regional expression](assets/06-gene-selection.png)

Nearest-cell linear regression identified a potential association between oligodendrocyte density and microglial TREM2/APOE expression. NicheNet and pseudobulk analyses were then used to reproduce ligand-receptor signals regulating the microglial program in scRNA-seq data.

![Nearest-cell regression for oligodendrocyte-microglia associations](assets/07-neighbor-inference.png)

![NicheNet and pseudobulk reproducibility analysis](assets/08-cci-reproducibility.png)

## Methods

- Spatial transcriptomics and MERFISH data analysis
- Reference label transfer with scVI and CCF mapping
- Conditional variational autoencoder and ridge regression
- Latent-space gene-set decoding and pathway enrichment
- Nearest-cell regression, NicheNet, and pseudobulk validation

## Key Finding

The posterior hippocampus showed higher APOE and TREM2 expression in microglia, consistent with greater AD-related regional vulnerability and a plaque-associated chronic inflammatory program.

