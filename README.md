# Predicting Hippocampal Regional Vulnerability in Alzheimer's Disease


<img width="1090" height="1155" alt="Predicting Hippocampal Regional Vulnerability in AD-1" src="https://github.com/user-attachments/assets/2cf32d47-b731-4d70-b512-e799c5c19b7c" />

## (1) Introduction

The hippocampus exhibits distinct functions across its anterior and posterior regions. This project inferred which region is more vulnerable to AD-related functional loss.

## (2) Annotation and (3) Compare region specific

Annotation was performed at two levels: deep learning based scRNA reference label transfer with high-confidence scores, and manual CCF region mapping.

Across mapped regions, Microglia showed significantly higher populations in AD samples across all regions, suggesting a disease-associated relationship between Microglia and AD.

**Keyword:** Spatial transcriptomics, Merscope, Label transfer (SCVI), CCF mapping, Public data

![Introduction, annotation, and region-specific comparison](assets/02-label-transfer-ccf.png)

## (4) VAE modeling and (5) Latent decoding to gene sets

Metadata were used to quantify plaque protein size and spatial coordinates to derive a plaque score. Ridge-regularized multiple linear regression was performed between the plaque score and VAE latent dimensions. Based on this concept, latent dimensions with negative coefficients were interpreted as representing genes associated with both plaque clearance and active inflammatory responses.

Microglia showed the largest absolute latent coefficients across cell types. Since persistent chronic inflammation during plaque clearance may accelerate Alzheimer's disease pathology, I focused on the plaque-downregulating latent z15. Decoding z15 revealed enrichment of lipid metabolism and chronic inflammation pathways, suggesting its role in Alzheimer's disease pathology.

**Keyword:** Variational Auto Encoder, Plaque score, Latent coefficient function, Ridge penalty, Latent decode

![VAE modeling and latent decoding](assets/04-cvae-plaque-model.png)

## (6) Select significant genes, (7) Causal Inference from CCI, and (8) CCI Reproducibility in scRNA

Genes decoded from z15 were selected based on expression density. APOE and TREM2 were retained despite sparse expression and confirmed as AD-associated genes in the literature. Both showed significantly higher expression in pHIP, suggesting increased AD vulnerability in this region.

Linear regression identified a potential causal relationship between higher Oligodendrocyte density and increased TREM2 and APOE expression in Microglia. At the scRNA-seq level, NicheNet was used to infer ligand-receptor interactions. Pseudobulk analysis confirmed that significant L/R responses regulate TREM2-APOE expression within Microglia, supporting reproducibility.

**Keyword:** TREM2, APOE, Alzheimer, Gaussian density, Cell-cell interaction, Linear regression, Nearest cell equation, Nichenet, Pseudobulk

![Significant genes, CCI causal inference, and reproducibility](assets/08-cci-reproducibility.png)
<img width="1111" height="1155" alt="Predicting Hippocampal Regional Vulnerability in AD-8" src="https://github.com/user-attachments/assets/cd77d161-8768-4ce0-b112-2a1adff44ba0" />

