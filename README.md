# Predicting Hippocampal Regional Vulnerability in Alzheimer's Disease


## (1) Introduction

<img width="1090" height="1155" alt="Predicting Hippocampal Regional Vulnerability in AD-1" src="https://github.com/user-attachments/assets/f59b08bc-5586-4e62-b4be-7e6adfef0b87" />

The hippocampus exhibits distinct functions across its anterior and posterior regions. This project inferred which region is more vulnerable to AD-related functional loss.

## (2) Annotation 

<img width="1048" height="1154" alt="Predicting Hippocampal Regional Vulnerability in AD-2" src="https://github.com/user-attachments/assets/88ee4c8b-fcb0-4096-9907-82af25204c97" />

Annotation was performed at two levels: deep learning based scRNA reference label transfer with high-confidence scores, and manual CCF region mapping.

@code : https://github.com/AiJayce/Soongsil.Univ/blob/Alzheimer-Vulnerability/MERSCOPE.ipynb

## (3) Compare region specific

<img width="1204" height="1154" alt="Predicting Hippocampal Regional Vulnerability in AD-3" src="https://github.com/user-attachments/assets/e0021f29-e78a-4887-9b5e-8ae4f6688711" />

Across mapped regions, Microglia showed significantly higher populations in AD samples across all regions, suggesting a disease-associated relationship between Microglia and AD.

## (4) VAE modeling 

<img width="1499" height="1155" alt="Predicting Hippocampal Regional Vulnerability in AD-4" src="https://github.com/user-attachments/assets/cfbaaaa0-7d15-49dc-98df-77028bb5e5ed" />

Metadata were used to quantify plaque protein size and spatial coordinates to derive a plaque score. Ridge-regularized multiple linear regression was performed between the plaque score and VAE latent dimensions. Based on this concept, latent dimensions with negative coefficients were interpreted as representing genes associated with both plaque clearance and active inflammatory responses.

## (5) Latent decoding to gene sets

<img width="1768" height="1155" alt="Predicting Hippocampal Regional Vulnerability in AD-5" src="https://github.com/user-attachments/assets/bcaf1bf1-ffad-4024-b43b-e97794c4d779" />

Microglia showed the largest absolute latent coefficients across cell types. Since persistent chronic inflammation during plaque clearance may accelerate Alzheimer’s disease pathology, I focused on the plaque-downregulating latent z15. Decoding z15 revealed enrichment of lipid metabolism and chronic inflammation pathways, suggesting its role in Alzheimer’s disease pathology.


## (6) Select significant genes

<img width="1090" height="1155" alt="Predicting Hippocampal Regional Vulnerability in AD-6" src="https://github.com/user-attachments/assets/d1a938ce-f56a-43a9-8a48-6d46f63e3b5b" />

Genes decoded from z15 were selected based on expression density. APOE and TREM2 were retained despite sparse expression and confirmed as AD-associated genes in the literature. Both showed significantly higher expression in pHIP, suggesting increased AD vulnerability in this region.

##  (7) Causal Inference from CCI 

<img width="1053" height="1151" alt="Predicting Hippocampal Regional Vulnerability in AD-7" src="https://github.com/user-attachments/assets/e52d7c21-bb76-4e58-a337-66a22d8cc872" />

Linear regression identified a potential causal relationship between higher Oligodendrocyte density and increased TREM2 and APOE expression in Microglia.

## (8) CCI Reproducibility in scRNA

<img width="1111" height="1155" alt="Predicting Hippocampal Regional Vulnerability in AD-8" src="https://github.com/user-attachments/assets/377a86c7-89b5-4a8a-8803-a33bef169524" />

At the scRNA-seq level, NicheNet was used to infer ligand–receptor interactions. Pseudobulk analysis confirmed that significant L/R responses regulate TREM2–APOE expression within Microglia, supporting reproducibility.


## Discussion

### Data and analysis context

This project integrated MERFISH spatial transcriptomics with a normal hippocampal formation scRNA-seq reference from the Allen Mouse Brain Atlas. Reference label transfer and manual CCF mapping were used to annotate spatial cells and distinguish posterior and anterior hippocampal regions. Plaque size and spatial coordinates were incorporated as metadata for plaque-score modeling.

### Interpretation

The regional comparison suggested increased microglial abundance in AD samples. Conditional VAE latent dimensions captured plaque-associated expression programs, and z15 was interpreted as a plaque-downregulating, inflammation-associated factor. Decoding z15 highlighted lipid metabolism and chronic inflammation, while APOE and TREM2 showed higher expression in the posterior hippocampus. Nearest-cell regression and NicheNet provided complementary evidence for oligodendrocyte-microglia interactions associated with TREM2-APOE regulation.

### Considerations

The regional vulnerability interpretation depends on spatial annotation accuracy, reference-transfer confidence, plaque-score construction, and the observational nature of cell-cell interaction inference. The results support a regional association and a biologically plausible mechanism, but causal validation would require perturbation experiments and independent AD datasets.

