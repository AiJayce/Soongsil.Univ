<img width="2153" height="633" alt="Immune Evasion in Colorectal Cancer-2" src="https://github.com/user-attachments/assets/c93b3148-d254-45ed-8fae-e78c80aa3a0f" /># Immune Evasion in Colorectal Cancer

## (1) Introduction

<img width="1098" height="964" alt="41568_2024_715_Fig2_HTML" src="https://github.com/user-attachments/assets/b67e4073-581e-4d81-9ec4-8d7cb34ba426" />

Immune checkpoint inhibitors (ICIs) have shown limited efficacy in MSS colorectal cancer, highlighting an urgent need to uncover the mechanisms of immune evasion. This project identified a TP53 loss-driven CEBPB signaling axis that promotes immune evasion by inducing CTLA4 upregulation in T cells, providing mechanistic insights into immunotherapy resistance and suggesting potential targeted therapeutic strategies for colorectal cancer.

## (2) Basic Analysis Pipeline

<img width="3224" height="1287" alt="Immune Evasion in Colorectal Cancer-1" src="https://github.com/user-attachments/assets/6ab430ae-afc5-4ee9-adba-193521e936f4" />

The basic analysis pipeline was performed using the Seurat pipeline in R. A positive correlation was observed between high CEBPB expression in epithelial cells and increased CTLA4 expression in the T-cell population. This finding suggests that elevated CEBPB transcription factor activity may contribute to CD4 T-cell exhaustion, potentially enabling immune checkpoint inhibitor (ICI) evasion.


## (3) ChIP-seq analysis

<img width="2153" height="633" alt="Immune Evasion in Colorectal Cancer-2" src="https://github.com/user-attachments/assets/a84fc2be-e477-4def-8b9a-8ca42aca0c9b" />

To validate the proposed CEBPB–CTLA4 regulatory relationship, I performed p53 target ChIP-seq analysis, which showed stronger p53 binding in TP53-mutant SW480 cells and supported the scRNA-seq findings.

## (4) miRNA analysis

<img width="2068" height="1247" alt="Immune Evasion in Colorectal Cancer-3" src="https://github.com/user-attachments/assets/3025ddc4-e35c-4b66-8ea0-a9778750ad35" />

To identify factors that may disrupt the proposed CEBPB–CTLA4 regulatory axis, I analyzed miRNAs potentially involved in CEBPB-mediated CTLA4 regulation under TP53-mutant conditions. Differential expression analysis identified TP53 mutation-specific miRNAs, highlighting candidate post-transcriptional regulators of the transcription factor gene regulatory network (GRN).

## (4) Spatial Spot-level coexpression

<img width="3224" height="1184" alt="Immune Evasion in Colorectal Cancer-4" src="https://github.com/user-attachments/assets/4fea96c7-d5a5-4467-a054-dfef030d1775" />

I investigated whether the regulatory relationship identified in the scRNA-seq analysis could be reproduced in spatial transcriptomics data. Pseudo-spots were defined by aggregating each spot with its 1-layer neighboring spots to represent the local spatial microenvironment. A sliding window algorithm was then applied to calculate the summed expression values for each pseudo-spot across the tissue, and Pearson correlation analysis demonstrated a significant spatial correlation between CEBPB and CTLA4, providing additional evidence supporting the robustness of the scRNA-seq findings.

## Discussion

### Data and analysis context

The project integrated colorectal cancer scRNA-seq, TCR clonotype data, bulk RNA-seq, cancer-genome information, p53 ChIP-seq, miRNA sequencing, and 10x Visium spatial transcriptomics. These complementary data types were used to move from immune-cell profiling to regulatory validation and spatial confirmation.

### Interpretation

The single-cell analysis nominated a TP53 loss-associated CEBPB-CTLA4 axis linked to T-cell immune evasion. ChIP-seq supported the TP53-dependent regulatory context, miRNA analysis identified candidate post-transcriptional regulators, and spatial pseudo-spot analysis reproduced positive CEBPB-CTLA4 coexpression across samples.

### Considerations

The analyses provide convergent evidence but do not by themselves establish that the CEBPB-CTLA4 relationship is causal. Differences among cell lines, public datasets, sequencing platforms, and spatial resolution may affect reproducibility. Perturbation of TP53, CEBPB, or CTLA4 in matched experimental systems would be needed to test the proposed mechanism directly.


