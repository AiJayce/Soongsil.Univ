# Immune Evasion in Colorectal Cancer

## (1) Introduction

Immune checkpoint inhibitors (ICIs) have shown limited efficacy in MSS colorectal cancer, highlighting an urgent need to uncover the mechanisms of immune evasion. This project identified a TP53 loss-driven CEBPB signaling axis that promotes immune evasion by inducing CTLA4 upregulation in T cells, providing mechanistic insights into immunotherapy resistance and suggesting potential targeted therapeutic strategies for colorectal cancer.

**Keyword:** scRNA seq, TCR clonal expansion, Bulk RNA seq, Spatial transcriptomics, GMM, TCGA, Cancer genome

![Single-cell immune-profiling workflow](assets/01-scrna-pipeline.png)

## (2) ChIP-seq analysis and (3) miRNA analysis

To validate the proposed CEBPB-CTLA4 regulatory relationship, I performed p53 target ChIP-seq analysis, which showed stronger p53 binding in TP53-mutant SW480 cells and supported the scRNA-seq findings.

Cell line - SW480 (p53 Mutation)  
Cell line - HCT116 (p53 Wildtype)

To identify factors that may disrupt the proposed CEBPB-CTLA4 regulatory axis, I analyzed miRNAs potentially involved in CEBPB-mediated CTLA4 regulation under TP53-mutant conditions. Differential expression analysis identified TP53 mutation-specific miRNAs, highlighting candidate post-transcriptional regulators of the transcription factor gene regulatory network (GRN).

**Keyword:** ChIP seq, Mutation analysis, Peak calling, Motif analysis (Homer), Public data, Cancer genome, Validation, miRNA, edgeR

![ChIP-seq comparison of TP53-mutant and TP53-wild-type cell lines](assets/02-chip-seq.png)

## (4) Spatial Spot-level coexpression

I investigated whether the regulatory relationship identified in the scRNA-seq analysis could be reproduced in spatial transcriptomics data. Pseudo-spots were defined by aggregating each spot with its 1-layer neighboring spots to represent the local spatial microenvironment. A sliding window algorithm was then applied to calculate the summed expression values for each pseudo-spot across the tissue, and Pearson correlation analysis demonstrated a significant spatial correlation between CEBPB and CTLA4, providing additional evidence supporting the robustness of the scRNA-seq findings.

**Keyword:** Spatial transcriptomics, 10X Visium, Sliding window algorithm, Public data, Validation, Co-expression

![Spatial pseudo-spot coexpression analysis of CEBPB and CTLA4](assets/04-spatial-coexpression.png)


