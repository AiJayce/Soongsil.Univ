# Immune Evasion in Colorectal Cancer

## Overview

Immune checkpoint inhibitors have limited efficacy in microsatellite-stable colorectal cancer. This project investigates mechanisms of immune evasion and identifies a TP53 loss-associated CEBPB signaling axis that is associated with CTLA4 upregulation in T cells.

## Workflow

### 1. Single-cell immune profiling

Colorectal cancer single-cell RNA-seq data were processed through quality control, integration, broad and sub-cell-type annotation, and downstream immune analyses. The workflow included copy-number variation, Gaussian mixture modeling, differential expression, correlation analysis, functional characterization, cell-cell interaction analysis, and T-cell receptor clonotype and clonal-expansion analyses.

![Single-cell immune-profiling workflow](assets/01-scrna-pipeline.png)

### 2. ChIP-seq validation of the regulatory relationship

p53-target ChIP-seq tracks were compared between TP53-mutant SW480 and TP53-wild-type HCT116 cells. Stronger p53 binding in the TP53-mutant context supported the proposed CEBPB-CTLA4 regulatory relationship from the single-cell analysis.

![ChIP-seq comparison of TP53-mutant and TP53-wild-type cell lines](assets/02-chip-seq.png)

### 3. miRNA analysis

Candidate miRNAs that could disrupt CEBPB-mediated CTLA4 regulation under TP53-mutant conditions were identified. The workflow included read preprocessing, differential-expression analysis across conditions, and target-miRNA selection with miRDB. This highlighted mutation-specific candidate post-transcriptional regulators of the gene regulatory network.

![miRNA preprocessing, differential expression, and target selection](assets/03-mirna-analysis.png)

### 4. Spatial coexpression validation

The CEBPB-CTLA4 relationship was evaluated in 10x Visium spatial transcriptomics data. Each pseudo-spot combined a spot with its one-layer neighbors to represent its local microenvironment. Pearson correlation analysis showed significant coexpression across four samples (PCC = 0.28-0.51; P < 0.001).

![Spatial pseudo-spot coexpression analysis of CEBPB and CTLA4](assets/04-spatial-coexpression.png)

## Methods

- scRNA-seq, TCR clonotyping, and clonal-expansion analysis
- Bulk RNA-seq and cancer-genome analysis
- p53 ChIP-seq, peak calling, and HOMER motif analysis
- miRNA differential expression with edgeR and miRDB prioritization
- 10x Visium spatial transcriptomics and sliding-window coexpression analysis

## Key Finding

The analyses support a TP53 loss-associated CEBPB-CTLA4 axis as a candidate mechanism of immune evasion and immunotherapy resistance in colorectal cancer.

