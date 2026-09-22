# Implantation Failure in Uterine Aging

## Overview

This project investigates transcriptomic mechanisms underlying implantation failure in the aged uterus. It identifies an aging-associated transition from a receptive endometrial network to a non-receptive, mucin-overexpression network driven by senescent immune-cell interactions.

![Implantation failure and senescent CD8 T-cell cause inference](assets/01-implantation-problem.png)

## Workflow

### 1. Define receptive and non-receptive endometrial networks

Gene regulatory network modules were built from luminal epithelial scRNA-seq data. Hierarchical clustering separated receptive and non-receptive networks, with mucin genes enriched in the non-receptive modules.

![Luminal epithelial receptive and non-receptive network modules](assets/02-network-modules.png)

### 2. Classify network states and prioritize transition factors

The same ensemble Random Forest classifier was used to annotate regulatory modules and assess feature importance. Non-receptive modules were predominantly enriched in aged uterus samples.

![Ensemble Random Forest classification of regulatory modules](assets/03-random-forest-classifier.png)

TENET was then applied to a trajectory from receptive to non-receptive epithelial states. Transcription factors regulating the transition were prioritized with graph-theoretic criteria, identifying Mxd1 as a candidate regulator of the mucin-associated program.

![TENET trajectory, network inference, and Mxd1 prioritization](assets/04-tenet-transition-factors.png)

### 3. Infer immune-to-epithelial signaling

NicheNet analysis identified CD8 T-cell-mediated signals associated with transcription-factor upregulation in luminal epithelial receiver cells. These interactions support a model in which senescent immune signaling promotes the non-receptive state.

![NicheNet inference of CD8 T-cell to luminal epithelial signaling](assets/05-nichenet-cci.png)

### 4. Validate across spatial and epigenomic data

Spatial transcriptomics identified aged samples enriched for receiver cells with elevated transition-factor and mucin expression. ATAC-seq peak calling and HOMER motif enrichment supported accessible regulatory regions and transcription-factor binding at mucin loci.

![Spatial validation of senescent-cell interactions and Muc4 expression](assets/06-spatial-validation.png)

![ATAC-seq peak calling and HOMER motif validation at the Muc4 locus](assets/07-atac-validation.png)

## Methods

- scRNA-seq and spatial transcriptomics
- Gene regulatory network construction and scWGCNA
- Ensemble Random Forest classification and feature importance
- Pseudotime and TENET network-transition inference
- NicheNet cell-cell communication analysis
- ATAC-seq peak calling and HOMER motif enrichment

## Key Finding

Senescent CD8 T-cell signaling is associated with an Mxd1-linked transition to a mucin-overexpression, non-receptive luminal epithelial network, providing a potential mechanism for implantation failure in uterine aging.

