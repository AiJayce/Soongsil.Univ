# Implantation Failure in Uterine Aging

## (1) Introduction

I elucidated the transcriptomic mechanisms underlying implantation failure in the aged uterus by demonstrating how aging drives a transition toward a non-receptive endometrial network. Because successful embryo implantation, including after IVF, requires low mucin expression, I discovered that senescent immune cell-mediated interactions promote a transcriptomic transition toward a mucin-overexpression network, thereby remodeling the uterine microenvironment into an implantation-unfavorable state.

**Keyword:** scRNA seq, Scanpy, HI-C, ML/DL, Spatial transcriptomics, ATAC, Epigenomics, Public data

![Implantation failure and senescent CD8 T-cell cause inference](assets/01-implantation-problem.png)

![Introduction: implantation failure and cause inference](assets/01-implantation-problem.png)

## (2) Build GRN, (3) Module annotation by ML, and (4) Network Transition TFs

Hierarchical clustering of gene regulatory network modules derived from scRNA-seq analysis revealed distinct receptive and non-receptive endometrial states, with mucin genes enriched in the non-receptive modules.

GRN of luminal epithelial cells was constructed and classified using a Random Forest classifier. Notably, the non-receptive group identified by the classifier was predominantly enriched in the aged uterus.

TENET was used to identify transcription factors driving the network transition toward the non-receptive state, suggesting that aging enhances their regulatory activity. Based on the same regulatory transition, key TFs were further prioritized using three graph-theoretic criteria.

**Keyword:** scRNA seq, scWGCNA, RandomForest Classifier, Network analysis, Pseudotime, TENET

![Module annotation by machine learning](assets/03-random-forest-classifier.png)

## (5) CCI - Nichnetwork analysis, (6) Validation - Spatial transcript, and (7) Validation - ATAC (peak & motif)

I identified that CD8 T cell mediated cell-cell interactions induce the upregulation of transcription factors in receptor (Luminal epithelial) cells, thereby promoting the transition to the non-receptive transcriptional network.

Spatial transcriptomic analysis showed that receptor cells interacting with CD8 T cells and exhibiting significant upregulation of key transcription factors were enriched in aged samples and displayed higher mucin expression than receptor cells without such transcription factor upregulation.

To support the regulatory mechanism identified by scRNA-seq analysis, ATAC-seq peak analysis revealed accessible regulatory regions associated with the key transcription factor at mucin genes, while motif enrichment analysis provided statistical evidence for its binding to these regulatory regions.

**Keyword:** Nichenet, Cell - Cell interaction, 10X visium, Network perturbation, ATAC, Motif, Transfactor, Epigenomics, Deconvolution

![Spatial transcriptomic validation](assets/06-spatial-validation.png)

