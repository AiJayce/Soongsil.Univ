# Implantation Failure in Uterine Aging

## (1) Introduction

<img width="3252" height="1184" alt="Implantation Failure in Uterine Aging-1" src="https://github.com/user-attachments/assets/497c4eed-3719-48cb-9ba0-863e06132f99" />

I elucidated the transcriptomic mechanisms underlying implantation failure in the aged uterus by demonstrating how aging drives a transition toward a non-receptive endometrial network. Because successful embryo implantation, including after IVF, requires low mucin expression, I discovered that senescent immune cell–mediated interactions promote a transcriptomic transition toward a mucin-overexpression network, thereby remodeling the uterine microenvironment into an implantation-unfavorable state.


## (2) Build GRN

<img width="1114" height="1123" alt="Implantation Failure in Uterine Aging-2" src="https://github.com/user-attachments/assets/fee190d8-62b4-4661-ab1c-3419d96e0aca" />

Hierarchical clustering of gene regulatory network modules derived from scRNA-seq analysis revealed distinct receptive and non-receptive endometrial states, with mucin genes enriched in the non-receptive modules.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/scWGCNA.R

## (3) Module annotation by ML

<img width="1114" height="1123" alt="Implantation Failure in Uterine Aging-3" src="https://github.com/user-attachments/assets/688c0235-2b3d-4f39-a1c8-2e48fddbe54c" />

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb

## (4) Network Transition TFs

<img width="1117" height="1123" alt="Implantation Failure in Uterine Aging-4" src="https://github.com/user-attachments/assets/c1e3ea80-1d43-4bcc-94a5-855252d5a57f" />

TENET was used to identify transcription factors driving the network transition toward the non-receptive state, suggesting that aging enhances their regulatory activity. Based on the same regulatory transition, key TFs were further prioritized using three graph-theoretic criteria.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Cell%20oracle.ipynb

## (5) CCI - Nichnetwork analysis

<img width="1114" height="1113" alt="Implantation Failure in Uterine Aging-5" src="https://github.com/user-attachments/assets/06766976-a30d-462a-aefb-5d9dfbc34cd2" />

I identified that CD8 T cell mediated cell–cell interactions induce the upregulation of transcription factors in receptor (Luminal epithelial) cells, thereby promoting the transition to the non-receptive transcriptional network.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/nichenet-1.R, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/nichenet-2.R, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/nichenet-3_LRnetwork.R

(6) Validation - Spatial transcript

<img width="1054" height="1117" alt="Implantation Failure in Uterine Aging-6" src="https://github.com/user-attachments/assets/6039bb7a-49a3-49e1-a50a-4e3b0421dcf0" />

Spatial transcriptomic analysis showed that receptor cells interacting with CD8 T cells and exhibiting significant upregulation of key transcription factors were enriched in aged samples and displayed higher mucin expression than receptor cells without such transcription factor upregulation.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/spatial/Visium(Old%26Young).ipynb

(7) Validation - ATAC (peak & motif)

<img width="1117" height="1118" alt="Implantation Failure in Uterine Aging-7" src="https://github.com/user-attachments/assets/434e0ab6-2fb5-4f2e-a6e6-3d840541c135" />

To support the regulatory mechanism identified by scRNA-seq analysis, ATAC-seq peak analysis revealed accessible regulatory regions associated with the key transcription factor at mucin genes, while motif enrichment analysis provided statistical evidence for its binding to these regulatory regions.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Binding%20region.R

### Data and analysis context

This project integrated luminal epithelial scRNA-seq, spatial transcriptomics, ATAC-seq, and public reference data. scWGCNA and trajectory analysis defined receptive and non-receptive regulatory states, the ensemble Random Forest classifier annotated network modules, and TENET prioritized transcription factors involved in the transition. NicheNet was used to investigate CD8 T-cell to epithelial communication.

### Interpretation

The results support a model in which uterine aging is associated with a shift toward a mucin-overexpression, non-receptive epithelial network. CD8 T-cell-mediated signals were associated with transcription-factor upregulation in luminal epithelial cells, while spatial data linked high interaction states with increased mucin expression. ATAC-seq peaks and HOMER motif enrichment provided epigenomic support for regulatory activity at mucin-associated loci.

### Considerations

Network modules, pseudotime, NicheNet communication scores, and spatial correlations are inferential rather than direct measurements of causality. The proposed transition should therefore be tested with perturbation experiments, independent uterine-aging cohorts, and functional implantation assays.

