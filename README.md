# Implantation Failure in Uterine Aging

## (1) Introduction

<img width="3252" height="1184" alt="Implantation Failure in Uterine Aging-1" src="https://github.com/user-attachments/assets/497c4eed-3719-48cb-9ba0-863e06132f99" />

I elucidated the transcriptomic mechanisms underlying implantation failure in the aged uterus by demonstrating how aging drives a transition toward a non-receptive endometrial network. Because successful embryo implantation, including after IVF, requires low mucin expression, I discovered that senescent immune cell–mediated interactions promote a transcriptomic transition toward a mucin-overexpression network, thereby remodeling the uterine microenvironment into an implantation-unfavorable state.

## (2) Data information

<img width="4280" height="1731" alt="그림1" src="https://github.com/user-attachments/assets/2ce24416-f462-4633-bf93-76fffe010235" />

The data analysis was conducted at two levels: mouse and human. First, the mouse dataset was divided into Young and Aged groups based on a 20-week age threshold. Comparative analysis between the two groups was performed to identify biomarkers associated with implantation impairment during aging. For the human analysis, a dataset comprising individuals aged 18–46 years was constructed, covering a range in which the characteristics of Young and Aged groups were less clearly distinguishable than in the mouse model. This approach aimed to minimize the influence of age-related differences and focus on the effects of selected genes on the Gene Regulatory Network (GRN). Rather than simply identifying age-associated differences, the analysis focused on understanding how age-independent changes at the gene level are reflected in the structure and regulatory properties of GRNs associated with implantation impairment.

Another major distinction between the mouse and human datasets is the periodicity of hormonal cycles.

## (3) Build GRN

<img width="1114" height="1123" alt="Implantation Failure in Uterine Aging-2" src="https://github.com/user-attachments/assets/fee190d8-62b4-4661-ab1c-3419d96e0aca" />

Hierarchical clustering of gene regulatory network modules derived from scRNA-seq analysis revealed distinct receptive and non-receptive endometrial states, with mucin genes enriched in the non-receptive modules.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/scWGCNA.R

## (3) Module annotation by ML

<img width="5319" height="1623" alt="그림2" src="https://github.com/user-attachments/assets/f994f4ec-2009-4cd0-840d-ce9cc65245b1" />

Among the 11 generated module networks, machine learning techniques were used to identify modules that significantly captured cellular network characteristics, following the steps below.

1. WCSS and Silhouette analyses were performed to evaluate the clustering characteristics of the 11 networks. The results identified two modules with distinct clustering properties.

2. PCA was applied to cell-level gene expression data to reduce dimensionality. Based on the previously determined optimal number of clusters (K = 2), the K-means algorithm was used to classify cells into two clusters through unsupervised learning.

3. To identify the modules that best represented the characteristics of the two clusters, a Random Forest model was trained and Feature Importance was analyzed. The results showed that Module 7 best characterized Cluster A, whereas Module 2 best characterized Cluster B.

4. To validate the biological characteristics of the two clusters identified through unsupervised learning, the expression levels of Hbegf, a gene associated with promoting implantation, were compared. Hbegf expression was statistically significantly higher in the cluster corresponding to Module 7 than in the cluster corresponding to Module 2. Furthermore, DEG-based pathway analysis revealed pathways associated with implantation impairment in Module 2 and pathways associated with implantation promotion in Module 7.

5. Based on the results of unsupervised learning and biological validation, Module 2 was designated as the Non-receptive state, while Module 7 was designated as the Receptive state.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb

## (4) Network Transition TFs

<img width="1117" height="1123" alt="Implantation Failure in Uterine Aging-4" src="https://github.com/user-attachments/assets/c1e3ea80-1d43-4bcc-94a5-855252d5a57f" />

TENET was used to identify transcription factors driving the network transition toward the non-receptive state, suggesting that aging enhances their regulatory activity. Based on the same regulatory transition, key TFs were further prioritized using three graph-theoretic criteria.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Cell%20oracle.ipynb

## (5) Pseudotime Validation & Gene perterbation Simulation

<img width="4584" height="1788" alt="그림3" src="https://github.com/user-attachments/assets/c8075665-4708-4f40-b481-b6c0613fbff4" />

1. To validate the biological plausibility of the pseudotime used in TENET, the similarity of cellular transition vector fields inferred by three different algorithms—Monocle, CellRank, and RNA velocity—was compared. Monocle was used as the baseline algorithm, inferring trajectories by specifying both the starting and terminal points. In contrast, CellRank estimates cellular state transition directions based on specified starting points, whereas RNA velocity infers transition directions by modeling temporal changes in gene expression using differential equations without explicitly specifying either starting or terminal points. Four different statistical tests based on the dot products between vectors derived from the three algorithms revealed statistically significant similarities in transition directions. These results indicate that the pseudotime trajectory established using Monocle is consistent with the trajectories inferred by other algorithms and exhibits a certain degree of biological plausibility.

2. After validating the pseudotime trajectory, CellOracle was used to evaluate the functional role of Mxd1, a transcription factor identified by TENET, through a gene knockout (KO) simulation. Specifically, the effect of Mxd1 KO on the transition of cells toward the Non-receptive state was investigated. The results showed that knockout of Mxd1 reversed the pseudotime transition direction that had previously progressed toward the Non-receptive state. This finding suggests that Mxd1 may be involved in cellular state transitions toward the Non-receptive state.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Cell%20oracle.ipynb


## (6) CCI - Nichnetwork analysis

<img width="3486" height="1144" alt="그림4" src="https://github.com/user-attachments/assets/7c0fce06-984b-48e4-aa01-42c755a2d407" />

I identified that CD8 T cell-mediated cell–cell interactions become more active during aging, inducing the upregulation of transcription factors in receptor (luminal epithelial) cells and thereby promoting the transition to the non-receptive transcriptional network.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/nichenet-1.R, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/nichenet-2.R, https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/nichenet-3_LRnetwork.R


## (7) Expression S100a8 in CD8 T cells

<img width="3582" height="1863" alt="그림9" src="https://github.com/user-attachments/assets/d5e01b57-d921-4cfa-beb4-9df7a17acc6f" />

According to the NicheNet analysis, S100a8 expression in CD8 T cells was identified as a potential key mechanism underlying implantation impairment. However, since S100a8 is generally expressed predominantly in the myeloid lineage, further investigation was required to determine whether it is also expressed in CD8 T cells. A literature review revealed that S100a8 can also be expressed in aged CD8 T cells. To validate this finding, S100a8 expression in CD8 T cells was compared between the Young and Aged groups. The results reproduced the expression pattern reported in the literature, suggesting that CD8 T cells may exhibit age-associated S100a8 expression in the aged uterus and potentially contribute to implantation impairment through this mechanism.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Uterus_Epithelial.ipynb

## (8) Validation - Spatial transcript

<img width="1054" height="1117" alt="Implantation Failure in Uterine Aging-6" src="https://github.com/user-attachments/assets/6039bb7a-49a3-49e1-a50a-4e3b0421dcf0" />

Spatial transcriptomic analysis showed that receptor cells interacting with CD8 T cells and exhibiting significant upregulation of key transcription factors were enriched in aged samples and displayed higher mucin expression than receptor cells without such transcription factor upregulation. Samples in which the S100a8–Tlr4 interaction was associated with a significant increase in Mxd1 expression were classified as High correlation, whereas those without a statistically significant increase were classified as Low correlation.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/spatial/Visium(Old%26Young).ipynb

## (9) Validation - Spatial Plot

<img width="3681" height="1863" alt="그림5" src="https://github.com/user-attachments/assets/1f4d4d8a-b0a4-45c9-be50-22cc9ba9bc43" />

To investigate whether the identified relationship between S100a8, Tlr4, and Mxd1 was preserved spatially, spatial transcriptomic data were visualized. The results showed that regions with stronger S100a8–Tlr4 interactions exhibited higher Mxd1 expression, and these Mxd1-high regions also showed elevated Muc4 expression, supporting the spatial preservation of this regulatory relationship.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/spatial/Visium(Old%26Young).ipynb

## (10) Validation - ATAC (peak & motif)

<img width="1117" height="1118" alt="Implantation Failure in Uterine Aging-7" src="https://github.com/user-attachments/assets/434e0ab6-2fb5-4f2e-a6e6-3d840541c135" />

To support the regulatory mechanism identified by scRNA-seq analysis, ATAC-seq peak analysis revealed accessible regulatory regions associated with the key transcription factor at mucin genes, while motif enrichment analysis provided statistical evidence for its binding to these regulatory regions.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/Binding%20region.R

## (11) Validation - Human level

Thus far, three markers associated with implantation impairment have been investigated at the mouse level. Their expression patterns and regulatory networks were then examined at the human level to determine whether these findings could be reproduced across species.

<img width="3837" height="1653" alt="그림6" src="https://github.com/user-attachments/assets/60158667-552f-40f6-a4d3-81f8be41fc7e" />

1. The hormone cycle-dependent changes in epithelial subtypes observed in mice were also reproduced in humans.
2. Consistent with the TENET results in mice, Mxd1 and Muc4 showed a positive correlation in humans, suggesting a potential regulatory relationship.


<img width="4378" height="1654" alt="그림7" src="https://github.com/user-attachments/assets/2db3a567-7be7-4b14-aadc-900f641dc2c8" />

Previous studies reporting S100a8 expression in aged CD8 T cells were conducted at the mouse level. Therefore, it was necessary to investigate whether this expression pattern could also be reproduced in humans. To address this, the following analyses were performed.

1. Lymphoid cells were divided into groups based on S100a8 expression, and CDKN2A expression, a marker of cellular senescence, was compared at the sample level using pseudobulk analysis. Most samples with detectable S100a8 expression in the lymphoid lineage exhibited higher CDKN2A expression. These results suggest that the association between S100a8 expression in lymphoid cells and aging observed in mice is also reproducible in humans.

2. As this project specifically focused on S100a8 expression in CD8 T cells, samples were classified into S100a8-high and S100a8-low groups based on S100a8 expression levels in CD8 T cells.

3. Comparison between the two groups revealed higher Muc4 and Mxd1 expression in Luminal Epithelium in the S100a8-high group, suggesting an association between S100a8 expression in CD8 T cells and gene expression changes in Luminal Epithelium.

4. Comparison of lymphoid cell population proportions showed that the proportion of CD8 T cells was significantly higher in the S100a8-high group. This finding indicates that elevated S100a8 expression in CD8 T cells is also associated with changes in the composition of the lymphoid cell population.

<img width="3284" height="1839" alt="그림8" src="https://github.com/user-attachments/assets/ec25f9d7-d97f-4a59-ba1e-02169d15f682" />

Pathway analysis revealed that stronger Tlr4–S100a8 interactions between Luminal Epithelium and CD8 T cells were associated with increased Mxd1 expression in receptor cells, potentially distinguishing receptive from non-receptive states and thereby influencing implantation outcomes.
These findings demonstrate that the results obtained from mouse-level analyses are reproducible in human-level analyses, supporting the biological relevance of the proposed regulatory interactions and the validity of the analytical approach.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Human/Human_sc.ipynb

## Discussion

### Data and analysis context

This project integrated luminal epithelial scRNA-seq, spatial transcriptomics, ATAC-seq, and public reference data. scWGCNA and trajectory analysis defined receptive and non-receptive regulatory states, the ensemble Random Forest classifier annotated network modules, and TENET prioritized transcription factors involved in the transition. NicheNet was used to investigate CD8 T-cell to epithelial communication.

### Interpretation

The results support a model in which uterine aging is associated with a shift toward a mucin-overexpression, non-receptive epithelial network. CD8 T-cell-mediated signals were associated with transcription-factor upregulation in luminal epithelial cells, while spatial data linked high interaction states with increased mucin expression. ATAC-seq peaks and HOMER motif enrichment provided epigenomic support for regulatory activity at mucin-associated loci.

### Considerations

Network modules, pseudotime, NicheNet communication scores, and spatial correlations are inferential rather than direct measurements of causality. The proposed transition should therefore be tested with perturbation experiments, independent uterine-aging cohorts, and functional implantation assays.

