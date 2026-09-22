# Macrophage Transition in Gastric Cancer TME


## (1) Introduction

<img width="1221" height="1186" alt="1" src="https://github.com/user-attachments/assets/fa22cdf7-5561-4ffe-89c1-946a704c804c" />

Macrophages in the TME are broadly classified into M1 and M2 states with distinct roles, and are strongly influenced by CAFs; therefore, this project aims to investigate how CAFs drive macrophages toward M1 or M2 states.

<img width="2071" height="1188" alt="2" src="https://github.com/user-attachments/assets/14cf8333-1bda-4faa-b25a-8acdd19d0054" />

Macrophages are highly heterogeneous in the TME, making accurate annotation challenging. Marker-based annotation provides higher resolution but can reduce accuracy, whereas automated annotation improves accuracy at the cost of resolution. To overcome this trade-off, the macrophage classification model developed in a previous project was applied, combining the advantages of both approaches.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/TME-ML-Classifier/GC%20mac%20final.ipynb

## (2) Annotation Validation

<img width="1221" height="1186" alt="3" src="https://github.com/user-attachments/assets/0265886b-3df2-4d46-9c3c-f5fdc81bfe92" />

Classification was validated using marker-based analysis. Since the training data were FACS-sorted, they provided reliable M1 and M2 labels. The classified populations showed consistent trends in both marker expression and module scores, supporting the validity of the classification.

<img width="2064" height="1190" alt="4" src="https://github.com/user-attachments/assets/67e500e2-5b3b-4a18-a590-e7a67ba5c28e" />

Pathway-based validation further confirmed that the M1 and M2 states were well captured. M1 macrophages showed activation of cytokine- and inflammation-related pathways, whereas M2 macrophages exhibited immune-suppressive and EMT-related pathways associated with a cold tumor TME.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/TME-ML-Classifier/Machinelearning.ipynb

## (3) Cell Cell interaction analysis

<img width="1787" height="1186" alt="5" src="https://github.com/user-attachments/assets/dece7ade-2c71-463b-ba5d-4eba3ed04ef2" />

Immune cells in the TME are strongly influenced by surrounding cells. To identify which cells most strongly influence macrophages, cell–cell communication analysis revealed that ECM-rich fibroblasts (ECM-rich FBs), a CAF subtype, showed the strongest interactions with macrophages as ligand-producing cells. Among these interactions, the CSF pathway ranked first on the receptor side.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Gastric-Cancer/GC%20datamerge.ipynb

## (4) CSF pathway

<img width="1512" height="1187" alt="6" src="https://github.com/user-attachments/assets/09c86b48-2941-48cc-aea3-380322e35228" />

Previous studies have reported that CSF pathway activation induces STAT3 activation, promoting M1-to-M2 repolarization. 
In this project, I observed macrophage dynamics driven by CAF interactions in the TME, suggesting that CAF-mediated M1-to-M2 transition may contribute to the development of a cold tumor environment.

## Discussion

### Data and analysis context

Public scRNA-seq datasets, including FACS-sorted macrophage reference data, were used to construct an expression-based training set for M0, M1, and M2-like macrophage classification. An ensemble Random Forest classifier was trained on the gene expression features and applied to human gastric cancer TME scRNA-seq data. The predicted states were evaluated using UMAP, score distributions, marker gene expression, module scores, and pathway-enrichment analysis. CellChat, pySCENIC, and Monocle were subsequently used to investigate CAF–macrophage interactions, transcriptional regulation, and macrophage state transitions.

### Interpretation

The classifier distinguished M0, M1, and M2-like macrophage states and transferred the learned transcriptional features to human gastric cancer TME data. Marker genes, module scores, and pathway-enrichment patterns supported the biological characteristics of the predicted states. CellChat identified ECM-rich fibroblasts as a major interacting CAF subtype and revealed active CSF1–CSF1R signaling toward macrophages. Integration with pySCENIC and pseudotime analysis suggested a potential CSF1–CSF1R–STAT3-mediated regulatory mechanism associated with M1-to-M2-like repolarization.

### Considerations

The predicted macrophage states depend on the quality, species composition, label definitions, and batch compatibility of the public training datasets. M0, M1, and M2-like states should be interpreted as transcriptional states rather than fully discrete biological cell types. Although the integrated CellChat, pySCENIC, and pseudotime analyses support a potential CSF1–CSF1R–STAT3-mediated transition, these findings represent computational associations and do not establish causality. Independent datasets and experimental validation are required to confirm the proposed mechanism and its relevance to cold tumor formation.

The predictions depend on the quality, species composition, label definitions, and batch compatibility of the public training datasets. M0, M1, and M2 should be interpreted as transcriptional states rather than fully discrete biological cell types. Independent datasets, held-out samples, and experimental validation are needed to assess generalization and causal relevance.

