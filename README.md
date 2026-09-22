# Macrophage Transition in Gastric Cancer TME


## (1) Introduction

<img width="1221" height="1186" alt="1" src="https://github.com/user-attachments/assets/fa22cdf7-5561-4ffe-89c1-946a704c804c" />

Macrophages in the TME are broadly classified into M1 and M2 states with distinct roles, and are strongly influenced by CAFs; therefore, this project aims to investigate how CAFs drive macrophages toward M1 or M2 states.

<img width="2071" height="1188" alt="2" src="https://github.com/user-attachments/assets/14cf8333-1bda-4faa-b25a-8acdd19d0054" />

Macrophages are highly heterogeneous in the TME, making accurate annotation challenging. Marker-based annotation provides higher resolution but can reduce accuracy, whereas automated annotation improves accuracy at the cost of resolution. To overcome this trade-off, the macrophage classification model developed in a previous project was applied, combining the advantages of both approaches.


## (2) Annotation Validation

<img width="1221" height="1186" alt="3" src="https://github.com/user-attachments/assets/0265886b-3df2-4d46-9c3c-f5fdc81bfe92" />

Classification was validated using marker-based analysis. Since the training data were FACS-sorted, they provided reliable M1 and M2 labels. The classified populations showed consistent trends in both marker expression and module scores, supporting the validity of the classification.

<img width="2064" height="1190" alt="4" src="https://github.com/user-attachments/assets/67e500e2-5b3b-4a18-a590-e7a67ba5c28e" />

Pathway-based validation further confirmed that the M1 and M2 states were well captured. M1 macrophages showed activation of cytokine- and inflammation-related pathways, whereas M2 macrophages exhibited immune-suppressive and EMT-related pathways associated with a cold tumor TME.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/TME-ML-Classifier/Machinelearning.ipynb

## (3) Cell Cell interaction analysis

<img width="1787" height="1186" alt="5" src="https://github.com/user-attachments/assets/dece7ade-2c71-463b-ba5d-4eba3ed04ef2" />

Immune cells in the TME are strongly influenced by surrounding cells. To identify which cells most strongly influence macrophages, cell–cell communication analysis revealed that ECM-rich fibroblasts (ECM-rich FBs), a CAF subtype, showed the strongest interactions with macrophages as ligand-producing cells. Among these interactions, the CSF pathway ranked first on the receptor side.

## (4) CSF pathway

<img width="1512" height="1187" alt="6" src="https://github.com/user-attachments/assets/09c86b48-2941-48cc-aea3-380322e35228" />

Previous studies have reported that CSF pathway activation induces STAT3 activation, promoting M1-to-M2 repolarization. 
In this project, I observed macrophage dynamics driven by CAF interactions in the TME, suggesting that CAF-mediated M1-to-M2 transition may contribute to the development of a cold tumor environment.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/TME-ML-Classifier/Machinelearning.ipynb

## (4) Validation


## Discussion

### Data and analysis context

The workflow used public scRNA-seq datasets, including macrophage reference data, to construct an expression-based training set. Feature distributions were evaluated before fitting an ensemble Random Forest classifier. The model was then applied to human TME single-cell data and assessed with UMAP, score distributions, marker genes, and pathway-enrichment results.

### Interpretation

The classifier separated M0, M1, and M2-like myeloid states and transferred the learned expression features to human single-cell data. Feature importance, confusion matrices, and ROC-AUC curves provided model-level support, while enrichment analysis and marker-gene patterns supplied biological support for the predicted labels.

### Considerations

The predictions depend on the quality, species composition, label definitions, and batch compatibility of the public training datasets. M0, M1, and M2 should be interpreted as transcriptional states rather than fully discrete biological cell types. Independent datasets, held-out samples, and experimental validation are needed to assess generalization and causal relevance.

