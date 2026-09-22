# Developing a TME Immune Cell Annotation Tool

## (1) Introduction & Pipeline

<img width="960" height="540" alt="tumor-heterogeneity-353504-960x540" src="https://github.com/user-attachments/assets/060c9110-a093-4b7e-84ee-7b546033f3cc" />

The tumor microenvironment (TME) is highly heterogeneous, making marker-based annotation of myeloid cells challenging. To overcome this limitation in scRNA seq, I trained machine learning models on public datasets and successfully annotated heterogeneous myeloid populations, enabling efficient characterization of the immune landscape.

## (2) Pipeline

<img width="3118" height="1222" alt="Developing a TME immune cell annotation tool-1" src="https://github.com/user-attachments/assets/a0f1e347-aebf-4a9d-a9f1-40145e2b3b49" />

Step1 : Collect training data : Query Sorting scRNA dataset from GEO (GSM3272967)

Step2 : Build Feature : Before training, perform EDA to identify and construct the optimal training features.

Step3 : Training and validation (Confusion matrix, ROC AUC curve)

### code : https://github.com/AiJayce/Soongsil.Univ/blob/TME-ML-Classifier/Machinelearning.ipynb

## (3) Prediction & Validation

<img width="1335" height="1212" alt="Developing a TME immune cell annotation tool-2" src="https://github.com/user-attachments/assets/c243768c-463d-49e1-a4ec-de09e035d145" />

Highly important features identified by the Random Forest classifier were selected and used for reinforcement learning. The resulting model was subsequently applied to human single-cell data, where it demonstrated strong predictive performance. 

### code : https://github.com/AiJayce/Soongsil.Univ/blob/TME-ML-Classifier/Machinelearning.ipynb

## (4) Validation

<img width="1782" height="1205" alt="Developing a TME immune cell annotation tool-3" src="https://github.com/user-attachments/assets/32ae5f17-c4ea-4904-afeb-be87b9a558ea" />

The model predictions were further validated through enrichment analysis and the evaluation of established marker genes. These findings suggest that the model has the potential to enable accurate and effective annotation of myeloid cells within the tumor microenvironment (TME).

### code : https://github.com/AiJayce/Soongsil.Univ/blob/TME-ML-Classifier/Machinelearning.ipynb

## Discussion

### Data and analysis context

The workflow used public scRNA-seq datasets, including macrophage reference data, to construct an expression-based training set. Feature distributions were evaluated before fitting an ensemble Random Forest classifier. The model was then applied to human TME single-cell data and assessed with UMAP, score distributions, marker genes, and pathway-enrichment results.

### Interpretation

The classifier separated M0, M1, and M2-like myeloid states and transferred the learned expression features to human single-cell data. Feature importance, confusion matrices, and ROC-AUC curves provided model-level support, while enrichment analysis and marker-gene patterns supplied biological support for the predicted labels.

### Considerations

The predictions depend on the quality, species composition, label definitions, and batch compatibility of the public training datasets. M0, M1, and M2 should be interpreted as transcriptional states rather than fully discrete biological cell types. Independent datasets, held-out samples, and experimental validation are needed to assess generalization and causal relevance.

