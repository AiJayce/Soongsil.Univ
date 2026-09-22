# Developing a TME Immune Cell Annotation Tool

## (1) Introduction & Pipeline

The tumor microenvironment (TME) is highly heterogeneous, making marker-based annotation of myeloid cells challenging. To overcome this limitation in scRNA seq, I trained machine learning models on public datasets and successfully annotated heterogeneous myeloid populations, enabling efficient characterization of the immune landscape.

Step1 : Collect training data  
Step2 : Build Feature  
Step3 : Training

**Keyword:** Machine Learning, Randomforest Classifier, Feature importance, Confusion matrix, ROC AUC curve, scRNA seq

![Introduction and pipeline](assets/01-training-pipeline.png)

## (2) Prediction & Validation

Highly important features identified by the Random Forest classifier were selected and used for reinforcement learning. The resulting model was subsequently applied to human single-cell data, where it demonstrated strong predictive performance.

The model predictions were further validated through enrichment analysis and the evaluation of established marker genes. These findings suggest that the model has the potential to enable accurate and effective annotation of myeloid cells within the tumor microenvironment (TME).

**Keyword:** Machine Learning, Enrichment validation, Pathway validation, Gene marker validation, UMAP characteristic, homologous genes

![Prediction and validation](assets/03-validation.png)

## Discussion

### Data and analysis context

The workflow used public scRNA-seq datasets, including macrophage reference data, to construct an expression-based training set. Feature distributions were evaluated before fitting an ensemble Random Forest classifier. The model was then applied to human TME single-cell data and assessed with UMAP, score distributions, marker genes, and pathway-enrichment results.

### Interpretation

The classifier separated M0, M1, and M2-like myeloid states and transferred the learned expression features to human single-cell data. Feature importance, confusion matrices, and ROC-AUC curves provided model-level support, while enrichment analysis and marker-gene patterns supplied biological support for the predicted labels.

### Considerations

The predictions depend on the quality, species composition, label definitions, and batch compatibility of the public training datasets. M0, M1, and M2 should be interpreted as transcriptional states rather than fully discrete biological cell types. Independent datasets, held-out samples, and experimental validation are needed to assess generalization and causal relevance.

