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

