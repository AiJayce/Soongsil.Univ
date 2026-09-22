# Developing a TME Immune Cell Annotation Tool

## Objective

Build a reproducible myeloid-cell annotation workflow for heterogeneous tumor microenvironment scRNA-seq data, where marker genes alone are insufficient to distinguish immune states.

## Analysis Workflow

### 1. Training-data collection and feature construction

Public macrophage scRNA-seq data were used as training data. The upper panels in the figure show the state-specific feature distributions, while the lower panels show the fitted ensemble Random Forest, confusion matrix, feature rankings, and ROC-AUC evaluation.

![Training-data collection, feature construction, and ensemble model training](assets/01-training-pipeline.png)

### 2. Select a classifier-ready feature set

The feature-importance and ROC-AUC panels identify the expression scores that best separate the reference myeloid states. These scores form the classifier-ready feature set transferred to human TME data.

### 3. Prediction in human single-cell data

The trained model labels human cells as M0, M1, M2, or unknown. The UMAP shows the predicted spatial separation, while the violin and dot plots verify that each label has the expected score and marker-gene pattern.

![Predicted M0, M1, and M2 myeloid states](assets/02-prediction.png)

### 4. Biological validation

Gene-set enrichment validates the predicted programs. M1 cells show inflammatory and interferon-gamma-associated activity, while M2 cells show distinct hypoxia and epithelial-mesenchymal-transition programs, providing biological support for the annotations.

![Pathway-enrichment validation of predicted myeloid states](assets/03-validation.png)

## Methods

- Public scRNA-seq training data integration
- Feature construction and importance ranking
- Ensemble Random Forest classification
- Confusion-matrix and ROC-AUC evaluation
- UMAP, marker-gene, pathway, and enrichment validation

## Key Finding

The classifier transferred informative myeloid-state features from public reference data to human TME data and produced biologically supported myeloid annotations.

