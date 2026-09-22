# Pipeline of Time-Series Pathway Identification

## Objective

Identify biological pathways that change continuously across aging rather than only between two selected groups.

## Analysis Workflow

### 1. Why pairwise DEG analysis is insufficient

Conventional bulk RNA-seq analysis first selects differentially expressed genes from pairwise contrasts and then performs pathway enrichment. As shown below, this approach finds group-specific genes but does not model whether a pathway changes progressively throughout the time course.

![Volcano plot and conventional pathway analysis limitation](assets/01-conventional-analysis-limitation.png)

### 2. Prepare counts and confirm a temporal trajectory

Raw count matrices were normalized with DESeq2. The PCA and sample-to-sample Pearson correlation matrix in the figure show that samples follow an ordered young-to-aged trajectory, providing the basis for a linear time-series model.

![DESeq2 preprocessing, PCA, and sample correlation across time points](assets/02-preprocess-and-linearity.png)

### 3. Select genes with linear temporal behavior

Each gene was fitted with simple linear regression across time points. Genes with residuals at or below the 75th-percentile threshold were retained as well-described linear trajectories. The slope separated increasing programs from decreasing programs.

![Positive and negative temporal trajectories, linear model, and residual threshold](assets/03-linear-regression.png)

### 4. Test the resulting pathway programs

The positive-slope gene set was enriched toward aged samples, while the negative-slope gene set was enriched toward young samples. Pathway analysis recovered expected immune and inflammation, epigenetic modification, fertilization, and cell-cycle programs, validating the time-series approach.

![Enrichment-score and pathway-level validation of slope-based gene sets](assets/04-validation.png)

## Methods

- Bulk RNA-seq preprocessing with DESeq2
- Pearson correlation and principal-component analysis
- Simple linear regression and residual filtering
- Positive- and negative-slope gene-set construction
- Gene-set enrichment and pathway validation

## Key Finding

Temporal linear modeling recovered biologically relevant age-associated pathways that are difficult to capture with conventional pairwise DEG-based pathway analysis.







