# Pipeline of Time-Series Pathway Identification

## (1) Introduction and (2) Pipeline

Conventional bulk RNA-seq pathway analysis identifies DEGs through pairwise comparisons (e.g., one-vs-one or one-vs-others), followed by gene set generation and pathway enrichment analysis. However, this approach cannot effectively identify pathways with progressively increasing activity across a time series.

To assess the linearity of the time-series data, I constructed a Pearson correlation matrix and evaluated its eigenvalue structure. PCA further confirmed a clear temporal trajectory, with samples sequentially aligned along the principal components. This supported the use of linear equation-based modeling to characterize the temporal dynamics of the data.

**Keyword:** Bulk RNAseq, DEG, DESeq2, Pathway enrichment, Linearity, Linear algebra, Public data

![Volcano plot and conventional pathway analysis limitation](assets/01-conventional-analysis-limitation.png)

## (3) Pipeline and (4) Validation

Gene expression was modeled using linear regression across time points, and genes with residuals below the 75 percentile were retained. By modeling positive and negative slopes separately, two temporally linear gene sets were generated, representing progressively upregulated and downregulated genes, respectively.

The resulting gene sets were validated by calculating pathway enrichment scores and comparing them with previously reported aging-associated pathways. Successful reproduction of known aging-related pathways demonstrated the validity of the proposed approach, supporting a novel methodology for constructing pathways through temporal modeling rather than conventional DEG-based analysis.

Positive slope by time : Young -> Aged pathway  
Negative slope by time : Young pathway  
Residual =< 0.602 (quantile 75%)

**Keyword:** Bulk RNAseq, Geneset modeling, Pathway enrichment, Linearity, Linear algebra, Validation

![Linear regressor modeling and validation](assets/04-validation.png)

## Discussion

### Data and analysis context

This project used bulk RNA-seq read-count data collected across multiple time points and processed with DESeq2. Pearson correlation, eigenvalue structure, and PCA were used to evaluate temporal organization before fitting gene-wise linear models. The resulting positive- and negative-slope gene sets were tested with pathway enrichment.

### Interpretation

The approach treats aging-associated expression change as a temporal trajectory rather than a collection of independent pairwise contrasts. Genes with low regression residuals were retained as reliable linear trends, and the two slope directions separated programs that increased toward aged samples from programs associated with younger samples. Recovery of known immune, inflammation, epigenetic, fertilization, and cell-cycle pathways supported the biological relevance of the gene sets.

### Considerations

Linear modeling is appropriate only when the time points adequately support a linear approximation. Limited time points, batch effects, normalization choices, residual thresholds, and sample-level confounding can influence the selected genes. Nonlinear or transient biological programs may be missed, so independent time-series datasets and nonlinear models would be useful for further validation.







