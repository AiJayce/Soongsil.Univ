# Pipeline of Time-Series Pathway Identification

## (1) Introduction 

<img width="1554" height="1155" alt="Pipeline of Time-Series Pathway Identification-1" src="https://github.com/user-attachments/assets/38e07236-ebf7-42b7-9fef-9b301d1f940c" />

Conventional bulk RNA-seq pathway analysis identifies DEGs through pairwise comparisons (e.g., one-vs-one or one-vs-others), followed by gene set generation and pathway enrichment analysis. However, this approach cannot effectively identify pathways with progressively increasing activity across a time series.

## (2) Data information

<img width="3722" height="1956" alt="그림1" src="https://github.com/user-attachments/assets/f22b542c-1416-4d67-a3ea-7597546f80c8" />

The data were obtained from publicly available datasets, and preprocessing was performed using FASTQ files to enable RNA-level analysis. The dataset consisted of three age groups classified by developmental stage, with three samples in each group. Using 25 weeks as the cutoff, the samples were broadly classified into Young and Aged groups. The goal of this project was to identify pathways associated with aging by analyzing differences in RNA expression between the two groups.

## (3) Data preprocess

### NGS data preprocess
<img width="4339" height="2216" alt="그림2" src="https://github.com/user-attachments/assets/3f1077a0-37b9-4fde-9c48-fd666a6fb69a" />

During NGS preprocessing and quality control (QC), some abnormalities were observed in the FastQC reports. However, these findings were considered potentially attributable to technical biases and the characteristics of the RNA-seq data rather than necessarily indicating poor data quality.

 1) Random Hexamer Priming Bias
Random hexamer priming bias occurs because random hexamers do not bind uniformly to all possible 6-mer sequences during reverse transcription. Differences in GC content and binding stability can cause certain hexamers to be preferentially selected, resulting in sequence composition bias near the priming sites. This bias is often most pronounced at the beginning of the forward reads (near the 5′ end), where elevated GC content was observed in the FastQC per-base sequence content profile. As reverse transcription proceeds, the influence of the initial priming bias generally diminishes, leading to a more normalized sequence composition toward the end of the reads. Reverse reads typically exhibit a weaker pattern of this specific bias, although minor fluctuations may still occur due to other technical factors.

 2) Sequence Duplication Levels
Elevated sequence duplication levels were observed in the FastQC reports. This may be partly explained by the high expression of housekeeping genes such as GAPDH, which can contribute a substantial proportion of the total sequencing reads. Highly expressed transcripts are more likely to generate repeatedly sequenced fragments, increasing the apparent duplication rate. However, high duplication levels alone do not establish that the observed duplication is exclusively attributable to highly expressed genes; PCR amplification and low library complexity can also contribute. Therefore, this QC metric should be interpreted in the context of gene expression profiles, library preparation methods, and other quality indicators.

### Gene Alignment (BAM)

<img width="4315" height="2278" alt="그림3" src="https://github.com/user-attachments/assets/e70c9f53-855b-462e-9ab7-466b2c08fa06" />

STAR (Spliced Transcripts Alignment to a Reference) was used to align RNA-seq reads to the reference genome. STAR is a fast, splice-aware aligner that accurately identifies canonical and non-canonical splice junctions, including novel splice sites, and can detect chimeric transcripts. Compared with TopHat2, STAR offers substantially faster alignment, although it requires considerable memory resources. It also provides high sensitivity for mapping reads across complex exon–intron structures.

After alignment, the proportions of forward and reverse reads were examined. Their relatively balanced distribution suggested that the alignment process did not exhibit an obvious strand-specific imbalance. However, this observation alone is insufficient to conclude that alignment was free from chromosome-specific bias. Further evaluation of mapping rates, chromosome-level read distributions, and strand specificity would be required to assess alignment quality comprehensively.

### Read count preprocess

<img width="1554" height="517" alt="12" src="https://github.com/user-attachments/assets/8100ccdc-fae4-4bb3-8499-c177257594e5" />
<img width="4403" height="1681" alt="그림4" src="https://github.com/user-attachments/assets/29b84269-805c-48e6-9bfb-b3e419e19490" />

Raw read counts obtained after alignment cannot be directly compared across samples because of differences in sequencing depth, library size, and gene expression variability. RNA-seq count data typically exhibit a mean–variance relationship in which highly expressed genes have greater variance than lowly expressed genes. Therefore, normalization and variance stabilization are essential for reliable sample comparisons and downstream analyses.

DESeq2 was used to normalize count data by estimating sample-specific size factors that account for differences in sequencing depth and library composition. Subsequently, the Variance Stabilizing Transformation (VST) was applied to reduce the dependence of variance on mean expression. Unlike a conventional logarithmic transformation, VST uses the estimated mean–dispersion relationship to stabilize variance across the expression range, making the transformed data more suitable for PCA, clustering, and pathway-level exploration.

QC results showed reduced variability in expression distributions after VST, indicating improved consistency across samples. Changes in the zero-count ratio were also examined; however, VST itself does not eliminate zero counts, so any observed reduction requires interpretation in the context of the preprocessing and calculation methods used.

Overall, DESeq2 normalization and VST provided a more consistent representation of gene expression for investigating age-associated transcriptional changes and identifying pathways potentially involved in aging. Although this process helps address technical variability, VST alone does not directly remove batch effects.

## (4) Data analysis

### Assess Time-Series Linearity

<img width="4191" height="2359" alt="그림5" src="https://github.com/user-attachments/assets/5d8b2748-e4dc-4b99-bce4-528a9c0d49dc" />

To assess the linearity of the time-series data, I compared pairwise Pearson correlation and cosine similarity between samples. The two metrics showed similar patterns, suggesting that the relationships between samples were largely consistent in terms of both expression magnitude and directional similarity. This supported the assumption that the temporal changes could be approximated using a linear framework. PCA further revealed a sequential arrangement of samples along the principal components, providing additional evidence of an ordered temporal trajectory. Based on these observations, linear equation–based modeling was considered a reasonable approach for inferring temporal relationships and characterizing the progression of gene expression over time.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/RNA-seq/class.ipynb

### Linear Regressor modeling

<img width="3659" height="1636" alt="그림6" src="https://github.com/user-attachments/assets/15e6c2a8-5455-47ff-a62a-b11719d8f522" />

Gene expression was modeled using linear regression across time points, and genes with residuals below the 75 percentile were retained. By modeling positive and negative slopes separately, two temporally linear gene sets were generated, representing progressively upregulated and downregulated genes, respectively.

<img width="4328" height="2401" alt="그림7" src="https://github.com/user-attachments/assets/b1f9c2cf-1d7d-4b1b-be16-833b56b1b97d" />

By comparing gene-level expression trends using the median expression values of three samples at each time point, I observed approximately linear expression patterns across time points. These results supported the feasibility of constructing gene set models to characterize temporal expression dynamics and suggested that the models could effectively capture time-dependent transcriptional trends.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/RNA-seq/class.ipynb

## (5) Validation

<img width="1691" height="1166" alt="Pipeline of Time-Series Pathway Identification-4" src="https://github.com/user-attachments/assets/33463ade-a324-49ff-8c0a-994721cbe4e6" />

The resulting gene sets were validated by calculating pathway enrichment scores and comparing them with previously reported aging-associated pathways. Successful reproduction of known aging-related pathways demonstrated the validity of the proposed approach, supporting a novel methodology for constructing pathways through temporal modeling rather than conventional DEG-based analysis.

## (6) Supplementary Analysis

### Nonlinear Time-Series Pathway Analysis

<img width="4318" height="2368" alt="그림8" src="https://github.com/user-attachments/assets/ad43f437-97f7-4782-a1c0-f458d132cc9e" />

Pathway activity does not always follow a linear trajectory over time. Some pathways may exhibit transient activation, increasing during a specific period before declining at later time points. The previously developed linear model cannot adequately capture these non-monotonic temporal patterns, potentially overlooking biologically meaningful transient responses.

To capture these patterns, a piecewise linear mixture model was constructed to represent a rise-and-fall expression trajectory.

Pivot Timepoint Selection: A candidate timepoint was selected as a pivot, representing the transition point between the increasing and decreasing phases of gene expression.

Piecewise Linear Regression: Two separate linear regression models were fitted around the pivot. The pre-pivot segment represented an increasing trend with a positive slope, while the post-pivot segment represented a decreasing trend with a negative slope. This enabled the model to approximate transient expression peaks using two linear components.

Significance-Based Filtering: A 25% cutoff was applied to exclude genes that did not meet the predefined significance criteria, retaining genes with stronger statistical support for the modeled temporal pattern.

After fitting the piecewise linear mixture model, genes matching the specified rise-and-fall pattern were selected to construct gene sets. Pathway enrichment analysis was then performed to determine whether these genes were significantly associated with specific biological pathways. This approach aimed to identify pathways exhibiting transient activation over time that conventional linear models might fail to detect.

### Validation

<img width="4434" height="2331" alt="그림9" src="https://github.com/user-attachments/assets/1e298492-1b84-4121-a639-73e8786a115c" />

To evaluate the biological relevance of the gene sets identified by the pivot-based linear mixture model, the enriched pathways were compared with findings reported in the literature, particularly in[Temporal transcriptomic landscape of postnatal mouse ovaries reveals dynamic gene signatures  associated with ovarian aging.pdf](https://github.com/user-attachments/files/33672241/Temporal.transcriptomic.landscape.of.postnatal.mouse.ovaries.reveals.dynamic.gene.signatures.associated.with.ovarian.aging.pdf). The reported findings were consistent with the pathway enrichment results obtained from the pivot-based gene sets, providing literature-based support for their biological relevance. This agreement suggests that the proposed model can identify biologically meaningful temporal expression patterns and supports its potential utility in detecting transient pathway-associated changes.

### code : https://github.com/AiJayce/Soongsil.Univ/blob/Uterus-Aging/Oval%20-%20Mouse/RNA-seq/class.ipynb

## (7) Discussion

### Data and Analysis Context

This project investigated age-associated transcriptional dynamics using publicly available bulk RNA-seq data collected across multiple developmental time points. Raw sequencing reads were processed and aligned to the reference genome, followed by count normalization and variance stabilization using DESeq2. Pearson correlation, cosine similarity, eigenvalue structure, and principal component analysis (PCA) were used to evaluate sample similarity and temporal organization before gene-level modeling.

Gene expression trajectories were initially characterized using linear regression. Genes were classified according to their regression slopes, distinguishing genes with increasing expression over time from those with decreasing expression. Genes with sufficiently low regression residuals were retained to prioritize expression patterns that were well represented by the fitted linear models. The resulting gene sets were subsequently evaluated using pathway enrichment analysis to identify biological processes associated with aging.

### Interpretation of Temporal Expression Patterns

The primary objective was to characterize aging-associated expression changes as temporal trajectories rather than as independent pairwise comparisons between age groups. This approach provides a direct representation of the direction and consistency of gene expression changes across time points.

Positive-slope gene sets captured genes whose expression increased over time, potentially reflecting progressive activation of age-associated biological programs. Conversely, negative-slope gene sets represented genes whose expression decreased over time, potentially indicating reduced activity of processes associated with younger developmental stages.

Pathway enrichment analysis identified biological processes involving immune regulation, inflammation, epigenetic regulation, fertilization, and cell-cycle activity. The recovery of biologically interpretable pathways provided an initial indication that the selected gene sets captured structured transcriptional patterns rather than arbitrary expression fluctuations. Comparisons with previously published studies further provided a means of evaluating whether the identified pathways were consistent with established biological observations.

### Nonlinear Temporal Dynamics and Pivot-Based Modeling

A single linear regression model cannot adequately represent all temporal expression patterns. In particular, some genes and pathways may undergo transient activation, increasing during an intermediate period before declining at later time points. Such non-monotonic patterns may be overlooked when genes are classified exclusively according to a single positive or negative slope.

To address this limitation, a pivot-based linear mixture model was developed to approximate rise-and-fall expression trajectories. A candidate pivot time point was selected to divide the trajectory into two segments: an increasing phase before the pivot and a decreasing phase after it. Separate linear regression models were fitted to the two segments, and predefined significance-based filtering was applied to prioritize genes that supported the expected temporal pattern.

The resulting gene sets were then subjected to pathway enrichment analysis. By identifying pathways associated with transient expression changes, this approach extends the initial linear framework and provides a way to investigate temporally restricted biological responses that may not be captured by monotonic models. Agreement between the enriched pathways and findings reported in the literature offers supporting evidence for their biological relevance, although independent validation is required to establish the robustness of the model.

### Limitations and Future Directions

Several considerations are important when interpreting these results. First, the reliability of temporal modeling depends on the number and distribution of time points, the number of biological replicates, and the consistency of sample processing. A limited number of time points can constrain the identification of complex trajectories, while sample-level confounding may produce apparent temporal trends unrelated to aging itself.

Second, normalization and variance stabilization improve the comparability of expression profiles but do not independently eliminate batch effects. Differences in sequencing depth, library composition, experimental batches, and other technical variables should therefore be evaluated explicitly. Regression residual thresholds and pivot selection criteria may also influence which genes are retained and which pathways are subsequently enriched.

Finally, pathway enrichment indicates that selected genes are overrepresented in annotated biological processes; it does not directly establish pathway activation, causal relationships, or functional changes. Future work should evaluate model stability across alternative pivot positions and filtering thresholds, incorporate independent time-series datasets, and compare the proposed approach with nonlinear methods such as generalized additive models and spline-based regression. These analyses would help determine whether the identified temporal patterns are reproducible and whether the pivot-based model provides additional biological insight beyond conventional linear approaches.

