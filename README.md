# LUAD-transcriptomics-DESeq2
Lung Adenocarcinoma (LUAD) Transcriptomic Analysis
Overview
This project presents a bulk RNA-seq differential expression and downstream data mining analysis of lung adenocarcinoma (LUAD) samples. Using publicly available GEO data, we compare tumor and normal lung tissues to identify differentially expressed genes, characterize underlying biological pathways, and evaluate the predictive power of transcriptomic features for tumor classification.
The workflow integrates DESeq2-based differential expression, unsupervised visualization, machine learning classification, and functional enrichment analyses (GO, KEGG, Reactome, and GSEA) to provide biological and computational insights into LUAD progression.

Dataset
Source: GEO (GSE288479) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288479
Data type: Bulk RNA-seq raw gene count matrix
Samples: Tumor (solid component) vs Normal lung tissue
Organism: Homo sapiens

Analysis Workflow
1. Data Preprocessing
Imported raw count matrix
Filtered low-expression genes (≥10 counts in ≥20% of samples)
Constructed sample metadata (Tumor vs Normal)
Normalized counts using DESeq2

2. Differential Expression Analysis
Differential expression performed using DESeq2
Contrast: Tumor vs Normal
Log2 fold changes shrunk using apeglm
Significant genes defined by:
Adjusted p-value < 0.05
|log2FoldChange| ≥ 1

3. Exploratory Data Analysis
PCA: Assessed global separation between tumor and normal samples
Sample distance heatmap: Evaluated clustering consistency
Volcano plot: Visualized effect size vs statistical significance

4. Feature Construction
Selected biologically significant DE genes as features
Variance-stabilized counts (VST) used as model input

5. Predictive Modeling
Logistic regression (GLM) with Leave-One-Out Cross-Validation
Random Forest model for feature importance analysis
Model performance evaluated using ROC

6. Functional Enrichment Analysis
Gene annotation using biomaRt
Enrichment analyses:
Gene Ontology (Biological Process)
KEGG pathways
Reactome pathways
GSEA performed on ranked gene list (log2FC)

Key Results
Differential Expression
Hundreds of genes significantly dysregulated between tumor and normal tissue
Strongly upregulated genes include:
COL11A1, CA9, ITGA11, CHGA
These genes are associated with extracellular matrix remodeling, hypoxia response, and tumor invasiveness

Visualization Insights
PCA and heatmap demonstrate clear separation between tumor and normal samples
Volcano plot highlights genes with both strong effect size and statistical significance
Predictive Modeling
Logistic regression classifier achieved:
ROC ≈ 0.91

Indicates that transcriptomic features robustly distinguish tumor from normal tissue
Random Forest identified key discriminatory genes consistent with DESeq2 results
Functional Enrichment
Enriched biological processes include:
Hypoxia response
Extracellular matrix organization
Cell migration and adhesion
KEGG and Reactome pathways suggest activation of tumor progression and invasion mechanisms
GSEA confirms coordinated pathway-level dysregulation rather than isolated gene effects

Biological Interpretation
The analysis identifies genes and pathways that are differentially expressed at the RNA level, distinguishing LUAD tumor tissue from normal lung tissue. Importantly:
These results do not imply inherited cancer risk
They reflect tumor-associated transcriptional changes
Identified genes are potential biomarkers or therapeutic targets, not deterministic causes of cancer
Interpretation of Key Figures

1. PCA Plot (Tumor vs Normal)
What the plot shows
Principal Component Analysis on variance-stabilized gene expression
Samples cluster primarily by condition (Tumor vs Normal)

Interpretation

The clear separation between tumor and normal samples indicates that global transcriptomic differences exist between the two conditions
This suggests that LUAD tumors have a distinct expression profile, not driven by random noise or batch effects
The PCA validates that downstream differential expression and machine learning analyses are biologically meaningful

Why it matters
Confirms that tumor status is a dominant source of variation
Justifies supervised modeling and DE analysis

2. Sample Distance Heatmap
What the plot shows
Pairwise distances between samples based on normalized expression
Hierarchical clustering of samples

Interpretation
Tumor samples cluster together and separately from normal samples
Within-group similarity is higher than between-group similarity
Why it matters
Demonstrates internal consistency of tumor and normal groups
Indicates low technical noise and good sample quality
Supports the reliability of differential expression results

3. Volcano Plot (Shrunk log2FC vs −log10 adjusted p-value)
What the plot shows
Effect size (log2 fold change) vs statistical significance
Highlights significantly upregulated and downregulated genes

Results & Interpretation
1. Global Transcriptomic Structure
1.1 Principal Component Analysis (PCA)
Purpose
To assess whether tumor and normal samples differ at the global gene-expression level.
Key Observation
Samples segregate primarily by condition (Tumor vs Normal) along the first principal component.
Interpretation
This separation indicates that LUAD tumors possess a distinct transcriptional program compared to normal lung tissue. The dominant variance explained by tumor status confirms that biological differences—not technical noise—drive expression patterns.
Conclusion
The dataset is suitable for downstream differential expression and predictive modeling.

1.2 Sample-to-Sample Distance Heatmap
Purpose
To evaluate similarity between samples and confirm group consistency.
Key Observation
Tumor samples cluster together.
Normal samples form a separate cluster.
Interpretation
Higher within-group similarity and clear between-group separation demonstrate strong internal consistency and minimal batch effects.
Conclusion
Sample quality and normalization are reliable, supporting confidence in subsequent analyses.

2. Differential Gene Expression Analysis
2.1 Volcano Plot (DESeq2, Shrunk log2FC)
Purpose
To visualize genes with both large effect sizes and strong statistical support.
Key Observation
Numerous genes are significantly upregulated or downregulated in tumors.
Notable upregulated genes include COL11A1, CA9, ITGA11, and CHGA.
Interpretation
Upregulated genes are associated with hypoxia, extracellular matrix remodeling, and tumor invasion, while downregulated genes reflect loss of normal lung functions.
Conclusion
LUAD tumors exhibit coordinated dysregulation of genes linked to malignant progression.

2.2 Log2 Fold Change Shrinkage
Purpose
To stabilize effect size estimates for low-count or noisy genes.
Key Observation
Shrunk log2 fold changes preserve directionality but reduce extreme values.
Interpretation
Shrinkage improves reliability of gene ranking and downstream interpretation without altering biological conclusions.
Conclusion
Shrunk estimates are better suited for visualization, feature selection, and interpretation.

3. Functional and Pathway-Level Interpretation
3.1 Gene Ontology (GO) Biological Processes
Purpose
To identify biological processes overrepresented among differentially expressed genes.
Key Observation
Enriched processes include:
Hypoxia response
Cell migration
Extracellular matrix organization
Immune-related signaling
Interpretation
These processes are hallmark features of tumor progression and invasion in LUAD.
Conclusion
Differential expression reflects coordinated disruption of cancer-related biological programs.

3.2 KEGG Pathway Enrichment
Purpose
To link gene-level changes to known signaling pathways.
Key Observation
Enriched pathways involve:
Cancer signaling
Metabolic reprogramming
Cell–cell communication
Interpretation
Tumor samples activate pathways that promote survival, proliferation, and invasion.
Conclusion
Transcriptomic changes align with established LUAD molecular mechanisms.

3.3 Reactome Pathway Enrichment
Purpose
To provide fine-grained pathway interpretation.
Key Observation
Enrichment in pathways related to:
Extracellular matrix interactions
Inflammatory and immune signaling
Interpretation
Supports a tumor microenvironment characterized by stromal remodeling and immune interaction.
Conclusion
Pathway-level analysis strengthens biological relevance of DE findings.

3.4 Gene Set Enrichment Analysis (GSEA)
Purpose
To detect coordinated pathway shifts without relying on hard cutoffs.
Key Observation
Multiple cancer-related biological programs are enriched across the ranked gene list.
Interpretation
Tumor-associated transcriptional changes affect entire gene programs, not only top DEGs.
Conclusion
GSEA confirms robustness and system-level dysregulation in LUAD tumors.

4. Predictive Modeling and Feature Importance
4.1 Logistic Regression (LOOCV)
Purpose
To assess whether expression profiles can classify tumor vs normal samples.
Key Observation
ROC AUC ≈ 0.91
High sensitivity and specificity
Interpretation
A limited set of differentially expressed genes can accurately distinguish tumor from normal tissue.
Conclusion
Transcriptomic features are highly predictive of tumor status.

4.2 Random Forest Feature Importance
Purpose
To identify genes with the strongest predictive power.
Key Observation
A small subset of genes contributes disproportionately to classification accuracy.
Interpretation
These genes represent potential biomarkers and drivers of tumor-associated expression patterns.
Conclusion
Machine learning results corroborate DESeq2 findings and highlight candidate genes for further study.

5. Overall Biological Insight
Key Takeaways
LUAD tumors exhibit distinct global transcriptional profiles
Differential expression reflects coordinated cancer-related pathways
Identified genes differentiate tumor vs normal tissue at both statistical and predictive levels
Findings support known LUAD biology while highlighting candidate markers for progression
Important Note
These results identify tumor-associated expression changes, not inherited genetic risk. Differential expression reflects changes in tumor tissue, not germline mutations.
Limitations
Bulk RNA-seq does not resolve cell-type–specific expression
Small sample size limits generalizability
Findings are correlative and require experimental validation

Future Work
Integrate single-cell RNA-seq data to resolve tumor microenvironment
Validate biomarkers using independent LUAD cohorts
Incorporate multi-omics data (WES, methylation)
Develop more robust predictive models with external validation
Perform survival and clinical outcome association analysis
