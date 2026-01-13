install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "GEOquery",
  "apeglm",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA"
), ask = FALSE, update = FALSE)

install.packages(c(
  "tidyverse",
  "caret",
  "pheatmap",
  "randomForest",
  "ggplot2"
))
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(caret)
library(randomForest)
library(apeglm)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(biomaRt)
counts_raw <- read.csv(
  "C:/Users/Admin/Downloads/GSE288479_raw_counts_All_samples.txt.gz",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(counts_raw) <- counts_raw$gene_id
counts <- counts_raw[, grep("_counts$", colnames(counts_raw))]
counts[] <- lapply(counts, as.numeric)
keep <- rowSums(counts >= 10) >= (0.2 * ncol(counts))
counts_filtered <- counts[keep, ]
metadata <- data.frame(
  sample = colnames(counts_filtered),
  condition = ifelse(
    grepl("_S_counts$", colnames(counts_filtered)),
    "Tumor",
    "Normal"
  )
)
rownames(metadata) <- metadata$sample
metadata$condition <- factor(metadata$condition)
metadata$condition <- relevel(metadata$condition, ref = "Normal")
counts_filtered <- as.matrix(counts_filtered)
storage.mode(counts_filtered) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res <- res[order(res$padj), ]
res_shrunk <- lfcShrink(
  dds,
  coef = "condition_Tumor_vs_Normal",
  type = "apeglm"
)
volcano_df <- as.data.frame(res_shrunk)
volcano_df$gene <- rownames(volcano_df)

volcano_df$significance <- "Not significant"
volcano_df$significance[
  volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 1
] <- "Upregulated"
volcano_df$significance[
  volcano_df$padj < 0.05 & volcano_df$log2FoldChange < -1
] <- "Downregulated"

ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not significant" = "grey"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal()
vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = "condition")

sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists))
res_df <- as.data.frame(res_shrunk)

feature_genes <- rownames(res_df)[
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1
]

X <- t(assay(vsd)[feature_genes, ])
y <- factor(metadata$condition)

ml_data <- data.frame(X)
ml_data$condition <- y
ctrl <- trainControl(
  method = "LOOCV",
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

cv_model <- train(
  condition ~ .,
  data = ml_data,
  method = "glm",
  family = binomial,
  metric = "ROC",
  trControl = ctrl
)

cv_model
cv_model$results$ROC
rf_fit <- randomForest(
  condition ~ .,
  data = ml_data,
  importance = TRUE
)

varImpPlot(rf_fit)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = sub("\\..*", "", rownames(res_shrunk)),
  mart = mart
)

res_df$ensembl_gene_id <- sub("\\..*", "", rownames(res_df))

res_annotated <- merge(
  res_df,
  gene_map,
  by.x = "ensembl_gene_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)
sig_genes <- subset(
  res_annotated,
  padj < 0.05 & abs(log2FoldChange) > 1
)

gene_entrez <- bitr(
  sig_genes$hgnc_symbol,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

entrez_ids <- unique(gene_entrez$ENTREZID)

ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
)

ekegg <- enrichKEGG(
  gene     = entrez_ids,
  organism = "hsa"
)

ereact <- enrichPathway(
  gene     = entrez_ids,
  organism = "human",
  readable = TRUE
)

dotplot(ego, showCategory = 10)
dotplot(ekegg, showCategory = 10)
dotplot(ereact, showCategory = 10)
gene_list <- res$log2FoldChange
names(gene_list) <- sub("\\..*", "", rownames(res))
gene_list <- sort(gene_list, decreasing = TRUE)

gene_map <- bitr(
  names(gene_list),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_list_entrez <- gene_list[gene_map$ENSEMBL]
names(gene_list_entrez) <- gene_map$ENTREZID

gsea_go <- gseGO(
  geneList = gene_list_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pvalueCutoff = 0.05
)
