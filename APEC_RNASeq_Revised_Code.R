# APEC Chicken RNA-Seq Analysis Pipeline — Revised Code
# Dataset: GSE282645 | Organism: Gallus gallus | Tissue: Liver | Day 5 dpi
# Revision: Updated PCA methodology and GO enrichment cutoffs
# Author: Peter Owolabi

# SECTION 1: LOAD LIBRARIES

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "EnhancedVolcano", "pheatmap",
                       "clusterProfiler", "org.Gg.eg.db"))
install.packages(c("dplyr", "ggplot2"))

library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Gg.eg.db)

# SECTION 2: LOAD AND CLEAN DESEQ2 RESULTS

# Load DESeq2 results from Galaxy
res <- read.table("C:/Users/PETER/Downloads/DESeq2_results.tabular",
                  header = TRUE,
                  sep = "\t")

# Rename columns
colnames(res) <- c("GeneID", "baseMean", "log2FoldChange",
                   "lfcSE", "stat", "pvalue", "padj")

# Verify
head(res)
colnames(res)

# SECTION 3: FILTER SIGNIFICANT DEGs

# Filter: padj < 0.05 AND |log2FC| > 1
sig_genes <- res %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)

# Summary
cat("Total significant DEGs:", nrow(sig_genes), "\n")
cat("Upregulated:", nrow(filter(sig_genes, log2FoldChange > 1)), "\n")
cat("Downregulated:", nrow(filter(sig_genes, log2FoldChange < -1)), "\n")

# Save all results
write.csv(sig_genes,
          "C:/Users/PETER/Desktop/significant_DEGs.csv",
          row.names = FALSE)

# SECTION 4: LOAD AND RENAME NORMALIZED COUNTS

# Load normalized counts from Galaxy DESeq2
norm_counts <- read.delim("C:/Users/PETER/Downloads/DESeq2_Normalized_counts.tabular",
                          header = TRUE,
                          row.names = 1)

# Rename columns (dataset order: 188,189,190,191,192,185,23,193,186,187)
# Matched to SRR numbers from Galaxy history
colnames(norm_counts) <- c(
  "SRR31454341_Infected",
  "SRR31454342_Infected",
  "SRR31454343_Infected",
  "SRR31454344_Infected",
  "SRR31454345_Infected",
  "SRR31454336_Control",
  "SRR31454337_Control",
  "SRR31454338_Control",
  "SRR31454339_Control",
  "SRR31454340_Control"
)

# Create sample annotation
condition <- c("Infected","Infected","Infected","Infected","Infected",
               "Control","Control","Control","Control","Control")

sample_annotation <- data.frame(
  condition = factor(condition),
  row.names = colnames(norm_counts)
)

# Define annotation colors
ann_colors <- list(
  condition = c(Control = "blue", Infected = "red")
)

# SECTION 5: VOLCANO PLOT

EnhancedVolcano(res,
                lab = res$GeneID,
                x = "log2FoldChange",
                y = "padj",
                title = "APEC Infected vs Non-infected",
                subtitle = "Chicken Liver - Day 5 Post Inoculation",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3,
                col = c("grey30", "forestgreen", "royalblue", "red2"))

# Save
ggsave("C:/Users/PETER/Desktop/volcano_plot.png",
       width = 12, height = 10, dpi = 300)

# SECTION 6: HEATMAP (TOP 50 DEGs)

# Get top 50 gene IDs
top50_genes <- head(sig_genes$GeneID, 50)

# Extract normalized counts for top 50
mat <- norm_counts[top50_genes, ]

# Z-score normalization
mat_scaled <- t(scale(t(mat)))

# Plot heatmap
pheatmap(mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = sample_annotation,
         annotation_colors = ann_colors,
         main = "Top 50 DEGs - APEC Infected vs Control",
         color = colorRampPalette(c("blue", "white", "red"))(100))

# Save using RStudio: Plots pane -> Export -> Save as Image (300 dpi)

# SECTION 7: CORRECTED PCA PLOT
# NOTE: PCA performed on top 500 most variable genes independent of
# differential expression status to avoid circularity bias
# (Revised per reviewer comment)

# Step 1: Remove zero-count genes
norm_counts_filtered <- norm_counts[rowSums(norm_counts) > 0, ]
cat("Genes after filtering:", nrow(norm_counts_filtered), "\n")

# Step 2: Log2 transform
log_counts <- log2(norm_counts_filtered + 1)

# Step 3: Select top 500 most variable genes
# (independent of group assignment)
gene_variance <- apply(log_counts, 1, var)
top500_genes <- names(sort(gene_variance, decreasing = TRUE)[1:500])
mat_var <- log_counts[top500_genes, ]

# Step 4: Run PCA with scaling
pca_data <- prcomp(t(mat_var), scale. = TRUE)

# Step 5: Check variance explained
summary(pca_data)$importance[, 1:5]

# Step 6: Create PCA dataframe
pca_df <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  condition = factor(condition),
  sample = colnames(norm_counts)
)

# Step 7: Plot corrected PCA
ggplot(pca_df, aes(x = PC1, y = PC2,
                   color = condition,
                   label = sample)) +
  geom_point(size = 4) +
  geom_text(aes(label = sample), hjust = 0.5, vjust = -0.8, size = 3) +
  theme_classic() +
  scale_color_manual(values = c("Control" = "blue",
                                "Infected" = "red")) +
  labs(title = "PCA Plot - APEC Infected vs Control",
       subtitle = "Top 500 Most Variable Genes (log2 normalized counts)",
       x = "PC1: 71.8% variance",
       y = "PC2: 15.1% variance")

# Save
ggsave("C:/Users/PETER/Desktop/PCA_plot_corrected.png",
       width = 10, height = 8, dpi = 300)

# SECTION 8: CORRECTED GO ENRICHMENT ANALYSIS
# NOTE: GO enrichment rerun with proper BH-corrected p-value thresholds
# Upregulated: padj < 0.1 | Downregulated: padj < 0.15 (nominal)
# (Revised per reviewer comment)

# --- 8.1 Prepare gene lists ---

up_genes <- sig_genes %>%
  filter(log2FoldChange > 1) %>%
  pull(GeneID)

down_genes <- sig_genes %>%
  filter(log2FoldChange < -1) %>%
  pull(GeneID)

cat("Upregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")

# --- 8.2 Convert Ensembl to Entrez IDs ---

gene_entrez <- bitr(up_genes,
                    fromType = "ENSEMBL",
                    toType = "ENTREZID",
                    OrgDb = org.Gg.eg.db)
cat("Upregulated converted:", nrow(gene_entrez),
    "out of", length(up_genes), "\n")

down_entrez <- bitr(down_genes,
                    fromType = "ENSEMBL",
                    toType = "ENTREZID",
                    OrgDb = org.Gg.eg.db)
cat("Downregulated converted:", nrow(down_entrez),
    "out of", length(down_genes), "\n")

# --- 8.3 GO enrichment - Upregulated (padj < 0.1) ---

go_up_final <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Gg.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2,
  minGSSize = 3,
  maxGSSize = 500
)

cat("Upregulated GO terms (padj < 0.1):",
    nrow(as.data.frame(go_up_final)), "\n")

# Visualize
barplot(go_up_final,
        showCategory = 7,
        title = "GO Biological Process - Upregulated Genes (padj < 0.1)",
        font.size = 11)

ggsave("C:/Users/PETER/Desktop/GO_upregulated_corrected.png",
       width = 10, height = 8, dpi = 300)

# Save results
write.csv(as.data.frame(go_up_final),
          "C:/Users/PETER/Desktop/GO_upregulated_final.csv",
          row.names = FALSE)

# --- 8.4 GO enrichment - Downregulated (nominal padj < 0.15) ---
# Note: No terms reached padj < 0.05 after BH correction
# This reflects broad non-specific suppression across diverse pathways

go_down_final <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Gg.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.15,
  qvalueCutoff = 0.15,
  minGSSize = 3,
  maxGSSize = 500
)

cat("Downregulated GO terms (nominal padj < 0.15):",
    nrow(as.data.frame(go_down_final)), "\n")

# Visualize
dotplot(go_down_final,
        showCategory = 10,
        title = "GO Biological Process - Downregulated Genes (padj < 0.15)",
        font.size = 10)

ggsave("C:/Users/PETER/Desktop/GO_downregulated_corrected.png",
       width = 10, height = 8, dpi = 300)

# Save results
write.csv(as.data.frame(go_down_final),
          "C:/Users/PETER/Desktop/GO_downregulated_final.csv",
          row.names = FALSE)

# SECTION 9: SUPPLEMENTARY TABLE S1 — COMPLETE DEG LIST

# Get gene symbols for all DEGs
deg_symbols <- bitr(sig_genes$GeneID,
                    fromType = "ENSEMBL",
                    toType = "SYMBOL",
                    OrgDb = org.Gg.eg.db)

# Merge with DEG results
deg_complete <- merge(sig_genes, deg_symbols,
                      by.x = "GeneID",
                      by.y = "ENSEMBL",
                      all.x = TRUE)

# Remove duplicates
deg_complete <- deg_complete %>%
  arrange(padj) %>%
  distinct(GeneID, .keep_all = TRUE)

# Add direction column
deg_complete$Direction <- ifelse(deg_complete$log2FoldChange > 0,
                                 "Upregulated", "Downregulated")

# Verify
cat("Total DEGs:", nrow(deg_complete), "\n")
cat("With gene symbols:", sum(!is.na(deg_complete$SYMBOL)), "\n")

# Save
write.csv(deg_complete,
          "C:/Users/PETER/Desktop/Supplementary_Table_S1_All_DEGs.csv",
          row.names = FALSE)

cat("All files saved successfully!\n")

# SESSION INFO

sessionInfo()

# End of script
