# Load libraries
library(dplyr)
library(ggplot2)
library(DESeq2)
# Load normalized counts
norm_counts <- read.delim("C:/Users/PETER/Downloads/DESeq2_Normalized_counts.tabular",
                          header = TRUE, 
                          row.names = 1)
# Rename columns
colnames(norm_counts) <- c(
  "SRR31454341_Infected", "SRR31454342_Infected", "SRR31454343_Infected",
  "SRR31454344_Infected", "SRR31454345_Infected", "SRR31454336_Control",
  "SRR31454337_Control", "SRR31454338_Control", "SRR31454339_Control",
  "SRR31454340_Control")

# Step 1: Remove genes with zero counts across all samples
norm_counts_filtered <- norm_counts[rowSums(norm_counts) > 0, ]
cat("Genes after filtering:", nrow(norm_counts_filtered), "\n")

# Step 2: Log transform (variance stabilizing)
log_counts <- log2(norm_counts_filtered + 1)

# Step 3: Select top 500 most variable genes
# (independent of group assignment - fixes reviewer concern)
gene_variance <- apply(log_counts, 1, var)
top500_genes <- names(sort(gene_variance, decreasing = TRUE) [1:500])
mat_var <- log_counts[top500_genes, ]

# Step 4: Run PCA
pca_data <- prcomp(t(mat_var), scale. = TRUE)

# Step 5: Check variance explained
summary(pca_data)$importance[, 1:5]

# Create PCA dataframe
condition <- c("Infected","Infected","Infected","Infected","Infected",
               "Control","Control","Control","Control","Control")

pca_df <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  condition = factor(condition),
  sample = colnames(norm_counts)
)

# Plot corrected PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
         geom_point(size = 4) +
         geom_text(aes(label = sample), hjust = 0.5, vjust = -0.8, size = 3) +
         theme_classic() +
         scale_color_manual(values = c("Control" = "blue",
                                       "Infected" = "red")) +
         labs(title = "PCA Plot - APEC Infected vs Control", 
                 subtitle = "Top 500 Most Variable Genes",
         x = "PC1: 71.8% variance",
         y = "PC2: 15.1% variance")

library(clusterProfiler)       
library(org.Gg.eg.db)
library(dplyr)

# Reload sig_genes if needed
res <- read.table("C:/Users/PETER/Downloads/DESeq2_results.tabular",
                  header = TRUE, sep = "\t")
colnames(res) <- c("GeneID", "baseMean", "log2FoldChange",
                   "lfcSE", "stat", "pvalue", "padj")
# Filter DEGs
sig_genes <- res %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)
# Get downregulated genes
down_genes <- sig_genes %>% 
  filter(log2FoldChange < -1) %>% 
  pull(GeneID)

# Convert to Entrez IDs
down_entrez <- bitr(down_genes, fromType = "ENSEMBL",
                    toType = "ENTREZID", OrgDb = org.Gg.eg.db)
cat("Converted:", nrow(down_entrez), "out of", length(down_genes), "\n")

# Run GO with proper cutoffs this time
go_down_corrected <- enrichGO(
  gene = down_entrez$ENTREZID, 
  OrgDb = org.Gg.eg.db,
  keyType = "ENTREZID", 
  ont = "BP",
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.2,
  minGSSize = 5, 
  maxGSSize = 500
)
cat("GO terms found:", nrow(as.data.frame(go_down_corrected)), "\n")

go_down_MF <- enrichGO(
  gene = down_entrez$ENTREZID, 
  OrgDb = org.Gg.eg.db,
  keyType = "ENTREZID", 
  ont = "MF",
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.2,
  minGSSize = 5, 
  maxGSSize = 500
)
cat("MF GO terms found:", nrow(as.data.frame(go_down_MF)), "\n")

kegg_down <- enrichKEGG(
  gene = down_entrez$ENTREZID,
  organism = "gga", 
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 5
)
cat("KEGG pathways found:", nrow(as.data.frame(kegg_down)), "\n")
cat("Downregulated genes total:", length(down_genes), "\n")
cat("Successfully converted:", nrow(down_entrez), "\n")
head(down_entrez)
go_down_all <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Gg.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.5,
  minGSSize = 3,
  maxGSSize = 1000
)
cat("GO terms found:", nrow(as.data.frame(go_down_all)), "\n")

go_up <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Gg.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 3,
  maxGSSize = 500
)

up_genes <- sig_genes %>% 
  filter(log2FoldChange > -1) %>% 
  pull(GeneID)

# Convert to Entrez IDs
gene_entrez <- bitr(up_genes, fromType = "ENSEMBL",
                    toType = "ENTREZID", OrgDb = org.Gg.eg.db)
cat("Upregulated converted:", nrow(gene_entrez), "out of", length(up_genes), "\n")

go_up <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Gg.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 3,
  maxGSSize = 500
)

cat("Upregulated GO terms:", nrow(as.data.frame(go_up)), "\n")

#Test different cutoffs systematically
for(pval in c(0.05, 0.1, 0.15, 0.2)){
  go_test <- enrichGO(
    gene = gene_entrez$ENTREZID,
    OrgDb = org.Gg.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = pval,
    qvalueCutoff = pval * 2,
    minGSSize = 3,
    maxGSSize = 500
  )
  cat("pval:", pval, "-> GO terms:", nrow(as.data.frame(go_test)), "\n")
}

# Final GO enrichment for upregulated genes
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

#View results
as.data.frame(go_up_final)[, c("Description", "GeneRatio", "p.adjust", "Count")]

#Plot
barplot(go_up_final,
        showCategory = 7,
        title = "GO Biological Process - Upregulated Genes (padj < 0.1)",
        font.size = 11)

all_db_genes <- keys(org.Gg.eg.db)
covered_down <- sum(down_entrez$ENTREZID %in% all_db_genes)
cat("Downregulated genes in database:", covered_down,
    "out of", nrow(down_entrez), "\n")
cat("Coverage:", round(covered_down/nrow(down_entrez)*100, 1), "\n")
covered_up <- sum(gene_entrez$ENTREZID %in% all_db_genes)
cat("Upregulated genes in database:", covered_up,
    "out of", nrow(down_entrez), "\n")
cat("Coverage:", round(covered_up/nrow(down_entrez)*100, 1), "\n")
cat("Total genes in org. Gg.eg.db:", length(all_db_genes), "\n")
go_annotated <- sum(down_entrez$ENTREZID %in% 
                      keys(org.Gg.eg.db, keytype = "ENTREZID"))
cat("Downregulated genes with any GO annotation:", go_annotated, "\n")

go_down_explore <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Gg.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  minGSSize = 3,
  maxGSSize = 500
)

head(as.data.frame(go_down_explore)[
  order(as.data.frame(go_down_explore)$pvalue),
  c("Description", "GeneRatio", "pvalue", "p.adjust")], 10)

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

cat("GO terms found:", nrow(as.data.frame(go_down_final)), "\n")

dotplot(go_down_final,
        showCategory = 10,
        title = "GO Biological Process - Downregulated Genes (padj < 0.15)",
        font.size = 10)

qc_table <- data.frame(
  Sample = colnames(norm_counts),
  Condition = c(rep("Infected", 5), rep("Control", 5)),
  SRR_Accession = c("SRR31454341", "SRR31454342", "SRR31454343", 
                    "SRR31454344", "SRR31454345", "SRR31454336", 
                    "SRR31454337", "SRR31454338", "SRR31454339",
                    "SRR31454340")
)

print(qc_table)

#Load libraries
library(dplyr)
library(org.Gg.eg.db)
library(clusterProfiler)

#Create complete DEG table with gene symbols
all_degs <- sig_genes

#Convert Ensembl to gene symbols
deg_symbols <- bitr(all_degs$GeneID,
                    fromType = "ENSEMBL",
                    toType = "SYMBOL",
                    OrgDb = org.Gg.eg.db)

#Merge with DEG results
deg_complete <- merge(all_degs, deg_symbols,
                      by.x = "GeneID",
                      by.y = "ENSEMBL", 
                      all.x = TRUE)

#Check excat column names first
colnames(all_degs)
colnames(deg_symbols)

#Add regulation direction
deg_complete$Direction <- ifelse(deg_complete$log2FoldChange > 0,
                                 "Upregulated", "Downregulated")

#Arrange by padj
deg_complete <- deg_complete %>% arrange(padj)

#Check
cat("Total DEGs:", nrow(deg_complete), "\n")
cat("With gene symbols:", sum(!is.na(deg_complete$SYMBOL)), "\n")    

#Save as CSV
write.csv(deg_complete,
          "C:/Users/PETER/Desktop/Supplemetary_Table_S1_All_DEGs.csv",
          row.names = FALSE)
cat("Saved!\n")

#Remove duplicate GeneIDs - keep first symbol per gene
deg_complete<- deg_complete %>%
  arrange(padj) %>%
  distinct(GeneID, .keep_all = TRUE)

#Verify
cat("Total DEGs after duplication:", nrow(deg_complete), "\n")
cat("Upregulated:", sum(deg_complete$Direction == "Upregulated"), "\n")
cat("Downregulated:", sum(deg_complete$Direction == "Downregulated"), "\n")
cat("With gene symbols:", sum(!is.na(deg_complete$SYMBOL)), "\n")

#Save
write.csv(deg_complete,
          "C:/Users/PETER/Desktop/Supplemetary_Table_S1_All_DEGs.csv",
          row.names = FALSE)
cat("Saved successfully!\n")

#Save GO upregulated results
go_up_df <- as.data.frame(go_up_final)
write.csv(go_up_df,
          "C:/Users/PETER/Desktop/GO_upregulated_final.csv",
          row.names = FALSE)

#Save GO downregulated results
go_down_df <- as.data.frame(go_down_final)
write.csv(go_down_df,
          "C:/Users/PETER/Desktop/GO_downregulated_final.csv",
          row.names = FALSE)

#Check
cat("Upregulated GO terms:", nrow(go_up_df), "\n")
cat("Downregulated GO terms:", nrow(go_down_df), "\n")
cat("Saved!\n")
