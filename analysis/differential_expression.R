# Subset the metadata and expression data to the two secretory epithelial cell types

secretorymeta<- subset(metadata, metadata$Celltype %in%c('Secretory Epithelial-1','Secretory Epithelial-2'))
secretoryexp<- subset(exp_data, rownames(exp_data) %in% secretorymeta$CellID)

# Sanity check that both have same number of cols. Actually easy to visualise in RStudio.
if (dim(secretorymeta)[1] == dim(secretoryexp)[1]) {
  print("same")
} else {
  print("not same")
}

# transforming the dataframe to rows are cols and vice versa which is required for the following steps
# a design matrix has observations in rows and predictors in columns
secretoryexp <- t(secretoryexp)
# sanity check
identical(colnames(secretoryexp), secretorymeta$CellID)

# differential gene expression test
# create design matrix using transformed dataset
design <- model.matrix(~ Celltype, data = secretorymeta)
# model the expression of gene expression across cell types
fit <- lmFit(secretoryexp, design)
fit <- eBayes(fit, trend = TRUE)

# extract the top 20 genes. These are the most variable expression levels.
top20genes <- topTable(fit, coef = 2, number = 20)

# do the same again to take the top 100
top_genes_df <- topTable(fit, coef = 2, number = 100)
# create a new column called significance, fill column per gene as significant or not significant based on p value and logFC
top_genes_df$significance <- ifelse(top_genes_df$P.Value < 0.05 & abs(top_genes_df$logFC) > 1,
                                    "Significant", "Not Significant")

# Plot all 100 genes based on significance and log2 fold change. 
# Highlight and label those labelled "significant" in the above step.
volcanoplot <- ggplot(top_genes_df, aes(x = logFC, y = -log10(P.Value), color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change", y = "-Log10(p-value)") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top") +
  geom_text_repel(aes(label = ifelse(significance == "Significant", rownames(top_genes_df), "")),
                  box.padding = 0.5, point.padding = 0.5,
                  max.overlaps = 10, size = 3)

volcanoplot
