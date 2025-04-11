library(tidyverse)
library(ggplot2)
library(GEOquery)
library(limma)
library(umap)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(org.Hs.eg.db)


###Original Processed DATA from Paper Git hub ######

df <- read.table('GSE232050_processed_data.txt', 
                 sep = '\t',
                 header = TRUE)

#### download series and platform data from GEO#####

gset <- getGEO("GSE231994", GSEMatrix=TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL27956", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#### Extract expression and phenotype data####

expr_data <- exprs(gset)  # Get expression matrix from gset (measured expression values of different genes and samples)
pheno_data <- pData(gset) # Get phenotype data from gset

#### FIND diseased group and store count results in a table####
diagnosis <- pheno_data$characteristics_ch1.1  # chooses the characteristics we are looking at, disease state
diagnosis <- gsub("disease state: ", "", diagnosis)
print(table(diagnosis))  # Verify groups (should show PD vs psPD)

# 4. Create design matrix for LIMMA Differential Expression Analysis
design <- model.matrix(~0 + diagnosis)
colnames(design) <- make.names(levels(factor(diagnosis))) # Ensure valid names
print(colnames(design))  # Should show [1] "PD"   "psPD"

# 5. LIMMA analysis
fit <- lmFit(expr_data, design)
contrast_matrix <- makeContrasts(PD - psPD, levels=design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)  ##Normalize the data use Bayesian distribution

# 6. Get results (top differential expressed genes)
results <- topTable(fit2, number=Inf, adjust.method="BH")
sig_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ] ###filter the results for significant genes only

# 7. Save results
write.csv(results, "DE_results_PD_vs_psPD.csv")  ## save the filtered, significant results in a csv for visualization


## bASIC Volcano Plot
ggplot(results, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color=abs(logFC) > 1 & adj.P.Val < 0.05)) +
  scale_color_manual(values=c("gray","red")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed")





















