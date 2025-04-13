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

# load series and platform data from GEO
gset <- getGEO("GSE231994", GSEMatrix=TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL27956", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


#assayData is expression data (not count data), therefore use edgeR/limma for DEA
summary(gset)
dim(gset)
dim(df)
expr_data <- exprs(gset)
expr_df <- as.data.frame(expr_data)

hist(expr_data)
qqnorm(expr_data)


head(expr_df)

# see disease groups being compared : PD + pseudo PD
pheno_data <- gset@phenoData@data
head(pheno_data)
diagnosis <- pheno_data$characteristics_ch1.1
diagnosis <- gsub("treatment: ","", diagnosis)
table(diagnosis)


# log2 transform ----------------------------
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

# expr_data_log <- log2(expr_data +1)


# 2. Extract expression and phenotype data
expr_data <- exprs(gset)  # Get expression matrix
pheno_data <- pData(gset) # Get phenotype data

# 3. Clean and verify group labels
diagnosis <- pheno_data$characteristics_ch1.1  # choses the characteristics we are looking at, disease state
diagnosis <- gsub("disease state: ", "", diagnosis)
print(table(diagnosis))  # Verify groups (should show PD vs psPD)

# 4. Create design matrix
design <- model.matrix(~0 + diagnosis)
colnames(design) <- make.names(levels(factor(diagnosis))) # Ensure valid names
print(colnames(design))  # Should show [1] "PD"   "psPD"

# 5. LIMMA analysis
fit <- lmFit(expr_data, design)
contrast_matrix <- makeContrasts(PD - psPD, levels=design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)  ##Normalize the data use bayesian distribution

# 6. Get results (top differentially expressed genes)
results <- topTable(fit2, number=Inf, adjust.method="BH")
sig_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]

# 7. Save results
write.csv(results, "DE_results_PD_vs_psPD.csv")


## bASIC Volcano Plot
ggplot(results, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color=abs(logFC) > 1 & adj.P.Val < 0.05)) +
  scale_color_manual(values=c("gray","red")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed")




# do we need a log2 transform? ----------------------------
# run hist on expression data to see distribution
hist(expr_data)

qx <- as.numeric(quantile(expr_data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

# determining the quantile breaks in the ex data at 0% (min), 25% quantile, etc
# stored in qx as vector 
# na.rm=T remove na values

# Conditions for performing log2 transform
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
# if 5th item in qx (99th percentile) value > 100 OR 
# if 6th item in qx - qx[1] greater than 50 &  qx[2] is positive
if (LogC) { expr_data[which(expr_data <= 0)] <- NaN
# if LogC is true (either OR condition is met), then assign ex values <= 0 as NaN

ex<- log2(ex) } 

# and transform ex to log2(ex)
# if LogC were False, data distribution is already acceptable, no need to log2





# # box-and-whisker plot ----------------------------
# dev.new(width=3+ncol(gset)/6, height=5)
# par(mar=c(7,4,2,1))
# title <- paste ("GSE231994", "/", annotation(gset), sep ="")
# boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
# dev.off()
# 
# # expression value distribution plot
# par(mar=c(4,4,2,1))
# title <- paste ("GSE231994", "/", annotation(gset), " value distribution", sep ="")
# plotDensities(ex, main=title, legend=F)
# 
# # mean-variance trend
# ex <- na.omit(ex) # eliminate rows with NAs
# plotSA(lmFit(ex), main="Mean variance trend, GSE231994")
# 
# # UMAP plot (multi-dimensional scaling)
# ex <- ex[!duplicated(ex), ]  # remove duplicates
# ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
# plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
# library("maptools")  # point labels without overlaps
# pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)




###PROCESSED DATA###

df <- read.table('GSE232050_processed_data.txt', 
                 sep = '\t',
                 header = TRUE)




