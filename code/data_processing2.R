library(tidyverse)           # Data wrangling (contains dplyr, ggplot2,...)
# library(ggplot2)           # For plotting, redundant as it's already in tidyverse
library(GEOquery)            # To download datasets from GEO
library(limma)               # DEA (linear models for microarrays but usable for RNA-seq)
# library(umap)              # UMAP (non-linear dimensionality reduction)
library(EnhancedVolcano)     # To make more beautiful volcano plots (not used below)
# library(clusterProfiler)   # Enrichment analysis (GO, KEGG)
# library(enrichplot)        # Plotting enrichment plot
# library(ggplot2)             
library(ggpubr)              # To make publication-ready plot
library(gridExtra)           # Combine ggplots
library(org.Hs.eg.db)        # Human gene annotations

# Load data set from GEO (not the platform info).
# In case the dataset is associated with multiple platforms, 
# we specify the GLP (platform).
gset <- getGEO("GSE231994", GSEMatrix=TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL27956", attr(gset, "names")) else idx <- 1 #if >1, idx = location of 'GPL27956' under gset 'names', otherwise idx=1 
gset <- gset[[idx]] # gset <- content of gset[[data folder of interest]]. basically removes one extraneous level of folders; if length(gset) > 1 , has multiple folders --> wanted the folder with GPL2795 in title. 

# 2. Extract expression and phenotype data
expr_data <- exprs(gset)  # Get expression matrix as df
pheno_data <- pData(gset) # Get phenotype data as df

# 3. Clean and verify group labels
diagnosis_df<- pheno_data %>% dplyr::select(characteristics_ch1.1) # chooses the characteristics we are looking at (diagnosis: PD - psPD).
diagnosis_df$characteristics_ch1.1 <- gsub("disease state: ", "", diagnosis_df$characteristics_ch1.1) #pattern, replacement, x ; delete all "disease state: " text from diagnosis observations
print(table(diagnosis_df)) # Verify groups (should show PD vs psPD)
diagnosis<-diagnosis_df$characteristics_ch1.1 # list for downstream use

# 4. Create design matrix WITHOUT INTERCEPT (each group gets its own column)
# First we will make a design matrix WITHOUT INTERCEPT (each group gets its own column). 
# The intercept is essentially the "reference" (or baseline) group, from which 
# we will base the comparisons on (example: evaluate how different group B is 
# compared to group A). In this case if we make a design matrix with intercept, 
# the 'PD' group will become the baseline and essentially "locks" the comparison 
# to only seeing the relative difference from 'psPD' to 'PD'.
# This will prevents any other down-stream customization for analysis 
# (psPD to PD, PD to psPD, psPD/2 to PD,...), hence why we need a design matrix WITHOUT intercept.
design <- model.matrix(~0 + diagnosis, diagnosis_df) # makes each diagnosis into a column, allocates patients to one or the other column
colnames(design) <- make.names(levels(factor(diagnosis))) # Ensure valid names; colnames as factor type (= categories) 
# Name the columns based on the values in 'diagnosis'.
# factor() to convert the character values into categorical values.
# levels() to extract the unique categorical values/ group labels.
# make.names() to clean up the values to make it appropriate to use as column names
# (remove white space, parenthesis,...).
# It assigns the column name order by alphabetical order.
# So this will work with the design matrix we made as the column is also assigned 
# based on alphabetical order.
# Or just this:
# colnames(design) <- c("PD", "psPD")
print(colnames(design))  # Should show [1] "PD"   "psPD"

# 5. LIMMA analysis
fit <- lmFit(expr_data, design)  # Linear model for large dataset.  ## is it cross matching the sample names?
contrast_matrix <- makeContrasts(PD - psPD, levels=design) # Defines the names of the parameters to be used downstream for contrast.
# Here we compare the average expression in PD, MINUS the average expression in psPD.
# AKA to compare the gene expression in PD relative to psPD.
# levels=design is telling the program to use the columns of the design matrix 
# as the basis for the contrast. 
# So the column names become the available LEVELS to use for comparison.
# We can also do multiple contrasts: 
# x = c("psPD-PD", "(psPD+PD)/2")
# contrast_matrix <- makeContrasts(contrast = x, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix) # Applies the contrast to the linear model.
# This will give us the log fold change, standard errors, t-statistics, and residual degrees of freedom.
fit2 <- eBayes(fit2)  ## Normalize the data use Bayesian distribution.
# This also give us the p-value.

# 6. Get results (top differentially expressed genes)
results <- topTable(fit2, number=Inf, adjust.method="BH")
sig_genes <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]

# 7. Save results
write.csv(results, "DE_results_PD_vs_psPD.csv")

## bASIC Volcano Plot
ggplot(results, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color=abs(logFC) > 1 & adj.P.Val < 0.05)) +
  scale_color_manual(values=c("gray","red")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed")









?getGEO

?makeContrasts

gset <- getGEO("GSE231994", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL27956", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#assayData is expression data (not count data), therefore use edgeR/limma for DEA
summary(gset)
dim(gset)
dim(df)
expr_data <- exprs(gset)
expr_df <- as.data.frame(expr_data)
head(expr_df)

# see disease groups being compared : PD + pseudo PD
pheno_data <- gset@phenoData@data
head(pheno_data)
diagnosis <- pheno_data$characteristics_ch1.1
diagnosis <- gsub("treatment: ","", diagnosis)
table(diagnosis)

# expr_data_log <- log2(expr_data +1)

# doing some modeling + whatnots
design <- model.matrix(~0 +diagnosis)
colnames(design) <- levels(factor(diagnosis))

fit <- lmFit(expr_data_log, design)
contrasts <- makeContrasts(
  contrasts ="disease_state._PD - disease_state._psPD",
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
results <- topTable(fit2, number = Inf, adjust.method = "BH")










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
expr_data <- log2(expr_data) } 
# and transform ex to log2(ex)
# if LogC were False, data distribution is already acceptable, no need to log2


# # box-and-whisker plot ----------------------------
# dev.new(width=3+ncol(gset)/6, height=5)
# par(mar=c(7,4,2,1))
# title <- paste ("GSE231994", "/", annotation(gset), sep ="")
# boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
# dev.off()
# 
# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE231994", "/", annotation(gset), " value distribution", sep ="")
plotDensities(expr_data, main=title, legend=F)
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
