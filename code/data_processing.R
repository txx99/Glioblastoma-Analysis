library(tidyverse)
library(ggplot2)
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO
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



