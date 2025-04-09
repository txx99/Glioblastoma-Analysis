library(tidyverse)
library(ggplot2)
library(DESeq2)
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO#####

gset <- getGEO("GSE231994", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL27956", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


summary(gset)
dim(gset)
dim(df)

expr_data <- exprs(gset)
pheno_data <- gset@phenoData@data
expr_df <- as.data.frame(expr_data)

head(expr_df)
head(pheno_data)



expr_data_log <- log2(expr_data +1)

diagnosis <- pheno_data$characteristics_ch1.1
diagnosis <- gsub("disease state ","", diagnosis)
diagnosis <- gsub(" ", "_", diagnosis) 
diagnosis <- make.names(diagnosis)


table(diagnosis)

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




#log2 transform





# # box-and-whisker plot
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




###PROCESSED DATA####
df <- read.table('GSE232050_processed_data.txt', 
                 sep = '\t',
                 header = TRUE)



