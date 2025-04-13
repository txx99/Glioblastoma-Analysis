#Analysis of nCounter Data
#Prior to analysis please pass RCC data through Rosalind to perform normalization and differential analysis
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(org.Hs.eg.db)

setwd("Your Working Directory")

#Read in DE results from Rosalind
# res <- read.csv("Volcano Plot CSV from Rosalind.csv") 
# DE analysis was performed in Rosalind and data can be found in GEO GSE231994
res <- read.csv("DE_results_PD_vs_psPD.csv")
rownames(res) <- res$Symbol
head(res)



#Gene ontology
gene_list <- as.numeric(res$logFC) 
colnames(res)[1] <- "Symbol"  # changes "X" to "Symbol"

hs=org.Hs.eg.db
IDs <- AnnotationDbi::select(hs, 
                             keys = res$Symbol, 
                             columns = c("ENTREZID", "SYMBOL"), 
                             keytype = "SYMBOL")

res_merged <- merge(res, IDs, by.x = "Symbol", by.y = "SYMBOL")
# Create named vector
gene_list <- res_merged$logFC
names(gene_list) <- res_merged$ENTREZID

# Clean: remove NAs and duplicates
gene_list <- gene_list[!is.na(names(gene_list))]
gene_list <- gene_list[!duplicated(names(gene_list))]

# Sort descending (required by fgsea)
gene_list <- sort(gene_list, decreasing = TRUE)

# names(gene_list) <- IDs$ENTREZID
# gene_list = na.omit(sort(gene_list, decreasing = TRUE))

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = hs, 
             pAdjustMethod = "none")

#Visualize gene enrichment
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

# # res$logFC <- as.numeric(res$logFC)
# 
# EnhancedVolcano(res,
#                 lab = rownames(res),
#                 x = 'logFC',
#                 y = 'P.Value',
#                 xlab = bquote(~Log[2]~ 'fold change'),
#                 title = NULL,
#                 FCcutoff = 1,
#                 pCutoff = .05,
#                 drawConnectors = T,
#                 widthConnectors = 0.1,
#                 colConnectors = 'black', 
#                 labSize = 5, 
#                 pointSize = 1,
#                 legendPosition = 'none',
#                 subtitle = "",
#                 ylim = c(0, -log10(10e-11))) + scale_x_reverse()

