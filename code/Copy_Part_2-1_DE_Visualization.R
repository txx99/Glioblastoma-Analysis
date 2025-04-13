# Differential gene expression analysis was completed using limma from GEO dataset (GSE231994).
# The DEA results used here are saved in 'DE_results_PD_vs_psPD.csv'.
# We now perform gene set enrichment analysis using clusterProfiler.
# We will specifically gseGO() to use Gene Ontology (GO) gene set.
# NOTE 1: In our class we use GSEA(), which is better if we're running enrichment on custom gene sets
# like Hallmark, MSigDB, KEGG,...
# NOTE 2: gseGO() for GO enrichment analysis uses 

#=============== Library =================
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(org.Hs.eg.db)

# ========== Work directory ==============
#setwd("Your Working Directory")

# ========= Load DEA results =============
#Read in DE results
# Source code used NanoString Rosalind for DE, but we did it manually with {limma}
res <- read.csv("DE_results_PD_vs_psPD.csv")
colnames(res)[1] <- "Symbol"  # changes colname of gene names from "X" to "Symbol"
rownames(res) <- res$Symbol # assigns the gene names as the proper rownames not just elements of the first column
head(res)


#===Gene ontology==========================================================
#finding EntrezID from database and appending it to DE df
# In our results df, our genes are recorded as gene symbols (like TP53, GAPDH)
# But downstream tools like gseGO() expect standardized IDs like EntrezID.
# So we need to map the gene symbols to entrez ID.
# First access the human gene annotation database from package org.Hs.eg.db
# This is a massive "table" containing gene symbols, Entrz ID, Ensembl ID,
# GO terms, and more.
hs=org.Hs.eg.db
# Extract the relevant EntrezID for each of our gene symbol.
IDs <- AnnotationDbi::select(hs, 
                             keys = res$Symbol, 
                             # provide a list of genes to search up.
                             columns = c("ENTREZID", "SYMBOL"), 
                             # save EZID and Gene Symbol for the matching keys in unique columns.
                             keytype = "SYMBOL") 
# tells the code what the 'type' of the keys is, ex. gene Symbol, RefSeq ID, Accession Number, etc
# In this case, it is gene SYMBOL (or names).
# Now joins the df together by gene symbols. 
# Note that the gene symbol column is named differently in the 2 df 
# (res$Symbol and IDs$SYMBOL) so we need both by.x and by.y.
res_merged <- merge(res, IDs, by.x = "Symbol", by.y = "SYMBOL") 
rownames(res_merged) <- res_merged$Symbol

# Create named vector-----------------------------------------------
# Now we create a ranked list
#gene_list <- dplyr::select(res_merged, logFC) #nvm this makes df of logFC values vs rownames(Symbol)
gene_list<- res_merged$logFC 
# makes a vector of logFC values
names(gene_list) <- res_merged$ENTREZID 
#assigns EZID as colname to each value, making it a named list

# Clean: -----------------------------------------------------------
# remove NAs and duplicates 
gene_list <- gene_list[!is.na(names(gene_list))]
# if colname in gene_list is NA, do not save it
gene_list <- gene_list[!duplicated(names(gene_list))]
# if colname appears twice, do not save it (do not save either of them...?)

# Sort descending (required by fgsea)
gene_list <- sort(gene_list, decreasing = TRUE)

# names(gene_list) <- IDs$ENTREZID
# gene_list = na.omit(sort(gene_list, decreasing = TRUE))
set.seed(1234)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = hs, 
             seed = TRUE,
             pAdjustMethod = "none")

#Visualize gene enrichment
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
## RATHER THAN ACTIVATED VS SUPPRESSED, want to make this psPD vs PD.=======
# and also can adjust scale to be 0-0.8 like in the source image 



#Analysis of nCounter Data (NanoString Rosalind platform)
#nCounter = NanoString's modified system of DNA microarray technology.
#"Prior to analysis please pass RCC data through Rosalind to perform normalization and differential analysis"




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

