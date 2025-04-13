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
library(DOSE)

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
# Now we create a named ranked list/ vector for tools like gseGO() or fgsea(), 
# which expects a named numeric vector of log fold change.
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
# if colname appears twice, save only the first instance.

# Sort descending (required by fgsea and clusterProfiler for GSEA)
gene_list <- sort(gene_list, decreasing = TRUE)


# GENE ONTOLOGY!
# gseGO() randomly shuffles the ranked gene list over and over and calculate 
# the enrichment score of those randomly generated list.
# So if the actual enrichment score is much higher than those random gene lists, 
# then it means our results is statistically significant.
# However, because of that random permutation, the analysis can differ slightly 
# everytime we run the code.
# To ensure the consistency, we'll run set.seed() and set the seed parameter in 
# gseGO() to TRUE.
set.seed(1234)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL",              # Ontology: biological process (BP), cellular components (CC), molecular function (MF) or all 3 (ALL)
             keyType = "ENTREZID",    # Gene list is named by Entrez IDs
             pvalueCutoff = 0.05,     # Only keep GO terms with p < 0.05
             verbose = TRUE,          # Print progress in console
             OrgDb = hs,              # Organism database: human
             seed = TRUE,             # Makes the permutation process reproducible
             pAdjustMethod = "none")  # "none" for raw p-value or "BH" (Benjamin-Hochberg) for FDR correction.


#Visualize gene enrichment
dotplot(gse,                # GSEA results
        showCategory=10,    # Show top 10 most significantly enriched gene sets
        split=".sign") +    # Splits the plot into 2 frame, one activated/up-regulated and one suppressed/down-regulated.
                            # If we don't split, the plot only shows ranking of enrichment score or stat significance, 
                            # without accounting for whether the gene set is up or down-regulated.
                            # Then we can't see which pathways are more active in PD or psPD.
  facet_grid(.~.sign)       # facet_grid(row~column): on 1 row, split the plot into columns based on the signs (activated vs suppressed). 



## RATHER THAN ACTIVATED VS SUPPRESSED, want to make this psPD vs PD.
# Nguyen: So I think the author just cropped out the "activated" and "suppressed" 
# part to make it more comprehensible for readers. 
# During DEA using limma, we set the contrast to PD - psPD, which is to see where
# PD stands relative to psPD.
# So if a gene set is up-regulated in PD compared to psPD, the pathway is "activated",
# while down-regulated gene set in PD relative to psPD is "suppressed"
# To explain this in a different way:
# the "activated" pathways correspond to gene sets that are more upregulated in PD, 
# while "suppressed" pathway correspond to gene sets that are more upregulated in psPD.

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

