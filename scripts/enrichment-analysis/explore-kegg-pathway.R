
#BiocManager::install("clusterProfiler")
#BiocManager::install("GSEABase")
#BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)

sig_ids_df <- read.csv(file = 'sig_genes_entrez_ids.csv')
sig_ids <-  unique(sig_ids_df[!is.na(sig_ids_df$entrezId),]$entrezId)

# KEGG Enrichment Analysis of a gene set. Given a vector of genes, this function will return the enrichment KEGG categories with FDR control.
keggEnrich<-enrichKEGG(gene= sig_ids,
                       organism= "hsa",
                       pAdjustMethod="BH",
                       pvalueCutoff = 0.05)

  
  
# found results
#hsa05203
#hsa04933
#hsa05131
#hsa05164
#hsa04110


browseKEGG(keggEnrich, 'hsa05203')
