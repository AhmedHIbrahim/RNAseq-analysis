library(clusterProfiler)
library(org.Hs.eg.db)


sig_ids_df <- read.csv(file = 'sig_genes_entrez_ids.csv')
sig_ids <-  unique(sig_ids_df[!is.na(sig_ids_df$entrezId),]$entrezId)


# GO Enrichment Analysis of a gene set. Given a vector of genes, this function will return the enrichment GO categories after FDR control.
goEnrich<-enrichGO(gene= sig_ids,
                   OrgDb= org.Hs.eg.db,
                   ont= "ALL",
                   pAdjustMethod="BH",
                   pvalueCutoff = 0.05,
                   readable= TRUE)



enrichGoAndKegg("sig-",sig_ids, 0.01)
enrichGoAndKegg("sig-",sig_ids, 0.05)