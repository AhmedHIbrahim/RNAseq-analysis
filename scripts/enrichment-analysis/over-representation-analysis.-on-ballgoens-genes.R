#________________GO_Classification_____________________#
library(org.Hs.eg.db)
library(clusterProfiler)
library("topGO")

# Ballgown's top significant genes
geneList <-scan('ballgown-genes.txt', character(), quote = "")


# convering gene symbol to gene ENTREZ ID
gene.df <- bitr(geneList, fromType = "SYMBOL",
                toType = c("ENTREZID" ),
                OrgDb = org.Hs.eg.db)

# GO classification,  ?groupGO <- help
## based on molecular function
ggo_mf <- groupGO(gene     = gene.df$ENTREZID,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "MF",            
                  level    = 3,
                  readable = TRUE)

## based on biological-process
ggo_bp <- groupGO(gene     = gene.df$ENTREZID,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "BP",            
                  level    = 3,
                  readable = TRUE)

## based on cellular-component
ggo_cc <- groupGO(gene     = gene.df$ENTREZID,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "CC",           
                  level    = 3,
                  readable = TRUE)


ggo_cc_df <- data.frame(ggo_cc)
ggo_mf_df <- data.frame(ggo_mf)
ggo_bp_df <- data.frame(ggo_bp)

ggo_cc_df = ggo_cc_df[ggo_cc_df$Count > 0, ]
ggo_mf_df = ggo_mf_df[ggo_mf_df$Count > 0, ]
ggo_bp_df = ggo_bp_df[ggo_bp_df$Count > 0, ]

#write.csv(ggo_cc_df, "ora-go-classification-based-on-cellular-component.csv")
#write.csv(ggo_mf_df, "ora-go-classification-based-on-molecular function.csv")
#write.csv(ggo_bp_df, "ora-go-classification-based-on-biological-process.csv")



## GO over-representation analysis
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)

ego_df = data.frame(ego)


#write.csv(ego_df, "ballgown-based-ora-go.csv")



# KEGG
# KEGG Enrichment Analysis of a gene set. Given a vector of genes, this function will return the enrichment KEGG categories with FDR control.
keggEnrich<-enrichKEGG(gene= gene.df$ENTREZID,
                       organism= "hsa",
                       pAdjustMethod="BH",
                       pvalueCutoff = 0.05)

keggDf<-data.frame(keggEnrich)
# browseKEGG(keggEnrich, 'hsa04330')
#write.csv(keggDf, "ballgown-based-ora-kegg.csv")

# A universal enrichment analyzer

c7 <- read.gmt("c7.all.v7.1.entrez.gmt")
c6 <- read.gmt("c6.all.v7.1.entrez.gmt")

ImmunSigEnrichC7 <- enricher(gene.df$ENTREZID,     # a vector of gene id
                           TERM2GENE=c7,         # user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
                           pvalueCutoff = 0.05)  # adjusted pvalue cutoff on enrichment tests to report

oncogSigEnricC6 <- enricher(gene.df$ENTREZID,     # a vector of gene id
                             TERM2GENE=c6,         # user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
                             pvalueCutoff = 0.05)  # adjusted pvalue cutoff on enrichment tests to report

# mapping geneID to gene Symbol
ImmunSigEnrichC7 <- setReadable(ImmunSigEnrichC7,
                              OrgDb = org.Hs.eg.db,
                              keyType = "ENTREZID")

oncogSigEnrichC6 <- setReadable(oncogSigEnricC6,
                                OrgDb = org.Hs.eg.db,
                                keyType = "ENTREZID")


write.csv(data.frame(ImmunSigEnrichC7), "ballgown-based-ora-immune.csv")
write.csv(data.frame(oncogSigEnrichC6), "ballgown-based-ora-oncog.csv")

