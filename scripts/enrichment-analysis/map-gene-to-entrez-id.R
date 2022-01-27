if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("mygene")

# high significant genes list from stringTie::ballgown and htseq_count::edgeR analysis

sig_file="sig-intersection-between-ballgown-and-edgeR-genes.txt"
significant_genes <- scan(sig_file, character(), quote = "")

# Few of the Annovar's found variants
annotated_genes_file="annovar-annotated-genes.txt"
annotated_genes <-unique(scan(annotated_genes_file, character(), quote = ""))

# get the corrosponding gene's entrez_accessions
sig_id_df <- mygene::queryMany(significant_genes, scopes="symbol", fields="entrezgene", species="human")
anno_id_df <- mygene::queryMany(annotated_genes, scopes="symbol", fields="entrezgene", species="human")


as.data.frame(sig_id_df)
str(sig_id_df)
write.csv(data.frame(gene=sig_id_df$query, entrezId=sig_id_df$entrezgene),
          "sig_genes_entrez_ids.csv", row.names = FALSE)


write.csv(data.frame(gene=anno_id_df$query, entrezId=anno_id_df$entrezgene),
          "anno_genes_entrez_ids.csv", row.names = FALSE)
