#!/bin/bash

# Dir Ref
main_folder="$HOME/workspace/rnaseq"
workspace="$HOME/workspace/seq-analysis-course"

ALIGN_DIR="$workspace/alignment-data"
RESULT_DIR="$workspace/expression/htseq_counts"
REF_RNA_GTF="$workspace/ref/Homo_sapiens.GRCh38.99.gtf"


REPLICATES=("Ad-Rbm10_rep1" "Ad-Rbm10_rep2" "Ad-Rbm10_rep3"\
            "Ad-Gfb_rep1" "Ad-Gfb_rep2" "Ad-Gfb_rep3")


mkdir -p $RESULT_DIR
cd $RESULT_DIR

#1. calculate gene-level counts using HTSEQ

for REP in ${REPLICATES[@]}
do
   htseq-count --format bam\
               --order pos\
               --mode intersection-strict\
               --stranded reverse\
               --minaqual 1\
               --type exon\
               --idattr gene_id\
               $ALIGN_DIR/"$REP.bam"\
               $REF_RNA_GTF > $REP"_gene.tsv"

done


# 2. Merge results files into a single matrix for use in edgeR

join Ad-Rbm10_rep1_gene.tsv Ad-Rbm10_rep2_gene.tsv | join - Ad-Rbm10_rep3_gene.tsv | join - Ad-Gfb_rep1_gene.tsv | join - Ad-Gfb_rep2_gene.tsv | join - Ad-Gfb_rep3_gene.tsv > gene_read_counts_table_all.tsv

echo "GeneID Ad-Rbm10_rep1 Ad-Rbm10_rep2 Ad-Rbm10_rep3 Ad-Gfb_rep1 Ad-Gfb_rep2 Ad-Gfb_rep3" > header.txt

cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > gene_read_counts_table_all_final.tsv

rm -f gene_read_counts_table_all.tsv header.txt

head gene_read_counts_table_all_final.tsv

exit 0
