#!/bin/bash

# Dir Ref
main_folder="$HOME/workspace/rnaseq"
workspace="$HOME/workspace/seq-analysis-course"


RESULT_DIR="$workspace/expression/stringtie/ref_only"
REF_RNA_GTF="$workspace/ref/Homo_sapiens.GRCh38.99.gtf"

# Script Ref
stringtie_expression_matrix="$workspace/scripts/stringtie_expression_matrix.pl"


RESULT_DIRS="Ad-Rbm10_rep1,Ad-Rbm10_rep2,Ad-Rbm10_rep3,Ad-Gfb_rep1,Ad-Gfb_rep2,Ad-Gfb_rep3"
EXPRESSION_METRIC=("TPM" "FPKM" "Coverage")


cd $RESULT_DIR

for METRIC in ${EXPRESSION_METRIC[@]}
do
   $stringtie_expression_matrix\
         --expression_metric=$METRIC\
         --result_dirs=$RESULT_DIRS\
         --transcript_matrix_file="transcript_${METRIC,,}_all_samples.tsv"\
         --gene_matrix_file="gene_${METRIC,,}_all_samples.tsv"

done

exit 0
