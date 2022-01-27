#!/bin/bash



# Dir Ref
main_folder="$HOME/workspace/rnaseq"
workspace="$HOME/workspace/seq-analysis-course"

OUT_DIR="$workspace/expression/stringtie/ref_only"
ALIGN_DIR="$workspace/alignment-data"
REF_RNA_GTF="$workspace/ref/Homo_sapiens.GRCh38.99.gtf"

mkdir -p $OUT_DIR

Replicates=("Ad-Rbm10_rep1" "Ad-Rbm10_rep2" "Ad-Rbm10_rep3"\
            "Ad-Gfb_rep1" "Ad-Gfb_rep2" "Ad-Gfb_rep3")

for REP in ${Replicates[@]}
do
   stringtie -p 16\
             -G $REF_RNA_GTF\
             -e\
             -B\
             -o $OUT_DIR/$REP/transcripts.gtf\
             -A $OUT_DIR/$REP/gene_abundances.tsv\
             $ALIGN_DIR/$REP.bam
done

exit 0
