#!/bin/bash

WORKSPACE=$HOME/workspace
HISAT2=$WORKSPACE/rnaseq/student_tools/hisat2-2.1.0
REF_DIR=$WORKSPACE/seq-analysis-course/ref
REF_GTF=$REF_DIR/GRCh38_latest_genomic.gff
REF_FASTA=$REF_DIR/GRCh38_latest_genomic.fna
REF_INDEX=$REF_DIR/hg38

$HISAT2/hisat2_extract_splice_sites.py $REF_GTF > $REF_DIR/splicesites.tsv
$HISAT2/hisat2_extract_exons.py $REF_GTF > $REF_DIR/exons.tsv


$HISAT2/hisat2-build -p 16\
                       --ss $REF_DIR/splicesites.tsv\
                       --exon $REF_DIR/exons.tsv\
                       $REF_FASTA\
                       $REF_INDEX


#$HISAT2/hisat2-build -p 8\
#                       --ss $REF_DIR/splicesites.tsv\
#                       --exon $REF_DIR/exons.tsv\
#                       $REF_FASTA\
#                       $REF_INDEX
