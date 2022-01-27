#!/bin/bash


# Dir Ref
main_folder="$HOME/workspace/rnaseq"
workspace="$HOME/workspace/seq-analysis-course"

IN_DIR="$workspace/alignment-data"
OUT_DIR="$workspace/post-qc"


# Tools Ref
samtools=$main_folder/student_tools/samtools-1.9/samtools
fastqc=$HOME/workspace/rnaseq/student_tools/FastQC/fastqc

# Create output dir if doesnt exist
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/single_qc
mkdir -p $OUT_DIR/multi_qc

###
## Basic summary of an alignment
#

#$samtools flagstat $IN_DIR/Ad-Gfb.bam > $OUT_DIR/Ad-Gfb--flagstat.logs
#$samtools flagstat $IN_DIR/Ad-Rbm10.bam > $OUT_DIR/Ad-Rbm10--flagstat.logs


###
## Fastqc and MultiQC on Bam files
#

$fastqc --outdir $OUT_DIR/single_qc $IN_DIR/Ad-*_*.bam

multiqc --outdir $OUT_DIR/multi_qc $OUT_DIR/single_qc/.

