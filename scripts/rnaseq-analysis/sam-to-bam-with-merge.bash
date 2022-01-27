#!/bin/bash

# Dir Ref
main_folder="$HOME/workspace/rnaseq"
workspace="$HOME/workspace/seq-analysis-course"
IO_DIR=$workspace/alignment-data

# Tools Ref
samtools=$main_folder/student_tools/samtools-1.9/samtools
picard=$main_folder/student_tools/picard.jar


declare -a SamFilesName=("Ad-Rbm10_rep1" "Ad-Rbm10_rep2" "Ad-Rbm10_rep3" "Ad-Gfb_rep1" "Ad-Gfb_rep2" "Ad-Gfb_rep3")


# Convert HISAT2 sam files to bam files and sort by aligned position
for file in ${SamFilesName[@]};
do
   $samtools sort -@ 8 -o $IO_DIR/$file.bam $IO_DIR/$file.sam
done


# Index Bam Files -> .bai index files
for file in ${SamFilesName[@]};
do
   $samtools index $IO_DIR/$file.bam $IO_DIR/$file.bam.bai
done


#exit 1

# Merge HISAT2 BAM files
java -Xmx2g -jar $picard MergeSamFiles\
                 OUTPUT=$IO_DIR/Ad-Rbm10.bam\
                 INPUT=$IO_DIR/Ad-Rbm10_rep1.bam\
                 INPUT=$IO_DIR/Ad-Rbm10_rep2.bam\
                 INPUT=$IO_DIR/Ad-Rbm10_rep3.bam

java -Xmx2g -jar $picard MergeSamFiles\
                 OUTPUT=$IO_DIR/Ad-Gfb.bam\
                 INPUT=$IO_DIR/Ad-Gfb_rep1.bam\
                 INPUT=$IO_DIR/Ad-Gfb_rep2.bam\
                 INPUT=$IO_DIR/Ad-Gfb_rep3.bam

exit 1
