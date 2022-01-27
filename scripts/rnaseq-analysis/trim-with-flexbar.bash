#!/bin/bash

main_folder="$HOME/workspace/rnaseq"
workspace="$HOME/workspace/seq-analysis-course"

data_folder=$workspace/data
fastq_data_folder=$workspace/fastq-data

data_file=$data_folder/SRR_Acc_List.txt
output_folder=$workspace/trimmed-data

flexbar=$HOME/workspace/rnaseq/student_tools/flexbar-3.4.0-linux/flexbar

while read p;
do
  $flexbar --adapter-min-overlap 7  \
           --adapter-trim-end RIGHT \
           --adapters $workspace/ref/illumina_multiplex.fa \
           --max-uncalled 300 \
           --pre-trim-left 13 \
           --min-read-length 25 \
           --threads 8 \
           --zip-output GZ \
           --reads $fastq_data_folder/$p"_1.fastq.gz" \
           --reads2 $fastq_data_folder/$p"_2.fastq.gz" \
           --target $output_folder/$p
done <"$data_file"


<<'###FULL-FLEXBAR-COMMENT'

$flexbar --adapter-min-overlap 7  \
           --adapter-trim-end RIGHT \
           --adapters $workspace/ref/illumina_multiplex.fa \
           --pre-trim-left 13 \
           --max-uncalled 300 \
           --min-read-length 25 \
           --threads 8 \
           --zip-output GZ \
           --reads $fastq_data_folder/$p"_1.fastq.gz" \
           --reads2 $fastq_data_folder/$p"_2.fastq.gz" \
           --target $output_folder/$p
###FULL-FLEXBAR-COMMENT
