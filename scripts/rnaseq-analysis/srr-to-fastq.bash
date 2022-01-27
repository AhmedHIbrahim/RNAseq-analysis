#!/bin/bash

main_folder="$HOME/workspace/rnaseq"
workspace="$HOME/workspace/seq-analysis-course"

data_dir=$workspace/data

data_file=$data_dir/SRR_Acc_List.txt

output_dir=$workspace/fastq-data

fastq_dump=$HOME/workspace/rnaseq/sratoolkit.2.11.2-ubuntu64/bin/fastq-dump

while read p;
do
  $fastq_dump --split-3 --gzip --outdir $output_dir $data_dir/$p/$p.sra
done <"$data_file"


# ref: https://edwards.flinders.edu.au/fastq-dump/
