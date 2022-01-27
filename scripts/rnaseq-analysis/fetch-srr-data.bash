#!/bin/bash

main_folder="$HOME/workspace/rnaseq"
workspace="$HOME/workspace/seq-analysis-course"
data_file=$workspace/data/SRR_Acc_List.txt

prefetch=$HOME/workspace/rnaseq/sratoolkit.2.11.2-ubuntu64/bin/prefetch

while read p;
do
  $prefetch $p
done <"$data_file"
