#!/bin/bash

fastqc=$HOME/workspace/rnaseq/student_tools/FastQC/fastqc
input_dir=$HOME/workspace/seq-analysis-course/alignment-data
output_dir=$HOME/workspace/seq-analysis-course/qc-results/alignment


mkdir -p $output_dir
mkdir -p $output_dir/single_qc


$fastqc --outdir $output_dir/single_qc $input_dir/*.bam

#multiqc --outdir $output_dir/multi_qc $output_dir/single_qc/.
