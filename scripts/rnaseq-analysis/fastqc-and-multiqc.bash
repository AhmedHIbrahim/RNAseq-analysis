#!/bin/bash

fastqc=$HOME/workspace/rnaseq/student_tools/FastQC/fastqc
input_dir=$HOME/workspace/seq-analysis-course/fastq-data
output_dir=$HOME/workspace/seq-analysis-course/qc-results

$fastqc --outdir $output_dir/single_qc $input_dir/*.fastq.gz
# any Memory problem refer to  http://asearchforsolutions.blogspot.com/2019/04/fastqc-javalangoutofmemoryerror-java.html


multiqc --outdir $output_dir/multi_qc $output_dir/single_qc/.
