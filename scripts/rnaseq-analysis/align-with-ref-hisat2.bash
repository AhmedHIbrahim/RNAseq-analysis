#!/bin/bash


# Tool Ref
hisat2=$HOME/workspace/rnaseq/student_tools/hisat2-2.1.0/hisat2


# Dirs and files Ref
main_folder="$HOME/workspace/rnaseq"
workspace="$HOME/workspace/seq-analysis-course"

SRR_LIST_FILE=$data_folder/SRR_Acc_List.txt

DATA_DIR=$workspace/trimmed-data
OUTPUT_DIR=$workspace/alignment-data

RNA_REF_IDX=$workspace/ref/grch38_snp_tran/genome_snp_tran

mkdir -p $OUTPUT_DIR

declare -A types

types["SRR16798080"]="Ad-Rbm10 rep1 Mix1"
types["SRR16798081"]="Ad-Rbm10 rep2 Mix1"
types["SRR16798082"]="Ad-Rbm10 rep3 Mix1"

types["SRR16798083"]="Ad-Gfb rep1 Mix2"
types["SRR16798084"]="Ad-Gfb rep2 Mix2"
types["SRR16798085"]="Ad-Gfb rep3 Mix2"


for key in "${!types[@]}";
do
  val=(`echo "${types[$key]}"`)

  echo $key
  echo ${val[0]}
  echo ${val[1]}

  $hisat2 -p 8\
        --rg-id="${val[0]}_${val[1]}"\
        --rg SM:"${val[0]}"\
        --rg LB:"${val[0]}_${val[1]}_${val[2]}"\
        --rg PL:ILLUMINA\
        -x $RNA_REF_IDX --dta\
        --rna-strandness RF\
        -1 $DATA_DIR/$key"_1.fastq.gz"\
        -2 $DATA_DIR/$key"_2.fastq.gz"\
        -S $OUTPUT_DIR/"${val[0]}_${val[1]}".sam

done



exit 1


# Ad-Rbm10 => rep1 | rep2 | rep3
# Ad-Gfp   => rep1 | rep2 | rep3

#$hisat2 -p 8\
#        --rg-id=Ad-Rbm10_Rep1\
#        --rg SM:Ad-Rbm\
#        --rg LB:Ad-Rbm10_Rep1_ERCC-Mix1\
#        --rg PL:ILLUMINA\
#        --rg PU:CXX1234-ACTGAC.1\  <====
#        -x $RNA_REF_IDX --dta\
#        --rna-strandness RF\
#        -1 $DATA_DIR/SRR16798080_1.fastq.gz\
#        -2 $DATA_DIR/SRR16798080_2.fastq.gz\
#        -S ./Ad-Rbm10_Rep1.sam

#while read p;
#do
#   $hisat2 -p 8\
#           --rg-id=UHR_Rep1\
#           --rg SM:UHR\
#           --rg LB:UHR_Rep1_ERCC-Mix1\
#           --rg PL:ILLUMINA\
#           --rg PU:CXX1234-ACTGAC.1\
#           -x $RNA_REF_INDEX --dta\
#           --rna-strandness RF\
#           -1 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz\
#           -2 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz\
#           -S ./UHR_Rep1.sam
#done <"$SRR_LIST_FILE"

