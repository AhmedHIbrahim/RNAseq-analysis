#!/bin/bash

GTF_REF=$HOME/workspace/seq-analysis-course/ref/Homo_sapiens.GRCh38.99.gtf
OUT_DIR=$HOME/workspace/seq-analysis-course/de/htseq_counts

mkdir -p $OUT_DIR
cd $OUT_DIR

perl -ne 'if ($_ =~ /gene_id\s\"(ENSG\S+)\"\;/) { $id = $1; $name = undef; if ($_ =~ /gene_name\s\"(\S+)"\;/) { $name = $1; }; }; if ($id && $name) {print "$id\t$name\n";} if ($_=~/gene_id\s\"(ERCC\S+)\"/){print "$1\t$1\n";}' $GTF_REF | sort | uniq > ENSG_ID2Name.txt
