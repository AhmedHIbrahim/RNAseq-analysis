#!/bin/bash


WORK_DIR=$HOME/workspace/seq-analysis-course/de
cd $WORK_DIR

cut -f 1 $WORK_DIR/ballgown/ref_only/DE_genes.txt | sort  > ballgown_DE_gene_symbols.txt
cut -f 2 $WORK_DIR/htseq_counts/DE_genes.txt | sort > htseq_counts_edgeR_DE_gene_symbols.txt
