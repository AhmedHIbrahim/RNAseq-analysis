#!/bin/bash


# 1) Create a ref file to your data stringtie output

WORKSPACE="$HOME/workspace/seq-analysis-course"
REF_DIR="$WORKSPACE/expression/stringtie/ref_only"
BALLGOWN_DIR="$WORKSPACE/de/ballgown/ref_only/"

mkdir -p $BALLGOWN_DIR
cd $BALLGOWN_DIR

printf "\"ids\",\"type\",\"path\"\n\"Ad-Rbm10_rep1\",\"Ad-Rbm10\",\"$REF_DIR/Ad-Rbm10_rep1\"\n\"Ad-Rbm10_rep2\",\"Ad-Rbm10\",\"$REF_DIR/Ad-Rbm10_rep2\"\n\"Ad-Rbm10_rep3\",\"Ad-Rbm10\",\"$REF_DIR/Ad-Rbm10_rep3\"\n\"Ad-Gfb_rep1\",\"Ad-Gfb\",\"$REF_DIR/Ad-Gfb_rep1\"\n\"Ad-Gfb_rep2\",\"Ad-Gfb\",\"$REF_DIR/Ad-Gfb_rep2\"\n\"Ad-Gfb_rep3\",\"Ad-Gfb\",\"$REF_DIR/Ad-Gfb_rep3\"\n" > Ad-Rbm10_vs_Ad-Gfb.csv

cat Ad-Rbm10_vs_Ad-Gfb.csv
