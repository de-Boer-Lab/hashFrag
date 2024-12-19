#!/bin/bash

FASTA_PATH=$1
BLAST_PATH=$2

OUTPUT_PATH=$( echo $BLAST_PATH | sed 's/\.tsv/\.pairwise_scores\.tsv\.gz/g')

python compute_pairwise_SW_scores.py \
-i $FASTA_PATH \
-b $BLAST_PATH \
-o $OUTPUT_PATH
