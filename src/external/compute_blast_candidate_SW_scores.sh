#!/bin/bash

FASTA_PATH=$1
BLAST_PATH=$2

OUTPUT_PATH=$BLAST_PATH.pairwise_scores.csv.gz

python compute_pairwise_SW_scores.py \
-i $FASTA_PATH \
-b $BLAST_PATH \
-o $OUTPUT_PATH
