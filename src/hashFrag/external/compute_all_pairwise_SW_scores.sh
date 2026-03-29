#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N bkSW
#$ -l s_vmem=4G
#$ -l arm
#$ -pe def_slot 16
#$ -j y

# qsub compute_all_pairwise_SW_scores.sh /home/brett/work/OrthogonalTrainValSplits/hashFrag/data/K562.sample_10000.fa.gz

CPU=16

source $HOME/setenv/miniconda_arm.sh

FASTA_PATH=$1
OUTPATH=$( echo $FASTA_PATH | sed 's/.fa.gz/.pairwise_scores.tsv.gz/g' )

python compute_pairwise_SW_scores.py \
-i $FASTA_PATH \
-n $CPU \
-o $OUTPATH
