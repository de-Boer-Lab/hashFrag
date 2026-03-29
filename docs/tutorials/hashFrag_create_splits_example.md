# hashFrag tutorial: Creating orthogonal splits

> This notebook refers to the case when users have a nucleotide sequence dataset and are interested in creating homology-aware train-test data splits for sequence-to-expression models.

This example workflow is performed on a subsampled MPRA dataset (K562) containing 10,000 sequences (provided in the `data` directory). When calling the `create_orthogonal_splits` pipeline, heuristic alignment scores derived from the `blastn` output are used to define similarity between sequences.

Example call of the complete pipeline (`lightning` mode):
```
hashFrag create_orthogonal_splits \
--fasta-path ../data/example_full_dataset.fa \
--word-size 7 \
--max-target-seqs 10000 \
--evalue 100 \
--threshold 60 \
--n-splits 10 \
--force \
--skip-revcomp \
--output-dir ../data/tutorial.create_orthogonal_splits.work
```

However, it may be desirable to instead use exact alignment scores (e.g., Smith-Waterman local alignment scores) for the homology search process. This notebook serves as a walkthrough for how users can use manually computed local alignment scores for the BLAST candidate pairs by calling the individual modules comprising the `create_orthogonal_splits` pipeline.

##  A note on the selected parameters for this tutorial 

Successful identification of cases of homology is paramount to effectively mitigate homology-based data leakage. As such, we configure the BLASTn parameters such that recall is maximized, even if it comes at the expense of increased false-positives. Here we consider the following parameters of BLASTn:

* `word_size`: smaller word sizes results in more exact word matches found between the query and sequences in the database, leading to more alignment score calculations being initialized.
* `max_target_seqs`: set to the size of the database to remove any constraints and allow for all possible candidate sequences to be returned for a given query.
* `evalue`: the e-value statistic is a measure of how likely you observe the alignment by chance (lower value corresponds to less likely to observe). By increasing the e-value threshold, less stringent matches that could be due to chance are returned.
* `dust`: by setting dust off, low-complexity (e.g., repetitive sequences) are no longer masked/filtered out.

An alignment score threshold of 60 was determined to be appropriate based on an analysis looking at alignment scores between dinucleotide shuffled (i.e., random) sequences.

# Section 1 - Identifying candidate similar sequences

The process of identifying candidate pairs of similar sequences involves first creating a BLAST database of the dataset, and then querying each sequence against the database. The BLASTn algorithm returns pairwise matches that represent potential cases of homology. 

Run the following command in terminal (e.g., Bash script):
```bash
FASTA_PATH=../data/example_full_dataset.fa
WORK_DIR=../data/tutorial.create_orthogonal_splits.work

hashFrag blastn_module \
--fasta-path $FASTA_PATH \
--max-target-seqs 10000 \
--word-size 7 \
--evalue 100 \
--blastdb-label "hashFrag" \
--skip-revcomp \
--output-dir $WORK_DIR
```
Output:

    2025-07-08 14:22:22 - blastn_module - INFO - Calling module...
    2025-07-08 14:22:22 - blastn_module - INFO - One FASTA files detected. Computing pairwise BLAST comparisons for all sequence-pairs...
    2025-07-08 14:22:34 - blastn_module - INFO - BLASTn output: 
    
    Building a new DB, current time: 07/08/2025 14:22:33
    New DB name:   /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/hashFrag.blastdb
    New DB title:  hashFrag
    Sequence type: Nucleotide
    Keep MBits: T
    Maximum file size: 1000000000B
    Adding sequences from FASTA; added 10000 sequences in 0.389443 seconds.
    
    2025-07-08 14:22:34 - blastn_module - INFO - BLAST DataBase construction finished and written to: /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/hashFrag.blastdb
    2025-07-08 14:28:27 - blastn_module - INFO - BLASTn process finished and written to: /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/hashFrag.blastn.out
    2025-07-08 14:28:27 - blastn_module - INFO - Module execution completed.
    

## Section 1.1 - Processing raw `blastn` output file

This processing step extracts the top-scoring alignment for each unique query-subject sequence pair and corrects the heuristic alignment score for subsequent steps. The processed tab-delimited file contains 3 columns (query sequence ID, subject sequence ID and their corrected heuristic alignment score).

Run the following command in terminal (e.g., Bash script):
```bash
WORK_DIR=../data/tutorial.create_orthogonal_splits.work
LABEL=hashFrag
BLASTN_PATH=$WORK_DIR/${LABEL}.blastn.out
PROCESSED_BLASTN_PATH=$WORK_DIR/${LABEL}.blastn.processed.tsv

hashFrag process_blast_results_module --blastn-path $BLASTN_PATH --processed-blastn-path $PROCESSED_BLASTN_PATH
```
Output:

    2025-07-08 14:34:56 - process_blast_results_module - INFO - Calling module...
    2025-07-08 14:34:58 - process_blast_results_module - INFO - Processed BLASTn results written to: ../data/tutorial.create_orthogonal_splits.work/hashFrag.blastn.processed.tsv
    2025-07-08 14:34:58 - process_blast_results_module - INFO - Module execution completed.
    


# Section 2: Filter false-positives based on a defined threshold

The next step involves filtering candidate pairings with alignment scores lower than the specified threshold. There are two different modes of hashFrag depending on what alignment score is selected.

1. `hashFrag-lightning` is the faster (and default) version where the alignment score computed from the BLAST output file. BLASTn is a heuristic method and the alignment scores were found to highly correlate with the optimal alignment scores; however, its underestimation of homology in some cases can lead to slightly worse recall. 

The following call performs the default behavior:
```
WORK_DIR=../data/tutorial.create_orthogonal_splits.work
INPUT_PATH=$WORK_DIR/hashFrag.blastn.processed.tsv
hashFrag filter_candidates_module -i $INPUT_PATH -t 60 -o $WORK_DIR
```

2. `hashFrag-pure` is the slower but more comprehensive method that is based on the optimal, Smith-Waterman local alignment scores between pairs of sequences. The calculation of optimal alignment scores incurs an additional cost to filtering.


## Section 2.1: hashFrag-pure mode

To limit memory usage, we'll start by partitioning the blast output file based on size. 

After completion of this step, all downstream steps will now be based on the homology identified using the exact alignment scores.

Run the following command in terminal (e.g., Bash script):
```bash
WORK_DIR=../data/tutorial.create_orthogonal_splits.work

cd $WORK_DIR
PROCESSED_BLASTN_PATH=$PWD/hashFrag.blastn.processed.tsv
BLAST_DIR=$PWD/blast_partitions
LABEL=$( basename -s ".tsv" $PROCESSED_BLASTN_PATH )

# Create directory for partitioned processed BLAST file
mkdir -p $BLAST_DIR
cd $BLAST_DIR

# Split the file based on number of lines
split -l 100000 -a 4 --additional-suffix=.tsv $PROCESSED_BLASTN_PATH ${LABEL}.partition_
ls -thor $BLAST_DIR
```
Output:

    total 4.0K
    -rw-r----- 1 brett 1.2M Jul  8 14:37 hashFrag.blastn.processed.partition_aaah.tsv
    -rw-r----- 1 brett 3.2M Jul  8 14:37 hashFrag.blastn.processed.partition_aaag.tsv
    -rw-r----- 1 brett 3.2M Jul  8 14:37 hashFrag.blastn.processed.partition_aaaf.tsv
    -rw-r----- 1 brett 3.2M Jul  8 14:37 hashFrag.blastn.processed.partition_aaae.tsv
    -rw-r----- 1 brett 3.2M Jul  8 14:37 hashFrag.blastn.processed.partition_aaad.tsv
    -rw-r----- 1 brett 3.2M Jul  8 14:37 hashFrag.blastn.processed.partition_aaac.tsv
    -rw-r----- 1 brett 3.2M Jul  8 14:37 hashFrag.blastn.processed.partition_aaab.tsv
    -rw-r----- 1 brett 3.2M Jul  8 14:37 hashFrag.blastn.processed.partition_aaaa.tsv


Run the following command in terminal (e.g., Bash script):
```bash
DATA_DIR=../data
cd $DATA_DIR

FASTA_PATH=$PWD/example_full_dataset.fa
WORK_DIR=$PWD/tutorial.create_orthogonal_splits.work
BLAST_DIR=$WORK_DIR/blast_partitions

cd ../src/external

echo "Computing exact alignment scores for partitioned files..."
for PARTITIONED_BLAST_PATH in $BLAST_DIR/*.partition_*.tsv
do
    echo $PARTITIONED_BLAST_PATH
    bash compute_blast_candidate_SW_scores.sh $FASTA_PATH $PARTITIONED_BLAST_PATH
done

echo "Concatenating partitioned files..."
cat $BLAST_DIR/*.pairwise_scores.tsv > $WORK_DIR/hashFrag_pure.blastn_candidates.sw_scores.tsv
cat $WORK_DIR/hashFrag_pure.blastn_candidates.sw_scores.tsv | head -n 10
```
Output:

    Computing exact alignment scores for partitioned files...
    /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/blast_partitions/hashFrag.blastn.processed.partition_aaaa.tsv
    /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/blast_partitions/hashFrag.blastn.processed.partition_aaab.tsv
    /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/blast_partitions/hashFrag.blastn.processed.partition_aaac.tsv
    /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/blast_partitions/hashFrag.blastn.processed.partition_aaad.tsv
    /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/blast_partitions/hashFrag.blastn.processed.partition_aaae.tsv
    /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/blast_partitions/hashFrag.blastn.processed.partition_aaaf.tsv
    /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/blast_partitions/hashFrag.blastn.processed.partition_aaag.tsv
    /home/brett/hashFrag/data/tutorial.create_orthogonal_splits.work/blast_partitions/hashFrag.blastn.processed.partition_aaah.tsv
    Concatenating partitioned files...
    BCL11A_1532_Reversed	peak24767_Reversed	14.0
    BCL11A_1532_Reversed	BCL11A_976	17.0
    BCL11A_1532_Reversed	ENSG00000139330_Reversed	18.0
    BCL11A_1532_Reversed	peak5045_Reversed	15.0
    BCL11A_1532_Reversed	peak70281_Reversed	18.0
    BCL11A_1532_Reversed	peak70514_Reversed	17.0
    BCL11A_1532_Reversed	peak67427_Reversed	15.0
    BCL11A_1532_Reversed	peak3273_Reversed	14.0
    BCL11A_1532_Reversed	HBA2_4771_Reversed	14.0
    BCL11A_1532_Reversed	peak13703	15.0

# Section 3: Determine groups of homology

There are often distinct groups of sequences exhibiting different cases of homology throughout the dataset. To determine such groups, we represent the "hits" (i.e., pairs of sequences with an alignment score greater than the threshold) as a sparse adjacency matrix. A graph can then be constructed, where nodes correspond to sequences and edges denote shared homology between the two sequences. The process of identifying groups of homology can readily be solved by identifying disconnected subgraphs. 

An efficient implementation for this graph-based task is provided in the `igraph` Python library.

Run the following command in terminal (e.g., Bash script):
```bash
WORK_DIR=../data/tutorial.create_orthogonal_splits.work
hashFrag identify_homologous_groups_module -i $WORK_DIR/hashFrag_pure.blastn_candidates.sw_scores.tsv -t 60 -o $WORK_DIR/homologous_groups.pure.csv
```
Output:

    2025-07-08 14:48:37 - identify_homologous_groups_module - INFO - Calling module...
    2025-07-08 14:48:40 - identify_homologous_groups_module - INFO - 4523 distinct groups.
    2025-07-08 14:48:40 - identify_homologous_groups_module - INFO - Homologous groups written to: ../data/tutorial.create_orthogonal_splits.work/homologous_groups.pure.csv
    2025-07-08 14:48:40 - identify_homologous_groups_module - INFO - Module execution completed.
    


# Section 4: Use case(s)

Upon identifying groups of sequences exhibiting high similarity (i.e., homology), we can create train-test data splits using a graph-based method. Specifically, by representing sequences as nodes and using edges to denote whether sequences were found to be homologous (yes or no), identifying homologous groups of sequences can be reduced to the task of identifying all disconnected subgraphs in the population. 

## Creating homology-aware data splits

Below we show how splits can be created based on the homologous groups identified from either the `hashFrag-lightning` or `hashFrag-pure` methods.

Run the following command in terminal (e.g., Bash script):
```bash
WORK_DIR=../data/tutorial.create_orthogonal_splits.work
HOMOLOGY_PATH=$WORK_DIR/homologous_groups.pure.csv # pure mode
OUT_DIR=$WORK_DIR
hashFrag create_orthogonal_splits_module -i $HOMOLOGY_PATH -n 10 -o $OUT_DIR
```
Output:

    2025-07-08 14:49:27 - create_orthogonal_splits_module - INFO - Calling module...
    2025-07-08 14:49:27 - create_orthogonal_splits_module - INFO - Creating 10 orthogonal splits in directory: ../data/tutorial.create_orthogonal_splits.work
    2025-07-08 14:49:27 - create_orthogonal_splits_module - INFO - Module execution completed.
    

# Creating homology-aware data folds

Run the following command in terminal (e.g., Bash script):
```bash
hashFrag create_orthogonal_folds_module -h
```
Output:

    usage: hashFrag create_orthogonal_folds_module [-h] -i HOMOLOGY_PATH
                                                   [-f FOLDS] [-s SEED] -o
                                                   OUTPUT_DIR
    
    Given the clustering of sequences based on homology
    ('identify_homologous_groups'), create homology-aware folds. Homologous groups
    are defined as disjoint sets over the population of sequences.
    
    optional arguments:
      -h, --help            show this help message and exit
      -i HOMOLOGY_PATH, --homology-path HOMOLOGY_PATH
                            The tab-delimited file containing the homologous group
                            labels that the relevant sequences in the test split
                            belong to (output file of
                            'identify_homologous_groups').
      -f FOLDS, --folds FOLDS
                            Number of folds to create.
      -s SEED, --seed SEED  Random seed to use for the creation of homology-aware
                            data splits (Default 21).
      -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                            The directory to write the created train-test splits.


Run the following command in terminal (e.g., Bash script):
```bash
WORK_DIR=../data/tutorial.create_orthogonal_splits.work
HOMOLOGY_PATH=$WORK_DIR/homologous_groups.pure.csv # pure mode
OUT_DIR=$WORK_DIR
hashFrag create_orthogonal_folds_module -i $HOMOLOGY_PATH -f 10 -o $OUT_DIR
```
Output:

    2025-07-08 14:53:48 - create_orthogonal_folds_module - INFO - Calling module...
    2025-07-08 14:53:48 - create_orthogonal_folds_module - WARNING - There exist(s) homologous groups of larger size than the expected fold size. Resulting folds may be imbalanced!
    2025-07-08 14:53:48 - create_orthogonal_folds_module - INFO - Creating 10 orthogonal folds...
    2025-07-08 14:53:48 - create_orthogonal_folds_module - INFO - Orthogonal folds written to ../data/tutorial.create_orthogonal_splits.work/hashFrag.10_orthogonal_folds.tsv
    2025-07-08 14:53:48 - create_orthogonal_folds_module - INFO - Module execution completed.

# Further details

Call the `help` command to list out all parameters.

Run the following command in terminal (e.g., Bash script):
```bash
hashFrag create_orthogonal_splits -h
```
Output:

    usage: hashFrag create_orthogonal_splits [-h] [-f FASTA_PATH] [-w WORD_SIZE]
                                             [-g GAPOPEN] [-x GAPEXTEND]
                                             [-p PENALTY] [-r REWARD]
                                             [-m MAX_TARGET_SEQS]
                                             [--exec-makeblastdb-only]
                                             [--skip-revcomp]
                                             [--xdrop-ungap XDROP_UNGAP]
                                             [--xdrop-gap XDROP_GAP]
                                             [--xdrop-gap_final XDROP_GAP_FINAL]
                                             [-e EVALUE] [-d DUST]
                                             [-b BLASTDB_ARGS]
                                             [--blastdb-label BLASTDB_LABEL]
                                             [-B BLASTN_ARGS] [-T THREADS] -t
                                             THRESHOLD [--p-train P_TRAIN]
                                             [--p-test P_TEST] [-n N_SPLITS]
                                             [-s SEED] [--force] [-o OUTPUT_DIR]
    
    Execute the full workflow of commands to create homology-aware train-test
    splits. This involves identifying identifying pairs of sequences sharing
    similarities with BLAST, filtering candidates based on a specified threshold,
    identifying all the different subgroups of sequences exhibiting a distinct
    case of homology, and creating train-test splits with no leakage.
    
    optional arguments:
      -h, --help            show this help message and exit
      -f FASTA_PATH, --fasta-path FASTA_PATH
                            Input FASTA file containing all sequences in the
                            dataset. All sequences will comprise the BLAST
                            database and each sequence will subsequently be
                            queried against it (supports unzipped or gzipped file
                            formats).
      -w WORD_SIZE, --word_size WORD_SIZE
                            Length of exact matching subsequences of initial match
                            (Default: 11).
      -g GAPOPEN, --gapopen GAPOPEN
                            Penalty (positive value) for opening gap in the
                            alignment (Default: 2).
      -x GAPEXTEND, --gapextend GAPEXTEND
                            Penalty (positive value) for extending an existing gap
                            in the alignment (Default: 1).
      -p PENALTY, --penalty PENALTY
                            Nucleotide mismatch in penalty (negative value) the
                            alignment (Default: -1).
      -r REWARD, --reward REWARD
                            Nucleotide match reward in the alignment (Default: 1).
      -m MAX_TARGET_SEQS, --max-target-seqs MAX_TARGET_SEQS
                            The maximum number of target sequences that can be
                            returned per query sequence (Default: 500).
      --exec-makeblastdb-only
                            Only run the makeblastdb command (default: False, set
                            to True when specified).
      --skip-revcomp        Skip generating reverse complement of sequences
                            comprising the BLAST database (Default: False,
                            generated if not skipped).
      --xdrop-ungap XDROP_UNGAP
                            X-drop threshold for ungapped alignment extension
                            (Permissible values: real numbers; Default: 20).
      --xdrop-gap XDROP_GAP
                            X-drop threshold for gapped alignment extension
                            (Permissible values: real numbers; Default: 30).
      --xdrop-gap_final XDROP_GAP_FINAL
                            X-drop threshold for final alignment extension
                            (Permissible values: real numbers; Default: 100).
      -e EVALUE, --evalue EVALUE
                            The likelihood threshold required to report sequences
                            as a match (Default: 10).
      -d DUST, --dust DUST  Filter low-complexity (e.g., repetitive) regions
                            (Default: 'no').
      -b BLASTDB_ARGS, --blastdb-args BLASTDB_ARGS
                            Pass additional arguments for makeblastdb call.
      --blastdb-label BLASTDB_LABEL
                            A label for the BLAST database.
      -B BLASTN_ARGS, --blastn-args BLASTN_ARGS
                            Pass additional arguments for blastn call.
      -T THREADS, --threads THREADS
                            The number of CPUs for database search (Default: 1).
      -t THRESHOLD, --threshold THRESHOLD
                            Alignment score threshold to discern a pair of
                            sequences as homologous or a false-positive candidate.
      --p-train P_TRAIN     The proportion of sequences to send to the train data
                            split (Default: 0.8).
      --p-test P_TEST       The proportion of sequences to send to the test data
                            split (Default: 0.2).
      -n N_SPLITS, --n-splits N_SPLITS
                            Number of split replicates to create (Default: 1).
      -s SEED, --seed SEED  Random seed to use for the creation of homology-aware
                            data splits (Default 21).
      --force               Force overwrite existing BLAST module output files
                            (Default: False, existing output files will not be
                            overwritten).
      -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                            The directory to write the created train-test splits.
