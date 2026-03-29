# hashFrag tutorial: stratify existing splits

> This notebook refers to the case when users have existing train-test splits and are interested in identifying and mitigating data leakage attributed to shared sequence homology across splits.

The example workflow is performed with two data splits from a subsampled MPRA dataset (K562): a 8,000-sequence train split and a 2,000-sequence test split (provided in the `data` directory).

This notebook serves a walkthrough of calling the individual modules comprising the `hashFrag stratify_test_split` command.

Example call of the complete pipeline (`lightning` mode):

```
hashFrag stratify_test_split \
--train-fasta-path ../data/example_train_split.f \
--test-fasta-path ../data/example_test_split.fa \
--word-size 7 \
--max-target-seqs 8000 \
--evalue 100 \
--step 10 \
--skip-revcomp \
--force \
--output-dir ../data/tutorial.stratify_test_split.work
```

However, it may be desirable to instead use exact alignment scores (e.g., Smith-Waterman local alignment scores) for the homology search process. This notebook serves as a walkthrough for how users can use manually computed local alignment scores for the BLAST candidate pairs by calling the individual modules comprising the `create_orthogonal_splits` pipeline.

## A note on the selected parameters for this tutorial 

Successful identification of cases of homology is paramount to effectively mitigate homology-based data leakage. As such, we configure the BLASTn parameters such that recall is maximized, even if it comes at the expense of increased false-positives. Here we consider the following parameters of BLASTn:

* `word_size`: smaller word sizes results in more exact word matches found between the query and sequences in the database, leading to more alignment score calculations being initialized.
* `max_target_seqs`: set to the size of the database to remove any constraints and allow for all possible candidate sequences to be returned for a given query.
* `evalue`: the e-value statistic is a measure of how likely you observe the alignment by chance (lower value corresponds to less likely to observe). By increasing the e-value threshold, less stringent matches that could be due to chance are returned.
* `dust`: by setting dust off, low-complexity (e.g., repetitive sequences) are no longer masked/filtered out.

## Section 1 - Identifying candidate similar sequences

When user-derived train-test splits are provided, comparisons are constrained to pairs of sequences across splits. The process of identifying candidate pairs of similar sequences involves first creating a BLAST database of sequences in the train split, and then querying each test split sequence against the database. The BLASTn algorithm returns pairwise matches that represent potential cases of homology. 

Run the following command in terminal (e.g., Bash script):
```bash
TRAIN_FASTA_PATH=../data/example_train_split.fa
TEST_FASTA_PATH=../data/example_test_split.fa
WORK_DIR=../data/tutorial.stratify_test_split.work

hashFrag blastn_module \
--train-fasta-path $TRAIN_FASTA_PATH \
--test-fasta-path $TEST_FASTA_PATH \
--word-size 7 \
--max-target-seqs 8000 \
--evalue 100 \
--blastdb-label "hashFrag" \
--skip-revcomp \
-o $WORK_DIR
```
Output:

    2025-03-02 08:30:33 - blastn_module - INFO - Calling module...
    2025-03-02 08:30:33 - blastn_module - INFO - Train and test FASTA files detected. Computing pairwise BLAST comparisons across splits...
    2025-03-02 08:30:39 - blastn_module - INFO - BLASTn output: 
    
    Building a new DB, current time: 03/02/2025 08:30:38
    New DB name:   /home/brett/hashFrag/data/tutorial.stratify_test_split.work/hashFrag.blastdb
    New DB title:  hashFrag
    Sequence type: Nucleotide
    Keep MBits: T
    Maximum file size: 1000000000B
    Adding sequences from FASTA; added 8000 sequences in 0.12439 seconds.
    
    2025-03-02 08:30:39 - blastn_module - INFO - BLAST DataBase construction finished and written to: /home/brett/hashFrag/data/tutorial.stratify_test_split.work/hashFrag.blastdb
    2025-03-02 08:32:31 - blastn_module - INFO - BLASTn process finished and written to: /home/brett/hashFrag/data/tutorial.stratify_test_split.work/example_test_split.blastn.out
    2025-03-02 08:32:31 - blastn_module - INFO - Module execution completed.
    
### Section 1.1 - Processing raw `blastn` output file

This processing step extracts the top-scoring alignment for each unique query-subject sequence pair and corrects the heuristic alignment score for subsequent steps. The processed tab-delimited file contains 3 columns (query sequence ID, subject sequence ID and their corrected heuristic alignment score).

Run the following command in terminal (e.g., Bash script):
```bash
WORK_DIR=../data/tutorial.stratify_test_split.work
LABEL=example_test_split
BLASTN_PATH=$WORK_DIR/${LABEL}.blastn.out
PROCESSED_BLASTN_PATH=$WORK_DIR/${LABEL}.blastn.processed.tsv

hashFrag process_blast_results_module --blastn-path $BLASTN_PATH --processed-blastn-path $PROCESSED_BLASTN_PATH
```
Output:

    2025-03-02 08:33:12 - process_blast_results_module - INFO - Calling module...
    2025-03-02 08:33:13 - process_blast_results_module - INFO - Processed BLASTn results written to: ../data/tutorial.stratify_test_split.work/example_test_split.blastn.processed.tsv
    2025-03-02 08:33:13 - process_blast_results_module - INFO - Module execution completed.

## Section 2: Use Case(s)

### Stratify test split based on homology

Another potentially useful feature is to stratify the test split based on each test sequence's maximum alignment score compared to all sequences in the train split. This can aid in studying the effects that homology has on model performance evaluation. The range of values is specified by the `step` parameter.

The default behavior (`lightning` mode), can be obtained with the following command:

```
TEST_FASTA_PATH=../data/example_test_split.fa
WORK_DIR=../data/tutorial.stratify_test_split.work
INPUT_PATH=$WORK_DIR/hashFrag.blastn.out
OUTPUT_PATH=$WORK_DIR/hashFrag_lightning.stratified_test.tsv
hashFrag stratify_test_split_module -f $TEST_FASTA_PATH -i $INPUT_PATH -s 10 -o $OUTPUT_PATH
```

### Section 2.1 - hashFrag-pure mode

Run the following command in terminal (e.g., Bash script):
```bash
WORK_DIR=../data/tutorial.stratify_test_split.work

cd $WORK_DIR
PROCESSED_BLASTN_PATH=$PWD/example_test_split.blastn.processed.tsv
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

    total 1.5K
    -rw-r----- 1 brett 1.4M Mar  2 08:35 example_test_split.blastn.processed.partition_aaac.tsv
    -rw-r----- 1 brett 3.2M Mar  2 08:35 example_test_split.blastn.processed.partition_aaab.tsv
    -rw-r----- 1 brett 3.2M Mar  2 08:35 example_test_split.blastn.processed.partition_aaaa.tsv

Run the following command in terminal (e.g., Bash script):
```bash
DATA_DIR=../data
cd $DATA_DIR

FASTA_PATH=$PWD/example_full_dataset.fa
WORK_DIR=$PWD/tutorial.stratify_test_split.work
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
    /home/brett/hashFrag/data/tutorial.stratify_test_split.work/blast_partitions/example_test_split.blastn.processed.partition_aaaa.tsv
    /home/brett/hashFrag/data/tutorial.stratify_test_split.work/blast_partitions/example_test_split.blastn.processed.partition_aaab.tsv
    /home/brett/hashFrag/data/tutorial.stratify_test_split.work/blast_partitions/example_test_split.blastn.processed.partition_aaac.tsv
    Concatenating partitioned files...
    BCL11A_1238	peak79295	15.0
    BCL11A_1238	peak79031	13.0
    BCL11A_1238	peak16460_Reversed	15.0
    BCL11A_1238	peak10840_Reversed	13.0
    BCL11A_1238	GATA1_4792	14.0
    BCL11A_1238	peak35733_Reversed	18.0
    BCL11A_1238	peak74971	14.0
    BCL11A_1238	peak4007_Reversed	14.0
    BCL11A_1238	peak27092_Reversed	13.0
    BCL11A_1238	peak48369	13.0


The SW scores for candidate pairs of sequences are concatenated into a single `.tsv` file, which is provided as input to the `stratify_test_split_module` command.

Run the following command in terminal (e.g., Bash script):
```bash
TEST_FASTA_PATH=../data/example_test_split.fa
WORK_DIR=../data/tutorial.stratify_test_split.work
INPUT_PATH=$WORK_DIR/hashFrag_pure.blastn_candidates.sw_scores.tsv
STEP=10
OUTPUT_PATH=$WORK_DIR/hashFrag_pure.stratified_test.tsv

hashFrag stratify_test_split_module -f $TEST_FASTA_PATH -i $INPUT_PATH -s 10 -o $OUTPUT_PATH
```
Output:

    2025-03-02 08:40:46 - stratify_test_split_module - INFO - Calling module...
    2025-03-02 08:40:46 - stratify_test_split_module - INFO - Stratifying based on the provided pairwise scores...
    2025-03-02 08:40:46 - stratify_test_split_module - INFO - Stratification results written to: ../data/tutorial.stratify_test_split.work/hashFrag_pure.stratified_test.tsv
    2025-03-02 08:40:46 - stratify_test_split_module - INFO - Module execution completed.

## Further details

Run the following command in terminal (e.g., Bash script):
```bash
hashFrag stratify_test_split -h
```
Output:

    usage: hashFrag stratify_test_split [-h] [--train-fasta-path TRAIN_FASTA_PATH] [--test-fasta-path TEST_FASTA_PATH]
                                        [-w WORD_SIZE] [-g GAPOPEN] [-x GAPEXTEND] [-p PENALTY] [-r REWARD]
                                        [-m MAX_TARGET_SEQS] [--exec-makeblastdb-only] [--skip-revcomp]
                                        [--xdrop-ungap XDROP_UNGAP] [--xdrop-gap XDROP_GAP]
                                        [--xdrop-gap_final XDROP_GAP_FINAL] [-e EVALUE] [-d DUST] [-b BLASTDB_ARGS]
                                        [--blastdb-label BLASTDB_LABEL] [-B BLASTN_ARGS] [-T THREADS] [--force]
                                        [-o OUTPUT_DIR] [-s STEP]
    
    Execute the full workflow of commands to stratify the test split based on their maximum pairwise alignment score
    to the train split sequences. This involves identifying pairs of sequences sharing similarities with BLAST, and
    then stratiyfing the test split into subsplits based on the corrected BLAST alignment score.
    
    optional arguments:
      -h, --help            show this help message and exit
      --train-fasta-path TRAIN_FASTA_PATH
                            Input FASTA file for the training data split, which will comprise the BLAST database.
                            (supports unzipped or gzipped file formats)
      --test-fasta-path TEST_FASTA_PATH
                            Each sequence will be queried against the train split BLAST database. (supports unzipped
                            or gzipped file formats)
      -w WORD_SIZE, --word_size WORD_SIZE
                            Length of exact matching subsequences of initial match (Default: 11).
      -g GAPOPEN, --gapopen GAPOPEN
                            Penalty (positive value) for opening gap in the alignment (Default: 2).
      -x GAPEXTEND, --gapextend GAPEXTEND
                            Penalty (positive value) for extending an existing gap in the alignment (Default: 1).
      -p PENALTY, --penalty PENALTY
                            Nucleotide mismatch in penalty (negative value) the alignment (Default: -1).
      -r REWARD, --reward REWARD
                            Nucleotide match reward in the alignment (Default: 1).
      -m MAX_TARGET_SEQS, --max-target-seqs MAX_TARGET_SEQS
                            The maximum number of target sequences that can be returned per query sequence (Default:
                            500).
      --exec-makeblastdb-only
                            Only run the makeblastdb command (default: False, set to True when specified).
      --skip-revcomp        Skip generating reverse complement of sequences comprising the BLAST database (Default:
                            False, generated if not skipped).
      --xdrop-ungap XDROP_UNGAP
                            X-drop threshold for ungapped alignment extension (Permissible values: real numbers;
                            Default: 20).
      --xdrop-gap XDROP_GAP
                            X-drop threshold for gapped alignment extension (Permissible values: real numbers;
                            Default: 30).
      --xdrop-gap_final XDROP_GAP_FINAL
                            X-drop threshold for final alignment extension (Permissible values: real numbers; Default:
                            100).
      -e EVALUE, --evalue EVALUE
                            The likelihood threshold required to report sequences as a match (Default: 10).
      -d DUST, --dust DUST  Filter low-complexity (e.g., repetitive) regions (Default: 'no').
      -b BLASTDB_ARGS, --blastdb-args BLASTDB_ARGS
                            Pass additional arguments for makeblastdb call.
      --blastdb-label BLASTDB_LABEL
                            A label for the BLAST database.
      -B BLASTN_ARGS, --blastn-args BLASTN_ARGS
                            Pass additional arguments for blastn call.
      -T THREADS, --threads THREADS
                            The number of CPUs for database search (Default: 1).
      --force               Force overwrite existing BLAST module output files (Default: False, existing output files
                            will not be overwritten).
      -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                            Directory to write BLASTn results (Default: '.').
      -s STEP, --step STEP  The step size refers to how large each alignment score range is (Default: 10).


