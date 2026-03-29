# Advanced usage

## Direct module execution

The basic usage commands in the above section are implemented as pipelines that execute a series of modules. To provide users with additional control and flexibility over this homology search process, users can directly call these modules.

| Pipeline                   | Modules            |
|----------------------------|-----------------------|
| `filter_existing_splits`   | `blastn_module`, `process_blast_results_module`, `filter_candidates_module`, `filter_test_split_module` |
| `stratify_test_split`      | `blastn_module`, `process_blast_results_module`, `stratify_test_split_module` |
| `create_orthogonal_splits` | `blastn_module`, `process_blast_results_module`, `filter_candidates_module`, `identify_homologous_groups_module`, `create_orthogonal_splits_module` |

Each module has a `-h` or `--help` command to print its respective arguments (e.g., `hashFrag blastn_module -h`).
 
One of the main advantages of calling modules individually is that it enables the use of manually computed pairwise scores. For example, after identifying candidate pairs of sequences with the BLAST algorithm, instead of using BLAST-derived alignment scores users can provide the exact [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) local alignment scores for candidate pairs. This was found to improve the recall of homologous sequences on the datasets tested at the expense of additional computational costs.

When providing precomputed pairwise scores to hashFrag, the expected format is a tab-delimited file with 3 columns: `id_i`, `id_j`, and `score`. 

`example.tsv`
```
seq_A	seq_B	60
seq_C	seq_D	85
seq_E	seq_A	100
...     ...     ...
```

For a full breakdown of available modules, please see the notebooks provided in the `/tutorial` directory.

## Creating orthogonal data folds

Users can also generate an arbitrary number of homology-aware data folds with the following module:

```
hashFrag create_orthogonal_folds_module -i $HOMOLOGY_PATH -f 10 -o $OUTPUT_DIR
```

where `-i` points to the file path containing the identified homologous groups (output of the `identify_homologous_groups_module`); `-f` corresponds to the number of folds to create; and `-o` defines the output directory to write the created folds.

This module implements a heuristic approach that aims to create folds of roughly equal size based on the defined groups. In particular, homologous groups are sorted according to the number of sequences before being greedily added to the smallest fold, where fold sizes are updated following the assignment of each homologous group. 

Note that performance will be unstable (i.e., variable fold sizes) if there are homologous groups substantially larger than the target fold size, which is determined by dividing the total number of sequences by the requested number of folds.

## hashFrag: High-Performance Computing (HPC) mode

Generating array job scripts for large-scale hashFrag runs on an HPC

For large-scale datasets, it is recommended to run hashFrag on an HPC cluster. The `blastn_array_module` can be used to generate an array job script that users can manually submit to carry out the BLAST similarity search process in partitions. Specifically, a BLAST database will be created over the complete input FASTA file, and the query FASTA file will be partitioned into smaller chunks based on the number of sequences, which can be set with the `--query-partition-size` argument. Each chunk of the query sequences is then setup to be executed independently as part of an array job.

> hashFrag currently supports the setup of array job scripts for the following job schedulers: SLURM, Sun Grid Engine (SGE)

Note that the `blastn_array_module` module expects the path to a shell script to carry out the setup of the HPC environment in each job (`--environment-path` argument). This script will be sourced at the start of each job. It should activate any virtual environments (if using), load relevant modules (e.g., `BLAST+`), and set the PATH variable.

After querying a partitioned FASTA file to the BLAST database constructed over the complete FASTA input, BLAST results will automatically processed to select the top alignment for each query-subject sequence pair with the `process_blast_results_module`.

### Example: SLURM

```
hashFrag blastn_array_module \
--train-fasta-path example_train_split.fa \
--test-fasta-path example_test_split.fa \
--query-partition-size 1000 \
--skip-revcomp \
--job-scheduler "slurm" \
--job-account "<account-name>" \
--job-time "6:00:00" \
--job-memory "16GB" \
--num-cpus 4 \
--environment-path "/<path_to>/set_env.sh" \
-o blastn_array_job.work
```

This will create the following job script file: `./blastn_array_job.work/hashFrag.blastn_array_module.array_jobs.sh`
```
#!/bin/bash
#SBATCH --account=<account-name>
#SBATCH --job-name=hashFrag
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --array=1-4
#SBATCH --output=blastn_array_job.work/job_stdout_files/hashFrag.jobid_%A_%a.out
#SBATCH --error=blastn_array_job.work/job_stderr_files/hashFrag.jobid_%A_%a.err

source /<path_to>/set_env.sh

BLASTDB_PATH=blastn_array_job.work/example_train_split.blastdb
BLAST_DIR=blastn_array_job.work/blast_results

# Each line is the path to the partitioned query FASTA file
QUERIES_PATH=blastn_array_job.work/query_fasta_paths.txt
QUERY_PATH=$( sed -n ${SLURM_ARRAY_TASK_ID}p $QUERIES_PATH )
LABEL=$( basename -s '.fa' $QUERY_PATH )

# Output files
BLASTN_PATH=$BLAST_DIR/${LABEL}.blastn.out
PROCESSED_BLASTN_PATH=$BLAST_DIR/${LABEL}.blastn.processed.tsv

blastn -query $QUERY_PATH -db $BLASTDB_PATH -out $BLASTN_PATH -word_size 11 -gapopen 2 -gapextend 1 -penalty -1 -reward 1 -max_target_seqs 500 -xdrop_ungap 20 -xdrop_gap 30 -xdrop_gap_final 100 -evalue 10 -dust no -num_threads 4 -outfmt '6 std score positive gaps' -strand plus

hashFrag process_blast_results_module --blastn-path $BLASTN_PATH --processed-blastn-path $PROCESSED_BLASTN_PATH
```
Submit this with the following command: `sbatch ./blastn_array_job.work/hashFrag.blastn_array_module.array_jobs.sh`

### Example: SGE

```
hashFrag blastn_array_module \
--train-fasta-path example_train_split.fa \
--test-fasta-path example_test_split.fa \
--query-partition-size 1000 \
--skip-revcomp \
--job-scheduler "sge" \
--job-memory "4G" \
--num-cpus 4 \
--environment-path "/<path_to>/set_env.sh" \
-o blastn_array_job.work
```