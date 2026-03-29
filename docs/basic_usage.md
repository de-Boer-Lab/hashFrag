# Basic usage

hashFrag supports 3 primary use cases.

## If you already have train/test data splits


Existing train-test data splits can be handled by providing two separate FASTA files to hashFrag as input. This limits the homology search process to inter-data split comparisons. Specifically, a BLAST database is constructed over the train sequences and the test sequences are queried against this database to identify pairs with high local alignment scores. 

1. Filter existing splits: remove test sequences exhibiting homology to the training set

2. Stratify the test split based on homology ot the training set.

## If you want to create new splits

When a single FASTA file is provided as input, hashFrag will characterize homology for all pairwise comparisons. This involves constructing a BLAST database over all sequences in the population and then subsequently querying each sequence to the databae. 

3. Create orthogonal, homology-aware train-test splits

## Example dataset

Example usage of hashFrag is provided below on an example dataset composed of 10,000 sequences (each 200 base pairs in length).

Note that the `filter_existing_splits` and `create_orthogonal_splits` pipelines require specification of a pairwise alignment score threshold to define homology between sequences.

> If users do not have a predefined threshold for homology, we recommend computing pairwise alignment scores between a set of random (e.g., dinucleotide shuffled) genomic sequences, and then defining a threshold above the distribution of values observed.

## Defining homology

Defining homology in terms of alignment scores requires the specification of scoring parameters (e.g., `penalty`, `gapopen`, `gapextend`, and `reward` values).

> To remain consistent with BLAST scoring parameter specification, hashFrag expects a negative value for `penalty` but a positive value for `gapopen` and `gapextend` arguments (gap penalties will be subtracted from the alignment score during calculation). A positive value is expected for `reward`.

Changing these parameters can drastically impact the identification process of homology. Please see permissible scoring parameter combinations for the BLASTn algorithm [here](https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/) (Table D1).

## hashFrag syntax

By default, hashFrag generates reverse complement sequences when creating the BLAST database. This ensures that query sequences are assessed for homology for both sequence orientations. If the input FASTA files already contain both forward and reverse sequence orientations, include the `--skip-revcomp` argument to skip this step.

> Currently, hashFrag expects reverse complementary sequences to be denoted with a `_Reversed` suffix in the sequence header.

If both sequence orentations already exist in the input FASTA file(s), make sure one of the orientations is denoted with the `_Reversed` suffix. For example, if a sequence has a header such as `seq_A`, ensure that its reverse complement has the header `seq_A_Reversed`.

* Note that a warning will be generated if there already exists sequences with the `_Reversed` suffix and the `--skip-revcomp` argument is NOT specified.

# hashFrag commands

## Filter existing splits

Filter sequences in the test split exhibiting homology with any sequences in the train split
```
hashFrag filter_existing_splits \
--train-fasta-path example_train_split.fa \
--test-fasta-path example_test_split.fa \
-t 60 \
--skip-revcomp \
-o filter_existing_splits.work
```

<details>
<summary>Expand to view the full table of arguments</summary>

| Argument | Description | Expected input |
|---|---|---|
| `--train-fasta-path` | Input file containing train split sequences. | FASTA file path (unzipped) |
| `--test-fasta-path` | Input file containing test split sequences. | FASTA file path (unzipped) |
| `-w`, `--word-size` | Length of exact match to intialize alignment score calculation (`blastn_module`). | integer (Default: 11) |
| `-g`, `--gapopen` | Penalty for opening a gap in the alignment (`blastn_module`) | positive integer (Default: 2) |
| `-x`, `--gapextend` | Penalty for extending an existing gap in the alignment (`blastn_module`). | positive integer (Default: 1) |
| `-p`, `--penalty` | Nucleotide mismatch penalty (`blastn_module`). | negative integer (Default: -1) |
| `-r`, `--reward` | Nucleotide match reward (`blastn_module`). | positive integer (Default: 1) |
| `-m`, `--max-target-seqs` | Maximum number of target sequences that can be returned per query sequence (`blastn_module`). | positive integer (Default: 500) |
| `--exec-makeblastdb-only` | Only run the makeblastdb command (`blastn_module`). | Boolean (Default: False, set to True when specified) |
| `--skip-revcomp` | Do not generate reverse complement of sequences comprising the BLAST database  (`blastn_module`). | Boolean (Default: False, generated if not skipped) |
| `--xdrop-ungap`| X-drop threshold (heuristic value in bits) for ungapped alignment extension (`blastn_module`). | real number (Default: 20) |
| `--xdrop-gap` | X-drop threshold (heuristic value in bits) for gapped alignment extension (`blastn_module`). | real number (Default: 30) |
| `--xdrop-gap-final` | X-drop threshold (heuristic value in bits) for final alignment extension (`blastn_module`) | real number (Default: 100) |
| `-e`, `--e-value` | Likelihood threshold required to report a sequence as a match (`blastn_module`). | real number (Default: 10.0) |
| `-d`, `--dust` | Filter for low-complexity (i.e., repetitive) regions (`blastn_module`). | Permissible values: {'yes', 'no'} (Default: 'no') |
| `--blastdb-label` | Label for the BLAST database (`blastn_module`). | string (Default: None) |
| `-T`, `--threads` | Number of threads to use for `blastn_module` execution. | positive integer (Default: 1) |
| `-t`, `--threshold` | Alignment score threshold to define a pair of sequences as similar, or homologous (`filter_candidates_module`). | all real numbers (*Required*) |
| `--force` | Force overwrite existing `blastn_module` output files. | Boolean (Default: False, set to True when specified) |
| `-o`, `--output-dir` | Directory to write intermediate results. | string (Default: '.') |

</details>

## Stratify test split

Stratify the test split sequences into an arbitrary number of levels based on their maximum alignment scores to the train split sequences. 
```
hashFrag stratify_test_split \
--train-fasta-path example_train_split.fa \
--test-fasta-path example_test_split.fa \
--skip-revcomp \
-o stratify_test_split.work
```

* Note that the sizes of each stratified level will not necessarily be balanced.
* This can be useful to better understand a model’s behaviour over test splits at varying levels of orthogonality to the sequences the model was trained on.


<details>
<summary>Expand to view the full table of arguments</summary>

| Argument | Description | Expected input |
|---|---|---|
| `--train-fasta-path` | Input file containing train split sequences. | FASTA file path (unzipped) |
| `--test-fasta-path` | Input file containing test split sequences. | FASTA file path (unzipped) |
| `-w`, `--word-size` | Length of exact match to intialize alignment score calculation (`blastn_module`). | integer (Default: 11) |
| `-g`, `--gapopen` | Penalty for opening a gap in the alignment (`blastn_module`) | positive integer (Default: 2) |
| `-x`, `--gapextend` | Penalty for extending an existing gap in the alignment (`blastn_module`). | positive integer (Default: 1) |
| `-p`, `--penalty` | Nucleotide mismatch penalty (`blastn_module`). | negative integer (Default: -1) |
| `-r`, `--reward` | Nucleotide match reward (`blastn_module`). | positive integer (Default: 1) |
| `-m`, `--max-target-seqs` | Maximum number of target sequences that can be returned per query sequence (`blastn_module`). | positive integer (Default: 500) |
| `--exec-makeblastdb-only` | Only run the makeblastdb command (`blastn_module`). | Boolean (Default: False, set to True when specified) |
| `--skip-revcomp` | Do not generate reverse complement of sequences comprising the BLAST database  (`blastn_module`). | Boolean (Default: False, generated if not skipped) |
| `--xdrop-ungap`| X-drop threshold (heuristic value in bits) for ungapped alignment extension (`blastn_module`). | real number (Default: 20) |
| `--xdrop-gap` | X-drop threshold (heuristic value in bits) for gapped alignment extension (`blastn_module`). | real number (Default: 30) |
| `--xdrop-gap-final` | X-drop threshold (heuristic value in bits) for final alignment extension (`blastn_module`) | real number (Default: 100) |
| `-e`, `--e-value` | Likelihood threshold required to report a sequence as a match (`blastn_module`). | real number (Default: 10.0) |
| `-d`, `--dust` | Filter for low-complexity (i.e., repetitive) regions (`blastn_module`). | Permissible values: {'yes', 'no'} (Default: 'no') |
| `--blastdb-label` | Label for the BLAST database (`blastn_module`). | string (Default: None) |
| `-T`, `--threads` | Number of threads to use for `blastn_module` execution. | positive integer (Default: 1) |
| `-s`, `--step` | Step size for how large each alignment score range is (`stratify_test_split_module`). | positive integer (Default: 10) |
| `--force` | Force overwrite existing `blastn_module` output files. | Boolean (Default: False) |
| `-o`, `--output-dir` | Directory to write the created train-test splits. | string (Default: '.') |

</details>

## Creating orthogonal splits

Create homology-aware (i.g., orthogonal) train-test data splits.

```
hashFrag create_orthogonal_splits \
-f example_full_dataset.fa \
-t 60 \
--skip-revcomp \
-o create_orthogonal_splits.work
```
The creation of orthogonal train-test splits involves determining disjoint sets of sequences with respect to homology. This is accomplished using the [union-find](https://en.wikipedia.org/wiki/Disjoint-set_data_structure) data structure (also referred to as disjoint-set or merge-find set). From the homology cluster information, splits with no leakage can be created proportionally. 

<details>
<summary>Expand to view the full table of arguments</summary>

| Argument | Description | Expected input |
|---|---|---|
| `-f`, `--fasta-path` | Input file containing all sequences in the dataset. | FASTA file path (unzipped) |
| `-w`, `--word-size` | Length of exact match to intialize alignment score calculation (`blastn_module`). | integer (Default: 11) |
| `-g`, `--gapopen` | Penalty for opening a gap in the alignment (`blastn_module`) | positive integer (Default: 2) |
| `-x`, `--gapextend` | Penalty for extending an existing gap in the alignment (`blastn_module`). | positive integer (Default: 1) |
| `-p`, `--penalty` | Nucleotide mismatch penalty (`blastn_module`). | negative integer (Default: -1) |
| `-r`, `--reward` | Nucleotide match reward (`blastn_module`). | positive integer (Default: 1) |
| `-m`, `--max-target-seqs` | Maximum number of target sequences that can be returned per query sequence (`blastn_module`). | positive integer (Default: 500) |
| `--exec-makeblastdb-only` | Only run the makeblastdb command (`blastn_module`). | Boolean (Default: False, set to True when specified) |
| `--skip-revcomp` | Do not generate reverse complement of sequences comprising the BLAST database  (`blastn_module`). | Boolean (Default: False, generated if not skipped) |
| `--xdrop-ungap`| X-drop threshold (heuristic value in bits) for ungapped alignment extension (`blastn_module`). | real number (Default: 20) |
| `--xdrop-gap` | X-drop threshold (heuristic value in bits) for gapped alignment extension (`blastn_module`). | real number (Default: 30) |
| `--xdrop-gap-final` | X-drop threshold (heuristic value in bits) for final alignment extension (`blastn_module`) | real number (Default: 100) |
| `-e`, `--e-value` | Likelihood threshold required to report a sequence as a match (`blastn_module`). | real number (Default: 10.0) |
| `-d`, `--dust` | Filter for low-complexity (i.e., repetitive) regions (`blastn_module`). | Permissible values: {'yes', 'no'} (Default: 'no') |
| `--blastdb-label` | Label for the BLAST database (`blastn_module`). | string (Default: None) |
| `-T`, `--threads` | Number of threads to use for `blastn_module` execution. | positive integer (Default: 1) |
| `-t`, `--threshold` | Alignment score threshold to define a pair of sequences as similar, or homologous (`filter_candidates_module`). | all real numbers (*Required*) |
| `--p-train` | Proportion of sequences for the newly-created train data split (`create_orthogonal_splits_module`). | float (Default: 0.8) |
| `--p-test` | Proportion of sequences for the newly-created test data split (`create_orthogonal_splits_module`). | float (Default: 0.2) |
| `-n`, `--n-splits` | Number of train-test split replicates to create (`create_orthogonal_splits_module`). | positive integer (Default: 1) |
| `-s`, `--seed` | Random seed for creation of homology-aware train-test splits (`create_orthogonal_splits_module`). | positive integer (Default: 21) |
| `--force` | Force overwrite existing `blastn_module` output files. | Boolean (Default: False) |
| `-o`, `--output-dir` | Directory to write the created train-test splits. | string (Default: '.') |

</details>
