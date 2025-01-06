<img src="./imgs/hashFrag_logo.png" width=800>

# Overview

Neural networks have emerged as powerful tools to understand the functional relationship between genomic sequences and various biological processes. However, current practices of training and evaluating models on genomic sequences may fail to account for the widespread homology that permeates the genome. Homology spanning test-train data splits can result in [data leakage](https://en.wikipedia.org/wiki/Leakage_(machine_learning)#:~:text=In%20statistics%20and%20machine%20learning,when%20run%20in%20a%20production), potentially leading to overestimation of model performance and a reduction in model reliability and generalizability.

<img src="./imgs/hashFrag_workflow_diagram.png">

hashFrag is a scalable command-line tool to help users address homology-based data leakage during model development. The general workflow involves identifying “candidate” pairs of sequences exhibiting high similarity with [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), filtering these candidates based on a specified similarity threshold, and then using the resulting homology information to mitigate the potential occurrences of data leakage in existing or newly-created splits. 

We utilize local alignment scores to quantify the degree of homology between a pair of sequences. By default, the alignment score will be derived from the top BLAST alignment result for a pair of sequences (see Basic usage). However, users also have the option to provide precomputed alignment scores for added control over the homology search process (see Advanced usage). 

Because the precise definition of homology depends on the dataset and/or hypothesis under test, we also aim to provide general guidelines for users to choose an appropriate similarity threshold for their purposes. Notably, defining homology in terms of alignment scores requires the specification of scoring parameters (e.g., mismatch, gap open, and gap extension penalty scores and match reward), and changing these parameters can drastically impact the identification process of homology. Please see permissible scoring parameter combinations for the BLASTn algorithm [here](https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/) (Table D1).

# Installation

It is recommended to execute hashFrag in a `conda` or `virtualenv` environment with Python version 3.10.

Clone the repository using the following command:
```
git clone https://github.com/de-Boer-Lab/hashFrag.git
```

Install dependencies located in the `requirements.txt` file with the folowing command:
```
pip install -r requirements.txt
```

Export the source directory to your `PATH` with the following command:
```
export PATH="$PATH:./hashFrag/src"
```

## BLAST+ download

To install the BLAST command-line tool, follow the NCBI BLAST Command Line Applications User Manual found [here](https://www.ncbi.nlm.nih.gov/books/NBK569861/).

Directly download the BLAST+ package executables for different operating systems at the following NCBI FTP page:
> https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Follow the instructions to extract the downloaded file and export binaries to your `PATH`. 

Verify that `blastn` and `makeblastdb` commands were installed succesfully:
```
blastn -version
makeblastdb -version
```

# Basic usage

## Existing data splits

Existing train-test data splits can be handled by providing two separate FASTA files to hashFrag as input. This limits the homology search process to inter-data split comparisons. Specifically, a BLAST database is constructed over the train sequences and the test sequences are queried against this database to identify pairs with high local alignment scores. 

> Filter sequences in the test split exhibiting homology with any sequences in the train split. This requires specification of an alignment score threshold to define homology between sequences.
```
hashFrag filter_existing_splits \
--train_fasta_path example_train_split.fa.gz \
--test_fasta_path example_test_split.fa.gz \
-t 60 \
-o filter_existing_splits.work
```

> Stratify the test split sequences into an arbitrary number of levels based on their maximum alignment scores to the train split sequences. 
```
hashFrag stratify_test_split \
--train_fasta_path example_train_split.fa.gz \
--test_fasta_path example_test_split.fa.gz \
-o stratify_test_split.work
```
Note that the sizes of each stratified level will not necessarily be balanced. This can be useful to better understand a model’s behaviour over test splits at varying levels of orthogonality to the sequences the model was trained on.

## Creating data splits

When a single FASTA file is provided as input, hashFrag will characterize homology for all pairwise comparisons. This involves constructing a BLAST database over all sequences in the population and then subsequently querying each sequence to the databae. 

> Create homology-aware (i.g., orthogonal) train-test data splits. This requires specification of an alignment score threshold to define homology between sequences.
```
hashFrag create_orthogonal_splits \
-f example_full_dataset.fa.gz \
-t 60 \
-o create_orthogonal_splits.work
```
The creation of orthogonal train-test splits involves encoding the homologous relationships between sequences as a sparse adjacency matrix (unweighted in accordance with the alignment score threshold). A graph representation of the adjancency matrix is constructed, and then distinct groups of homologous sequences can be identified by finding disconnected subgraphs. From the homology cluster information, splits with no leakage can be created proportionally. 

# Advanced usage

The basic usage commands are implemented as pipelines that execute a series of modules. To provide users with additional control and flexibility over this homology search process, users can directly call these modules.

| Pipeline                   | Modules            |
|----------------------------|-----------------------|
| `filter_existing_splits`   | `blastn_module`, `filter_candidates_module`, `filter_test_split_module` |
| `stratify_test_split`      | `blastn_module`, `stratify_test_split_module` |
| `create_orthogonal_splits` | `blastn_module`, `filter_candidates_module`, `identify_homologous_groups_module`, `create_orthogonal_splits_module` |

The main advantage of calling modules individually is that it enables the use of precomputed pairwise scores. For example, after identifying candidate pairs of sequences with the BLAST algorithm, instead of using BLAST-derived alignment scores users can provide the optimal [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) local alignment scores for candidate pairs. This may be particularly important for improving recall in cases where the degree of homology between a pair of sequences is "modest". 

When providing precomputed pairwise scores to hashFrag, the expected format is a tab-delimited file with 3 columns: `id_i`, `id_j`, and `score`. 

`example.tsv`
```
seq_A	seq_B	60
seq_C	seq_D	85
seq_E	seq_A	100
...     ...     ...
```

For a full breakdown of available modules, please see the notebooks provided in the `/tutorials` directory.

# Paper

placeholder