<img src="./imgs/hashFrag_logo.png" width=800>

# Overview

Neural networks have emerged as powerful tools to understand the relationship between genomic sequences and various biological processes. However, current practices of training and evaluation of model performance on genomic sequences may fail to account for the widespread homology that permeates the genome, which can manifest as [data leakage](https://en.wikipedia.org/wiki/Leakage_(machine_learning)#:~:text=In%20statistics%20and%20machine%20learning,when%20run%20in%20a%20production). 

<img src="./imgs/hashFrag_workflow_diagram.png">

hashFrag represents a scalable tool that leverages [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to allow users to account for homology-based data leakage when developing models. The general workflow involves calling the BLAST algorithm on input sequences to identify candidate pairs exhibiting high similarity, filtering candidates based on a specified threshold of similarity, and then using the homology information to mitigate the potential bias caused by data leakage.

The primary use cases of hashFrag include:

1. Removal of homology spanning existing data splits (e.g., chromosomal splits) provided by the user.
2. Stratification of the test split based on their maximal alignment score with respect to the training split.
3. Creation of orthogonal train-test splits from input sequences.

Note that when existing data splits are provided, comparisons are only assessed for pairs of sequences across the data splits. 

We use the alignment score to quantify the degree of homology between a pair of sequences. By default, the alignment score will be derived from the top BLAST alignment result for a pair of sequences. However, users also have the option to precompute alternative alignment scores, such as Smith-Waterman.

Given that the precise definition of homology may depend on the dataset or hypothesis under test, we aim to provide general guidelines for users to choose an appropriate threshold for their purposes.

# Installation

placeholder

# Basic usage

## Existing data splits

Filter all test sequences exhibiting homology with any sequences in the train split according to a specified alignment score threshold.
```
hashFrag filter_test_split \
--train_fasta_path K562.sample_8000.train.fa.gz \
--test_fasta_path K562.sample_2000.test.fa.gz \
--threshold 60
```

Stratify the test split into an arbitrary number of levels based on their maximum alignment scores to the train split. Note that the sizes of each level will not necessarily be balanced. 
```
hashFrag stratify_test_split 
--train_fasta_path K562.sample_8000.train.fa.gz \
--test_fasta_path K562.sample_2000.test.fa.gz
```

## Creating data splits

Create homology-aware (i.e., orthogonal) train-test data splits according to a specified alignment score threshold.
```
# hashFrag create_orthogonal_splits \
--fasta_path K562.sample_10000.fa.gz \
--threshold 60
```

# Advanced usage

The basic usage commands are composed of modules for each step in the process. For added flexibility, users can directly call these modules. This can enable the use of precomputed pairwise scores (such as the optimal [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) alignment score), rather than the heuristic BLASTn alignment scores.

For a full breakdown of available modules, please see the notebooks provided in the `/tutorials` directory in this repository.


# Paper

placeholder