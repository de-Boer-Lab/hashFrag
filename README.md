<img src="./imgs/hashFrag_logo.png" width=800>

---
# Overview

Neural networks have emerged as powerful tools to understand the relationship between genomic sequences and various biological processes. However, current practices of training and testing model performance on genomic sequences may fail to account for the widespread homology that permeates the genome, which can lead to data leakage.

<img src="./imgs/hashFrag_workflow_diagram.png">

hashFrag represents a scalable tool that leverages [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to allow users to account for homology-based data leakage when training and evaluating models. Homology can either be removed from existing data splits (e.g., chromosomal splits) provided by the user, or orthogonal train-test splits can be created from input sequences.  

---
# Installation

placeholder

---
# Basic usage

## Existing data splits

Filter all test sequences exhibiting homology with any sequences in the train split according to a specified alignment score threshold.
```
# placeholder
```

Stratify the test split into an arbitrary number of levels based on their maximum alignment scores to the train split.
```
# placeholder
```

## Creating data splits

Create homology-aware (i.e., orthogonal) train-test data splits according to a specified alignment score threshold.
```
# placeholder
```

---
# Advanced usage

In addition to the basic usage commands, users have the option of added flexibility by calling the subcomponents directly.

---
# Paper

placeholder