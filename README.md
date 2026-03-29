<img src="https://raw.githubusercontent.com/de-Boer-Lab/hashFrag/main/imgs/hashFrag_logo.png" width=800>

[![Documentation Status](https://readthedocs.org/projects/hashfrag/badge/?version=latest)](https://hashfrag.readthedocs.io/en/latest/)
[![PyPI version](https://img.shields.io/pypi/v/hashfrag)](https://pypi.org/project/hashfrag/)

# Overview

Neural networks have emerged as powerful tools to understand the functional relationship between genomic sequences and various biological processes. However, current practices of training and evaluating models on genomic sequences may fail to account for the widespread homology that permeates the genome. Homology spanning train-test data splits can result in [data leakage](https://en.wikipedia.org/wiki/Leakage_(machine_learning)#:~:text=In%20statistics%20and%20machine%20learning,when%20run%20in%20a%20production), potentially leading to overestimation of model performance and a reduction in model reliability and generalizability.

<img src="https://raw.githubusercontent.com/de-Boer-Lab/hashFrag/main/imgs/hashFrag_workflow_diagram.png">

hashFrag is a scalable command-line tool to help users address homology-based data leakage during model development. The general workflow involves identifying “candidate” pairs of sequences exhibiting high similarity with [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), filtering these candidates based on a specified similarity threshold, and then using the resulting homology information to mitigate the potential occurrences of data leakage in existing or newly-created splits. 

# Documentation

Full documentation is available on [Read the Docs](https://hashfrag.readthedocs.io/en/latest/).

# Paper

Check out our [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.22.634321v1) titled, "*Detecting and avoiding homology-based data leakage in genome-trained sequence models*", for more details.
