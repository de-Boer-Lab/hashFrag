import os
import gzip
import random
import logging
from collections import defaultdict
import utils.helper_functions as helper

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")

    random.seed(args.seed)

    """
    homologous_groups will be a list of sublists. A distinct group of sequences sharing 
    homology make up each sublist. 
    """
    N = 0
    homologous_groups = defaultdict(set)
    with open(args.homology_path,"r") as handle:
        for line in handle:
            N += 1
            sample_id,group_id = line.strip().split("\t")
            homologous_groups[group_id].add(sample_id)
    homologous_groups = list(homologous_groups.values())
    homologous_groups = sorted(homologous_groups, key=len, reverse=True)
    max_group_size = max(len(g) for g in homologous_groups)

    if N//max_group_size > args.folds:
        logger.warning(f"There exist(s) homologous groups of larger size than the expected fold size. Resulting folds may be imbalanced!")

    """
    Greedy approach to construct folds from homologous groups. 
    Homologous groups (from largest to smallest in terms of group size)
    will be greedily added to smallest folds. This should obtain roughly
    equal sizes of folds.
    """
    logger.info(f"Creating {args.folds} orthogonal folds...")
    folds = [[] for _ in range(args.folds)]
    foldsizes = [0]*args.folds
    for group in homologous_groups:
        i = foldsizes.index(min(foldsizes))
        folds[i].extend(group)
        foldsizes[i] += len(group)

    filename = f"hashFrag.{args.folds}_orthogonal_folds.tsv"
    outpath  = os.path.join(args.output_dir,filename)
    helper.write_folds_to_tsv(folds,outpath)

    logger.info(f"Orthogonal folds written to {outpath}")
    logger.info(f"Module execution completed.\n")
