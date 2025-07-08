import os
import gzip
import random
import logging
from collections import defaultdict
import utils.helper_functions as helper

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")

    if args.p_train+args.p_test != 1:
        raise Exception("Intestid probabilities specified!")

    random.seed(args.seed)

    """
    homologous_groups will be a list of sublists. A distinct group of sequences sharing 
    homology make up each sublist. 

    hit_idset will be used to efficiently identify sequences are orthogonal or share
    homology with another sequence in the population.
    """
    N = 0
    homologous_groups = defaultdict(set)
    with open(args.homology_path,"r") as handle:
        for line in handle:
            N += 1
            sample_id,group_id = line.strip().split("\t")
            homologous_groups[group_id].add(sample_id)
    homologous_groups = list(homologous_groups.values())

    n_train = int(round(N*args.p_train,0))
    n_test  = int(round(N*args.p_test,0))

    if n_train+n_test != N:
        raise Exception(f"Combined train ({n_train}) and test ({n_test}) size doesn't add up to expected total ({N}).")

    logger.info(f"Creating {args.n_splits} orthogonal splits in directory: {args.output_dir}")
    for i in range(args.n_splits):
        split = f"split_{i+1:0>3}"

        groups_ = homologous_groups.copy()
        test_split = set()
        while len(test_split) < n_test:
            idx = random.randint(0,len(groups_)-1)
            group_ = groups_.pop(idx)
            test_split.update(group_)

        train_split = set()
        for group_ in groups_:
            train_split.update(group_)
        
        test_split = sorted(test_split)
        train_split = sorted(train_split)

        filename = f"hashFrag.train_{len(train_split)}.test_{len(test_split)}.{split}.tsv"
        outpath  = os.path.join(args.output_dir,filename)
        helper.write_splits_to_tsv(train_split,test_split,outpath)

    logger.info(f"Module execution completed.\n")
