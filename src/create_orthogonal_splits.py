import os
import gzip
import random
from collections import defaultdict
import utils.helper_functions as helper

def run(args):

    if args.p_train+args.p_test != 1:
        raise Exception("Intestid probabilities specified!")

    random.seed(args.seed)

    hit_idset = set()
    homologous_groups = defaultdict(set)

    with open(args.homology_path,"r") as handle:
        for line in handle:
            sample_id,group_id = line.strip().split(",")
            hit_idset.add(sample_id)
            homologous_groups[group_id].add(sample_id)

    homologous_groups = list(homologous_groups.values())

    ids = helper.load_fasta_ids(args.fasta_path)
    n = len(ids)

    n_train = int(round(n*args.p_train,0))
    n_test  = int(round(n*args.p_test,0))

    if n_train+n_test != n:
        raise Exception(f"Combined train ({n_train}) and test ({n_test}) size doesn't add up to expected total ({n}).")

    print("Writing splits...")
    for i in range(args.n_splits):
        split = f"split_{i+1:0>3}"

        ids_ = ids.copy()

        train_split = set()
        while len(train_split) < n_train:
            sample_id = random.choice(ids_)
            complementary_sample_id = helper.get_complementary_id(sample_id)

            train_split.update([sample_id,complementary_sample_id])
            ids_.remove(sample_id)
            ids_.remove(complementary_sample_id)

            if sample_id in hit_idset:
                for homologous_group in homologous_groups:
                    if sample_id in homologous_group:
                        ids_ = [id_ for id_ in ids_ if id_ not in homologous_group]
                        train_split.update(homologous_group)

        train_split = sorted(train_split)
        test_split  = sorted(ids_)

        filename = f"hashFrag.train_{len(train_split)}.test_{len(test_split)}.{split}.csv.gz"
        outpath  = os.path.join(args.out_dir,filename)
        print(" ",outpath,flush=True)
        with gzip.open(outpath,"wt") as handle:
            handle.write("id,split\n")
            for sample_id in train_split:
                handle.write(f"{sample_id},train\n")
            for sample_id in test_split:
                handle.write(f"{sample_id},test\n")
