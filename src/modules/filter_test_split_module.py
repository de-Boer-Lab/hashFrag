import gzip
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import utils.helper_functions as helper
import logging

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")
    hits_dict = defaultdict(set)
    for chunk_df in pd.read_csv(args.hits_path,sep="\t",chunksize=250_000):
        for qseqid, sseqid in zip(chunk_df["id_i"],chunk_df["id_j"]):
            hits_dict[qseqid].add(sseqid)
            hits_dict[sseqid].add(qseqid) # symmetric dictionary

    train_idset = set(helper.load_fasta_as_dictionary(args.train_fasta_path))
    test_fasta_dict = helper.load_fasta_as_dictionary(args.test_fasta_path)
    test_ids = list(test_fasta_dict)
    n_test = len(test_ids)
    filtered_test_ids = []
    while test_ids:
        sample_id = test_ids.pop()

        """
        Any intersection indicates homology spanning the train-test splits. Thus,
        we only include test sequences that are either orthogonal (empty set in hits dict)
        or whose group of homologous sequences is disjoint from the train split.
        """
        if hits_dict[sample_id].isdisjoint(train_idset):
            filtered_test_ids.append(sample_id)
    
    logger.info(f"{n_test-len(filtered_test_ids)} sequences filtered from test split.")

    filtered_test_fasta_path = args.test_fasta_path.replace(".fa",".filtered.fa")
    with gzip.open(filtered_test_fasta_path,"wt") as handle:
        for sample_id in filtered_test_ids:
            handle.write(f">{sample_id}\n{test_fasta_dict[sample_id]}\n")

    logger.info(f"Filtered results written to: {filtered_test_fasta_path}")
    logger.info(f"Module execution completed.\n")
