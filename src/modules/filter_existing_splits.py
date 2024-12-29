import os
import logging
import argparse

from modules.blastn_module import run as hashFrag_blastn
from modules.filter_candidates_module import run as hashFrag_filter_candidates
from modules.filter_test_split_module import run as hashFrag_filter_test_split

def run(args):

    logger = logging.getLogger("pipeline")
    logger.info("Initializing `filter_existing_splits` pipeline.\n")

    label = "hashFrag"

    blastn_args = argparse.Namespace(
        fasta_path=None,
        train_fasta_path=args.train_fasta_path,
        test_fasta_path=args.test_fasta_path,
        word_size=args.word_size,
        gapopen=args.gapopen,
        gapextend=args.gapextend,
        penalty=args.penalty,
        reward=args.reward,
        max_target_seqs=args.max_target_seqs,
        e_value=args.e_value,
        dust=args.dust,
        blastdb_args=args.blastdb_args,
        blastdb_label=label,
        blastn_args=args.blastn_args,
        threads=args.threads,
        output_dir=args.output_dir
    )
    hashFrag_blastn(blastn_args)
    blast_path = os.path.join(args.output_dir,f"{label}.blastn.out")

    filter_candidates_args = argparse.Namespace(
        input_path=blast_path,
        mode="lightning",
        gapopen=args.gapopen,
        gapextend=args.gapextend,
        penalty=args.penalty,
        reward=args.reward,
        threshold=args.threshold,
        output_dir=args.output_dir
    )
    hashFrag_filter_candidates(filter_candidates_args)
    hits_path = os.path.join(args.output_dir,f"hashFrag_lightning.similar_pairs.tsv.gz")

    filter_test_split_args = argparse.Namespace(
        train_fasta_path=args.train_fasta_path,
        test_fasta_path=args.test_fasta_path,
        hits_path=hits_path
    )
    hashFrag_filter_test_split(filter_test_split_args)

    logger.info("Completed execution of the pipeline to filter existing data splits!")
