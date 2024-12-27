import os
import argparse

from modules.blastn import run as hashFrag_blastn
from modules.filter_candidates import run as hashFrag_filter_candidates
from modules.filter_test_split import run as hashFrag_filter_test_split
    

def run(args):

    label = "hashFrag_lightning"

    print("Running blastn...")
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

    print("Running filter_candidates...")
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

    print("Running filter_test_split...")
    filter_test_split_args = argparse.Namespace(
        train_fasta_path=args.train_fasta_path,
        test_fasta_path=args.test_fasta_path,
        hits_path=hits_path
    )
    hashFrag_filter_test_split(filter_test_split_args)

    print("Master command completed successfully.")
