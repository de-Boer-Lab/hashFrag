import os
import argparse

from modules.blastn import run as hashFrag_blastn
from modules.stratify_test_split import run as hashFrag_stratify_test_split
    

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

    print("Running stratify_test_split...")
    stratified_path = os.path.join(args.output_dir,f"{label}.stratified_test_split.tsv.gz")
    stratify_test_split_args = argparse.Namespace(
        test_fasta_path=args.test_fasta_path,
        input_path=blast_path,
        mode="lightning",
        gapopen=args.gapopen,
        gapextend=args.gapextend,
        penalty=args.penalty,
        reward=args.reward,
        step=args.step,
        output_path=stratified_path
    )
    hashFrag_stratify_test_split(stratify_test_split_args)

    print("Master command completed successfully.")
