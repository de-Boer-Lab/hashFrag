import os
import logging
import argparse

import utils.helper_functions as helper

from modules.blastn_module import run as hashFrag_blastn
from modules.process_blast_results_module import run as hashFrag_process_blast_results
from modules.stratify_test_split_module import run as hashFrag_stratify_test_split
    

def run(args):

    logger = logging.getLogger("pipeline")
    logger.info("Initializing `stratify_test_split` pipeline.\n")

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
        exec_makeblastdb_only=args.exec_makeblastdb_only,
        skip_revcomp=args.skip_revcomp,
        xdrop_ungap=args.xdrop_ungap,
        xdrop_gap=args.xdrop_gap,
        xdrop_gap_final=args.xdrop_gap_final,
        evalue=args.evalue,
        dust=args.dust,
        blastdb_args=args.blastdb_args,
        blastdb_label=label,
        blastn_args=args.blastn_args,
        threads=args.threads,
        force=args.force,
        output_dir=args.output_dir
    )
    hashFrag_blastn(blastn_args)

    query_label  = helper.extract_fasta_file_label(args.test_fasta_path)
    blast_path  = os.path.join(args.output_dir,f"{query_label}.blastn.out")
    processed_blast_path = os.path.join(args.output_dir,f"{query_label}.blastn.processed.tsv")
    process_blast_results_args = argparse.Namespace(
        blastn_path=blast_path,
        word_size=args.word_size,
        gapopen=args.gapopen,
        gapextend=args.gapextend,
        penalty=args.penalty,
        reward=args.reward,
        processed_blastn_path=processed_blast_path
    )
    hashFrag_process_blast_results(process_blast_results_args)

    stratified_path = os.path.join(args.output_dir,f"{label}.stratified_test_split.tsv")
    stratify_test_split_args = argparse.Namespace(
        test_fasta_path=args.test_fasta_path,
        input_path=processed_blast_path,
        step=args.step,
        output_path=stratified_path
    )
    hashFrag_stratify_test_split(stratify_test_split_args)

    logger.info("Completed execution of the pipeline to stratify the test split by alignment score!")
