import os
import logging
import argparse

from modules.blastn_module import run as hashFrag_blastn
from modules.process_blast_results_module import run as hashFrag_process_blast_results
from modules.filter_candidates_module import run as hashFrag_filter_candidates
from modules.identify_homologous_groups_module import run as hashFrag_identify_groups
from modules.create_orthogonal_splits_module import run as hashFrag_create_orthosplits

def run(args):

    logger = logging.getLogger("pipeline")
    logger.info("Initializing `create_orthogonal_splits` pipeline.\n")

    label = "hashFrag"

    blastn_args = argparse.Namespace(
        fasta_path=args.fasta_path,
        train_fasta_path=None,
        test_fasta_path=None,
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

    blast_path = os.path.join(args.output_dir,f"{label}.blastn.out")
    processed_blast_path = os.path.join(args.output_dir,f"{label}.blastn.processed.tsv")
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

    # filter_candidates_args = argparse.Namespace(
    #     input_path=processed_blast_path,
    #     threshold=args.threshold,
    #     output_dir=args.output_dir
    # )
    # hashFrag_filter_candidates(filter_candidates_args)
    # hits_path = os.path.join(args.output_dir,f"hashFrag.similar_pairs.tsv")
    homology_path = os.path.join(args.output_dir,f"hashFrag.homologous_groups.tsv")
    identify_groups_args = argparse.Namespace(
        hits_path=processed_blast_path,
        threshold=args.threshold,
        output_path=homology_path
    )
    hashFrag_identify_groups(identify_groups_args)

    create_orthosplits_args = argparse.Namespace(
        fasta_path=args.fasta_path,
        homology_path=homology_path,
        p_train=args.p_train,
        p_test=args.p_test,
        n_splits=args.n_splits,
        seed=args.seed,
        output_dir=args.output_dir
    )
    hashFrag_create_orthosplits(create_orthosplits_args)

    logger.info("Completed execution of the pipeline to create homology-aware data splits!")
