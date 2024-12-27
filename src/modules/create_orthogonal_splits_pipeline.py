import os
import argparse

from modules.blastn import run as hashFrag_blastn
from modules.filter_false_positives import run as hashFrag_filter_candidates
from modules.identify_homologous_groups import run as hashFrag_identify_groups
from modules.create_orthogonal_splits import run as hashFrag_create_orthosplits

def run(args):

    label = "hashFrag_lightning"

    print("Running blastn...")
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

    print("Running filter_false_positives...")
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
    homology_path = os.path.join(args.output_dir,f"hashFrag_lightning.homologous_groups.tsv")

    print("Running identify_homologous_groups...")
    identify_groups_args = argparse.Namespace(
        hits_path=hits_path,
        output_path=homology_path
    )
    hashFrag_identify_groups(identify_groups_args)

    print("Running create_orthogonal_splits...")
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

    print("Master command completed successfully.")
