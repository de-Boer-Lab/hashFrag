import os
import gzip
import pandas as pd
import multiprocessing as mp

from functools import partial
from Bio import SeqIO,pairwise2
from itertools import combinations
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description="")
    parser.add_argument('-i','--fasta_path_i',type=str,required=True,help="")
    parser.add_argument('-j','--fasta_path_j',type=str,default=None,help="")
    parser.add_argument('-b','--blast_path',type=str,default=None,help="")
    parser.add_argument('-m','--match',type=int,default=1,help="")
    parser.add_argument('-x','--mismatch',type=int,default=1,help="")
    parser.add_argument('-g','--gap_open',type=int,default=2,help="")
    parser.add_argument('-e','--gap_extend',type=int,default=1,help="")
    parser.add_argument('-n','--n_processes',type=int,default=1,help="")
    parser.add_argument('-v','--verbose',type=bool,default=False,help="")
    parser.add_argument('-o','--outpath',type=str,required=True,help="")
    return parser.parse_args()

def load_fasta(path):
    fasta_dict = {}
    with gzip.open(path,'rt') as handle:
        for record in SeqIO.parse(handle,'fasta'):
            fasta_dict[record.description] = str(record.seq)
    return fasta_dict

def compute_sw_score(pair,match,mismatch,gap_open,gap_extend):
    """ Bio.pairwise2 faster than Bio.Align.PairwiseAligner if only computing alignment score
    """
    (id_i,seq_i),(id_j,seq_j) = pair
    score = pairwise2.align.localms(seq_i,seq_j,match,-mismatch,-gap_open,-gap_extend,score_only=True)
    return (id_i,id_j,score)

def compute_pairwise_SW_scores():
    args = parse_arguments()

    if args.verbose:
        print("Computing pairwise Smith-Waterman alignment scores:")
        print(f"\tFASTA file (i): {args.fasta_path_i}")
        print(f"\tFASTA file (j): {args.fasta_path_j}")
        print(f"\tBLAST output file: {args.blast_path}")
        print(f"\tSmith-Waterman alignment score parameters:")
        print(f"\t  Match: {args.match}")
        print(f"\t  Mismatch: {args.mismatch}")
        print(f"\t  Gap open: {args.gap_open}")
        print(f"\t  Gap extend: {args.gap_extend}")
        print(f"\tNumber of CPUs: {args.n_processes}")
        print(flush=True)

    if args.fasta_path_j is not None and os.path.exists(args.fasta_path_j):
        if args.verbose:
            print("[MODE] Computing pairwise scores between the two sequence files provided.",flush=True)
        def pairwise_generator(fasta_path_i,fasta_path_j):
            fasta_dict_i = load_fasta(fasta_path_i)
            fasta_dict_j = load_fasta(fasta_path_j)
            for id_i,seq_i in fasta_dict_i.items():
                for id_j,seq_j in fasta_dict_j.items():
                    yield ((id_i,seq_i),(id_j,seq_j))

        pairs = pairwise_generator(args.fasta_path_i,args.fasta_path_j)

    elif args.blast_path is not None and os.path.exists(args.blast_path):
        if args.verbose:
            print("[MODE] Computing pairwise scores of BLAST candidate pairs.",flush=True)
        def pairwise_generator(fasta_path,blast_path,chunksize=250_000):
            fasta_dict = load_fasta(fasta_path)
            columns = [
                "qseqid","sseqid","pident","length","mismatch",
                "gapopen","qstart","qend","sstart","send","evalue",
                "bitscore","score","positive","gaps"
            ]
            for chunk_df in pd.read_csv(blast_path,sep="\t",names=columns,chunksize=chunksize):
                for query_id,candidates in chunk_df.groupby("qseqid")["sseqid"].apply(set).items():
                    for subject_id in candidates:
                        yield ((query_id,fasta_dict[query_id]),(subject_id,fasta_dict[subject_id]))

        pairs = pairwise_generator(args.fasta_path_i,args.blast_path)

    else:
        if args.verbose:
            print("[MODE] Computing pairwise scores for the provided sequence file.",flush=True)
        def parse_fasta_to_generator(path):
            with gzip.open(path, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    yield record.description, str(record.seq)

        seqs  = parse_fasta_to_generator(args.fasta_path_i)
        pairs = combinations(seqs,2)

    with mp.Pool(args.n_processes) as p:
        compute_sw_score_partial = partial(
            compute_sw_score,
            match=args.match,
            mismatch=args.mismatch,
            gap_open=args.gap_open,
            gap_extend=args.gap_extend
        )
        results = p.map(compute_sw_score_partial,pairs)

    with gzip.open(args.outpath,"wt") as handle:
        for (id_i,id_j,sw) in results:
            handle.write(f"{id_i}\t{id_j}\t{sw}\n")

    if args.verbose:
        print(f"SW scores written file: {args.outpath}",flush=True)

if __name__ == "__main__":
    compute_pairwise_SW_scores()
