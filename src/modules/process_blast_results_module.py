import os
import pandas as pd
import utils.helper_functions as helper
import logging
from collections import defaultdict

blast_columns = [
    "qseqid","sseqid","pident","length","mismatch",
    "gapopen","qstart","qend","sstart","send","evalue",
    "bitscore","uncorrected_blast_score","positive","gaps"
]

score_columns = ["id_i","id_j","score"]

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")

    top_scores_dict = defaultdict(float)
    for blast_df in pd.read_csv(args.blastn_path,names=blast_columns,sep="\t",chunksize=25_000):
        helper.blast_file_validation(blast_df,blast_columns)
        blast_df["corrected_blast_score"] = ( # corrected BLAST alignment score
            args.reward*blast_df["positive"] +
            args.penalty*blast_df["mismatch"] -
            args.gapopen*blast_df["gapopen"] -
            args.gapextend*(blast_df["gaps"]-blast_df["gapopen"])
        )
        for qseqid,sseqid,score in zip(blast_df['qseqid'],blast_df['sseqid'],blast_df['corrected_blast_score']):
            pair = (qseqid,sseqid)
            if score > top_scores_dict[pair]:
                top_scores_dict[pair] = score

    with open(args.processed_blastn_path,"w") as handle:
        for (qseqid,sseqid),score in top_scores_dict.items():
            handle.write(f"{qseqid}\t{sseqid}\t{score}\n")

    logger.info(f"Processed BLASTn results written to: {args.processed_blastn_path}")
    logger.info(f"Module execution completed.\n")