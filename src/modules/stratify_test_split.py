import sys
import math
import numpy as np
import pandas as pd
import utils.helper_functions as helper
import logging

blast_columns = [
    "qseqid","sseqid","pident","length","mismatch",
    "gapopen","qstart","qend","sstart","send","evalue",
    "bitscore","uncorrected_blast_score","positive","gaps"
]

score_columns = ["id_i","id_j","score"]

def round_to_lower_bound(number,step):
    return int(math.floor(number/step)*step)

def round_to_upper_bound(number,step):
    return int(math.ceil(number/step)*step)

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")

    ids = helper.load_fasta_ids(args.test_fasta_path)
    idset = set(ids)
    scores_dict = { sample_id:[] for sample_id in ids }

    if args.mode == "lightning": # Expects a BLAST output file
        logger.info("Stratifying based on corrected BLAST alignment scores (lightning mode).")
        try:
            for blast_df in pd.read_csv(args.input_path,names=blast_columns,sep="\t",chunksize=50_000):
                helper.blast_file_validation(blast_df,blast_columns)
                blast_df["corrected_blast_score"] = (
                    args.reward*blast_df["positive"] +
                    args.penalty*blast_df["mismatch"] -
                    args.gapopen*blast_df["gapopen"] -
                    args.gapextend*(blast_df["gaps"]-blast_df["gapopen"])
                )
                for sample_id,score in zip(blast_df["qseqid"],blast_df["corrected_blast_score"]):
                    if sample_id not in idset:
                        raise Exception("Error: Sequence ID not found in test split.")
                    scores_dict[sample_id].append(score)
        except:
            raise Exception("Error: Unable to read BLAST file (lightning mode).")
    
    elif args.mode == "pure": # Expects a custom tsv file
        logger.info("Stratifying based on precomputed alignment scores (pure mode).")
        try:
            for score_df in pd.read_csv(args.input_path,names=score_columns,sep="\t",chunksize=50_000):
                for sample_id,score in zip(score_df["id_i"],score_df["score"]):
                    if sample_id not in idset:
                        raise Exception("Error: Sequence ID not found in test split.")
                    scores_dict[sample_id].append(score)
        except:
            raise Exception("Error: Unable to read custom scores (pure mode).")
    
    else:
        raise Exception("Error: invalid 'mode' specified. Permissible values include: {'lightning', 'pure'}.")

    scores_dict = { sample_id:np.max(scores) for sample_id,scores in scores_dict.items() }
    scores_df   = pd.DataFrame({"id":scores_dict.keys(),"score":scores_dict.values()})
    scores_df["stratification"] = "orthogonal"

    min_bound = round_to_lower_bound(scores_df["score"].min(),args.step)
    max_bound = round_to_upper_bound(scores_df["score"].max(),args.step)
    for lower_bound in range(min_bound,max_bound,args.step):
        upper_bound = lower_bound + args.step
        label = f"[{lower_bound},{upper_bound})"
        stratified_idset = set(scores_df[(scores_df["score"] >= lower_bound) & (scores_df["score"] < upper_bound)]["id"].tolist())
        scores_df.loc[scores_df["id"].isin(stratified_idset),"stratification"] = label

    scores_df.to_csv(args.output_path,compression="gzip",sep="\t",index=False)

    logger.info(f"Stratification results written to: {args.output_path}")
    logger.info(f"Module execution completed.\n")
