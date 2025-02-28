import sys
import math
import numpy as np
import pandas as pd
import utils.helper_functions as helper
import logging

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

    logger.info("Stratifying based on the provided pairwise scores...")
    try:
        for df in pd.read_csv(args.input_path,names=score_columns,sep="\t",chunksize=50_000):
            for id_i,id_j,score in zip(df["id_i"],df["id_j"],df["score"]):
                if id_i in idset:
                    scores_dict[id_i].append(score)
                if id_j in idset:
                    scores_dict[id_j].append(score)
    except:
        raise Exception(f"Unable to parse the provided input file ({args.input_path}).")

    scores_dict = { sample_id:np.max(scores) for sample_id,scores in scores_dict.items() if scores }
    scores_df   = pd.DataFrame({"id":scores_dict.keys(),"score":scores_dict.values()})
    scores_df["stratification"] = "orthogonal" # instantiate all sequences as orthogonal

    min_bound = round_to_lower_bound(scores_df["score"].min(),args.step)
    max_bound = round_to_upper_bound(scores_df["score"].max(),args.step)
    for lower_bound in range(min_bound,max_bound,args.step):
        upper_bound = lower_bound + args.step
        label = f"[{lower_bound},{upper_bound})"
        stratified_idset = set(scores_df[(scores_df["score"] >= lower_bound) & (scores_df["score"] < upper_bound)]["id"].tolist())
        scores_df.loc[scores_df["id"].isin(stratified_idset),"stratification"] = label

    scores_df.to_csv(args.output_path,sep="\t",index=False)

    logger.info(f"Stratification results written to: {args.output_path}")
    logger.info(f"Module execution completed.\n")
