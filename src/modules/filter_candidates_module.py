import os
import pandas as pd
import utils.helper_functions as helper
import logging

score_columns = ["id_i","id_j","score"]

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")

    """
    Parse through tabular input file in chunks to keep memory manageable.
    """
    filtered_dfs = []
    try:
        for df in pd.read_csv(args.input_path,names=score_columns,sep="\t",chunksize=50_000):
            df = df[(df["score"] >= args.threshold) & (df["id_i"] != df["id_j"])]
            filtered_dfs.append(df)
    except:
        raise Exception(f"Unable to filter the provided input file ({args.input_path}).")
    
    filtered_df = pd.concat(filtered_dfs)
    filtered_df['pair'] = filtered_df[['id_i','id_j']].apply(lambda x: tuple(sorted(x)),axis=1)
    filtered_df = filtered_df.drop_duplicates(subset='pair')
    filtered_df = filtered_df[["id_i","id_j","score"]]

    if filtered_df.shape[0] == 0:
        logger.warning(f"No sequence pairs with an alignment score above the specified threshold: {args.threshold}")

    filtered_path = os.path.join(args.output_dir,f"hashFrag.similar_pairs.tsv")
    filtered_df.to_csv(filtered_path,sep="\t",index=False)

    logger.info(f"Filtered results written to: {filtered_path}")
    logger.info(f"Module execution completed.\n")