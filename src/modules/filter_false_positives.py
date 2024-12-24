import os
import pandas as pd
from glob import glob
import utils.helper_functions as helper

blast_columns = [
    "qseqid","sseqid","pident","length","mismatch",
    "gapopen","qstart","qend","sstart","send","evalue",
    "bitscore","uncorrected_blast_score","positive","gaps"
]

score_columns = ["id_i","id_j","score"]

def run(args):

    filtered_dfs = []

    if args.mode == "lightning": # Expects a BLAST output file

        try:
            for blast_df in pd.read_csv(args.input_path,names=blast_columns,sep="\t",chunksize=50_000):
                helper.blast_file_validation(blast_df,blast_columns)
                blast_df["corrected_blast_score"] = (
                    args.reward*blast_df["positive"] +
                    args.penalty*blast_df["mismatch"] -
                    args.gapopen*blast_df["gapopen"] -
                    args.gapextend*(blast_df["gaps"]-blast_df["gapopen"])
                )
                blast_df = blast_df[
                    (blast_df["corrected_blast_score"] >= args.threshold) & \
                    (blast_df["qseqid"] != blast_df["sseqid"])
                ]
                blast_df = blast_df[["qseqid","sseqid","corrected_blast_score"]]
                blast_df.columns = ["id_i","id_j","score"]
                filtered_dfs.append(blast_df)
        except:
            raise Exception("Error: Unable to read BLAST file (lightning mode).")
    
    elif args.mode == "pure": # Expects a custom tsv file
    
        try:
            for score_df in pd.read_csv(args.input_path,names=score_columns,sep="\t",chunksize=50_000):
                score_df = score_df[
                    (score_df["score"] >= args.threshold) & \
                    (score_df["id_i"] != score_df["id_j"])
                ]
                filtered_dfs.append(score_df)
        except:
            raise Exception("Error: Unable to read custom scores (pure mode).")
    
    else:
        raise Exception("Error: invalid 'mode' specified. Permissible values include: {'lightning', 'pure'}.")

    filtered_df = pd.concat(filtered_dfs)
    filtered_df['pair'] = filtered_df[['id_i','id_j']].apply(lambda x: tuple(sorted(x)),axis=1)
    filtered_df = filtered_df.drop_duplicates(subset='pair')
    filtered_df = filtered_df[["id_i","id_j","score"]]
    filtered_blast_path = os.path.join(args.out_dir,f"hashFrag_{args.mode}.similar_pairs.tsv.gz")
    filtered_df.to_csv(filtered_blast_path,compression="gzip",sep="\t",index=False)

    print(f"Filtered results written to: {filtered_blast_path}")
