import os
import pandas as pd
from glob import glob
# import utils.helper_functions as helper

def run(args):

    score_columns = ["id_i","id_j","score"]
    score_dfs = []
    pattern = os.path.join(args.score_dir,"*")
    for partitioned_score_path in sorted(glob(pattern)):
        score_df = pd.read_csv(partitioned_score_path,names=score_columns)
        score_df = score_df[
            (score_df["score"] >= args.threshold) & \
            (score_df["id_i"] != score_df["id_j"])
        ]
        score_dfs.append(score_df)
    score_df = pd.concat(score_dfs)
    del score_dfs

    score_df['pair'] = score_df[['id_i', 'id_j']].apply(lambda x: tuple(sorted(x)), axis=1)
    score_df = score_df.drop_duplicates(subset='pair')

    blast_columns = [
        "qseqid","sseqid","pident","length","mismatch",
        "gapopen","qstart","qend","sstart","send","evalue",
        "bitscore","score","positive","gaps"
    ]
    # blast_dfs = []
    pattern = os.path.join(args.blast_dir,"*")
    for partitioned_blast_path in sorted(glob(pattern)):
        partition_label = os.path.basename(partitioned_blast_path)
        if partition_label.endswith(".augmented.tsv.gz"): continue

        partitioned_blast_df = pd.read_csv(partitioned_blast_path,names=blast_columns,sep="\t")
        partitioned_blast_df['pair'] = partitioned_blast_df[['qseqid','sseqid']].apply(lambda x: tuple(sorted(x)), axis=1)
        partitioned_blast_df = partitioned_blast_df.merge(score_df[['pair','smith_waterman_score']],on='pair',how='left')
        partitioned_blast_df = partitioned_blast_df[~partitioned_blast_df["smith_waterman_score"].isna()]

        augmented_blast_path = partitioned_blast_path+".augmented.tsv.gz"
        partitioned_blast_df.to_csv(augmented_blast_path,compression="gzip",sep="\t",index=False)
        del partitioned_blast_df

    pattern = os.path.join(args.blast_dir,"*.augmented.tsv.gz")
    augmented_blast_dfs = []
    for augmented_blast_path in sorted(glob(pattern)):
        blast_df = pd.read_csv(augmented_blast_path,sep="\t")
        augmented_blast_dfs.append(blast_df[["pair","smith_waterman_score"]].copy())
        del blast_df
    
    augmented_blast_df = pd.concat(augmented_blast_dfs)
    augmented_blast_df = augmented_blast_df.drop_duplicates(subset='pair')
    filtered_blast_path = os.path.join(args.blast_dir,"blastn_results.filtered_candidates.tsv.gz")
    augmented_blast_df.to_csv(filtered_blast_path,compression="gzip",sep="\t",index=False)

    print(f"Filtered collision results written to: {args.filtered_blast_path}")
