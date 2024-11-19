import numpy as np
import pandas as pd
import utils.helper_functions as helper

def run(args):

    collisions_dict = helper.load_collisions_dictionary(args.collisions_path)
    collisions_dict = { seq_id:set(collisions) for seq_id,collisions in collisions_dict.items() }

    df = pd.read_csv(args.score_matrix_path,compression="gzip",index_col=0)
    df = df.stack().reset_index()
    df.columns = ["id_i","id_j","sw"]
    collided = []
    for _,row in df.iterrows():
        id_i = row["id_i"]
        id_j = row["id_j"]
        collided.append(int((id_i in collisions_dict[id_j]) or (id_j in collisions_dict[id_i])))
    df["collided"] = collided

    results_dict = {
        "threshold":[],
        "TP":[],
        "FP":[],
        "TN":[],
        "FN":[],
        "recall":[],
        "FPR":[]
    }
    for threshold in np.arange(df["sw"].min(),df["sw"].max()+args.step,args.step):
        above_threshold = df["sw"] >= threshold
        below_threshold = df["sw"] < threshold
        tp = np.sum(above_threshold & (df["collided"] == 1))  # True Positives: number of hits LSH captures
        fn = np.sum(above_threshold & (df["collided"] == 0))  # False Negatives: number of hits LSH misses
        fp = np.sum(below_threshold & (df["collided"] == 1))  # False Positives: number of non-hits LSH successfully misses
        tn = np.sum(below_threshold & (df["collided"] == 0))  # True Negatives:  number of non-hits that LSH detects
        results_dict["threshold"].append(threshold)
        results_dict["TP"].append(tp)
        results_dict["FN"].append(fn)
        results_dict["FP"].append(fp)
        results_dict["TN"].append(tn)
        results_dict["recall"].append(tp/(tp+fn) if (tp+fn) > 0 else 0)
        results_dict["FPR"].append(fp/(fp+tn) if (fp+tn) > 0 else 0)

    results_df = pd.DataFrame(results_dict)
    results_df.to_csv(args.output_path,compression="gzip",index=False)

    print(f"Recall calculation results written to: {args.output_path}")
