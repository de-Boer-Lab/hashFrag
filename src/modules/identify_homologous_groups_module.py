import sys
import pandas as pd
from collections import defaultdict
import utils.helper_functions as helper
import logging

def get_symmetric_adjacency_dict(adjacency_dict):
    symmetric_adjacency_dict = defaultdict(set)
    for query_id,candidates in adjacency_dict.items():
        complementary_query_id = helper.get_complementary_id(query_id)

        query_set = {query_id,complementary_query_id}
        symmetric_adjacency_dict[query_id].update(query_set)
        symmetric_adjacency_dict[query_id].update(query_set)

        for subject_id in candidates:
            complementary_subject_id = helper.get_complementary_id(subject_id)

            subject_set = {subject_id,complementary_subject_id}
            symmetric_adjacency_dict[query_id].update(subject_set)
            symmetric_adjacency_dict[complementary_query_id].update(subject_set)
    return symmetric_adjacency_dict

class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self,x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self,x,y):
        root_x = self.find(x)
        root_y = self.find(y)
        if root_x != root_y:
            self.parent[root_x] = root_y

def get_disjoint_sets(setlist):
    uf = UnionFind()
    for s in setlist:
        for x in s:
            if x not in uf.parent:
                uf.parent[x] = x

    for reverse_id in [x for x in uf.parent if x.endswith("_Reversed")]:
        forward_id = reverse_id.replace("_Reversed","")
        if forward_id in uf.parent:
            uf.union(reverse_id,forward_id)

    for s in setlist:
        s_list = list(s)
        for i in range(1,len(s_list)):
            uf.union(s_list[0],s_list[i])

    groups = defaultdict(set)
    for x in uf.parent:
        groups[uf.find(x)].add(x)

    return list(groups.values())

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")
    
    """
    This is a list of pairs or singletons expected from the process_blast_results_module.
    The pairset of sequence IDs (i.e., {<id_i>,<id_j>}) is returned if their corrected BLAST
    score >= the specified homology threshold; otherwise, each sequence ID is appended as a
    singleton set (e.g., {<id_i>}, {id_j})
    """
    hitset_list = []
    for chunk_df in pd.read_csv(args.hits_path,sep="\t",chunksize=250_000,header=None,names=["qseqid","sseqid","score"]):
        for qseqid,sseqid,score in zip(chunk_df["qseqid"],chunk_df["sseqid"],chunk_df["score"]):
            if score >= args.threshold:
                hitset_list.append({qseqid,sseqid})
            else:
                hitset_list.append({qseqid})
                hitset_list.append({sseqid})

    homologous_groups = get_disjoint_sets(hitset_list)
    logger.info(f"{len(homologous_groups)} distinct groups.")

    with open(args.output_path,"w") as handle:
        for group_id,homologous_group in enumerate(homologous_groups):
            for sample_id in homologous_group:
                handle.write(f"{sample_id}\t{group_id}\n")

    logger.info(f"Homologous groups written to: {args.output_path}")
    logger.info(f"Module execution completed.\n")
