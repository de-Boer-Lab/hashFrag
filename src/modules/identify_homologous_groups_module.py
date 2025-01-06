import sys
import igraph as ig
import pandas as pd
from scipy.sparse import coo_matrix
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

def determine_homologous_communities(adjacency_dict):
    """
    From the adjacency dictionary representation of shared pairwise homology
    between sequences, we can first convert it into a sparse matrix (scipy dependency). 
    Next we construct a graph over the sparse matrix (igraph dependency), and then
    homologous groups can be simply identified by collecting the connected components
    in the graph.
    """

    candidate_idset = set()
    for candidates in adjacency_dict.values():
        candidate_idset.update(candidates)
    candidate_ids = list(candidate_idset)

    n = len(candidate_ids)
    id_map = {seq_id:i for i,seq_id in enumerate(candidate_ids)}

    row,col = [],[]
    for query_id,candidates in adjacency_dict.items():
        i = id_map[query_id]
        for subject_id in candidates:
            j = id_map[subject_id]
            row.append(i)
            col.append(j)
    weights = [1]*len(row)
    sparse_mat = coo_matrix((weights,(row,col)),shape=(n,n)).tocsr()

    g = ig.Graph.Weighted_Adjacency( # only Weighted_Adjancency seems to support sparse inputs
        matrix=sparse_mat,
        mode="undirected",
        loops=False,
        attr="candidates"
    )

    communities = g.connected_components(mode='weak')
    homologous_groups = []
    for community in communities:
        homologous_groups.append({ candidate_ids[i] for i in community })
    
    return homologous_groups

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")
    
    hits_dict = defaultdict(set)
    for chunk_df in pd.read_csv(args.hits_path,sep="\t",chunksize=250_000):
        for qseqid, sseqid in zip(chunk_df["id_i"],chunk_df["id_j"]):
            hits_dict[qseqid].add(sseqid)
            hits_dict[sseqid].add(qseqid)

    symmetric_hits_dict = get_symmetric_adjacency_dict(hits_dict)
    symmetric_hits_dict = {
        query_id:candidates
        for query_id,candidates in symmetric_hits_dict.items() 
        if candidates != {query_id,helper.get_complementary_id(query_id)}
    }
    if len(symmetric_hits_dict) == 0:
        logger.info("No homology detected!\n\tExiting.")
        sys.exit(0)
    else:
        logger.info(f"{len(symmetric_hits_dict)} sequences exhibiting homology.")

    homologous_groups = determine_homologous_communities(symmetric_hits_dict)
    logger.info(f"{len(homologous_groups)} distinct groups.")

    with open(args.output_path,"w") as handle:
        for group_id,homologous_group in enumerate(homologous_groups):
            for sample_id in homologous_group:
                handle.write(f"{sample_id}\t{group_id}\n")

    logger.info(f"Homologous groups written to: {args.output_path}")
    logger.info(f"Module execution completed.\n")
