from collections import defaultdict
import utils.helper_functions as helper

def run(args):

    collisions_dict = helper.load_collisions_dictionary(args.collisions_path)
    comparisons = helper.get_pairwise_comparisons(collisions_dict)
    print(f"Number of pairwise comparisons to compute: {len(comparisons)}",flush=True)

    fasta_dict = helper.load_fasta_as_dictionary(args.fasta_path)

    score_dict = {}
    for (id_i,id_j) in comparisons:
        score = helper.compute_sw_score(fasta_dict[id_i],fasta_dict[id_j])
        score_dict[(id_i,id_j)] = score
        score_dict[(id_j,id_i)] = score

    filtered_collisions_dict = defaultdict(list)
    for query_id,collisions in collisions_dict.items():
        for candidate_id in collisions:
            if query_id == candidate_id: continue
            if score_dict[(query_id,candidate_id)] >= args.threshold:
                filtered_collisions_dict[query_id].append(candidate_id)
    print(f"{len(filtered_collisions_dict)} sequences remain after filtering according to a threshold of {args.threshold}",flush=True)
    helper.write_collisions_dictionary(filtered_collisions_dict,args.output_path)
    print(f"Filtered collision results written to: {args.output_path}")
