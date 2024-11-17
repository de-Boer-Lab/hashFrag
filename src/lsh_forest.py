import utils.helper_functions as helper
import utils.lsh_forest_helpers as hasher

def run(args):

    k = args.kmer_size
    p = args.permutations
    l = args.n_trees
    n = args.neighbours

    fasta_dict     = helper.load_fasta_as_dictionary(args.fasta_path)
    minhash_dict   = hasher.MinHash_sequences(fasta_dict,p,k,args.seed)
    lsh_object     = hasher.LSH_forest_wrapper(minhash_dict,p,l)
    collision_dict = hasher.approximate_nearest_neighbors_queries(lsh_object,minhash_dict,n)
    helper.write_collisions_dictionary(collision_dict,args.output_path)
    print(f"LSH Forest process finished and written to: {args.output_path}")
