import bz2
import pickle
from datasketch import MinHash,MinHashLSHForest

def get_sequence_shingles(seq : str, k : int):
    """ 
    Decompose a sequence, `seq`, into its shingle set, where `k` refers to shingle (or k-mer) size.
    Note: the shingle set does not account for frequency with which k-mers appear.
    """
    shingles = set()
    for i in range(0,len(seq)-k+1,1):
        shingle = seq[i:i+k].encode('utf-8')
        shingles.add(shingle)
    return shingles

def MinHash_sequence_helper(seq : str, p : int, k : int, seed : int):
    """
    Given a sequence, `seq`, get its shingle set and subject each shingle to the MinHashing.
    This technique represents the input (shingle set or arbitrary size) in terms of a smaller, fixed-size
    output (referred to as the "signature vector").
    `p`: number of permutations (or hash functions) to use for MinHash. Higher p increases accuracy at the
    expense of increased time and memory cost (i.e., signature vector has a length of `p`).
    `k`: shingle (or k-mer) size.
    """
    minhash  = MinHash(num_perm=p,seed=seed)
    shingles = get_sequence_shingles(seq,k)
    for shingle in shingles:
        minhash.update(shingle)
    return minhash

def MinHash_sequences(fasta_dict,p,k,seed):
    """
    A wrapper to MinHash a list of sequences - `seqs` - with the specified number of permutations `p`
    and shingle (or k-mer size) `k`. The corresponding list of sequence IDs are provided (`ids`).
    This outputs a dictionary where keys correspond to sequence IDs and values refer to the
    MinHash object (signature vector) created by the MinHashing technique.
    """
    minhash_dict = {}
    for seq_id,seq in fasta_dict.items():
        minhash_dict[seq_id] = MinHash_sequence_helper(seq,p,k,seed)
    return minhash_dict

def LSH_forest_wrapper(minhash_dict,p,l):
    """
    Given the dictionary of MinHash objects corresponding to the input sequences, create LSH Forest
    object. 
    `p`: number of permutations in the LSH method
    `l`: number of LSH prefix Trees in the LSH Forest. Note that this cannot exceed the value of `p`.
    """
    forest = MinHashLSHForest(num_perm=p,l=l)
    for seq_id,minhash in minhash_dict.items():
        forest.add(seq_id,minhash)
    forest.index()
    return forest

def approximate_nearest_neighbors_queries(lsh_forest,minhash_dict,n):
    """
    Given the LSH Forest object, collect results into a dictionary where there is a key for
    each sequence ID (query) eand the corresponding sequence ID is collides with during
    the LSH procedure. These collisions denote 'candidate pairings' that can be subsequently
    be validated directly by computing the actual SW scores
    """
    collision_dict = {}
    for seq_id,minhash in minhash_dict.items():
        collision_dict[seq_id] = lsh_forest.query(minhash=minhash,k=n)
    return collision_dict
