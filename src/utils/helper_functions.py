import bz2
import gzip
import pickle
from Bio import SeqIO,pairwise2

def load_fasta_as_dictionary(path,idset=None):
    fasta_dict = {}
    if path.endswith(".gz"):
        handle = gzip.open(path,"rt")
    else:
        handle = open(path,"r")
    if idset is None:
        for record in SeqIO.parse(handle,"fasta"):
            fasta_dict[record.id] = str(record.seq)
    else:
        for record in SeqIO.parse(handle,"fasta"):
            if record.id in idset:
                fasta_dict[record.id] = str(record.seq)   
    handle.close()
    return fasta_dict

def load_fasta_as_list(path):
    ids  = []
    seqs = []
    with gzip.open(path,"rt") as handle:
        for record in SeqIO.parse(handle,"fasta"):
            ids.append(record.description)
            seqs.append(str(record.seq))
    return ids,seqs

def parse_fasta_to_dictionary(path):
    if path.endswith(".gz"):
        handle = gzip.open(path,"rt")
    else:
        handle = open(path,"rt")
    fasta_dict = {}
    for record in SeqIO.parse(handle,"fasta"):
        fasta_dict[record.id] = str(record.seq)
    handle.close()
    return fasta_dict

def load_collisions_dictionary(path):
    """ Loads collision dictionary from compressed pickle
    """
    with bz2.BZ2File(path,"rb") as handle:
        collision_dict = pickle.load(handle)
    return collision_dict

def write_collisions_dictionary(collision_dict,path):
    """ Writes collision dictionary to a compressed pickle
    """
    with bz2.BZ2File(path,"wb") as handle:
        pickle.dump(collision_dict,handle,protocol=pickle.HIGHEST_PROTOCOL)
    return

def parse_fasta_to_generator(path):
    """ Computing all-all SW scores for the provided sequence file
    """
    with gzip.open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield record.description, str(record.seq)

def pairwise_generator_from_fastas(fasta_path_i,fasta_path_j):
    """ Computing pairwise SW scores between the two sequence files provided.
    """
    fasta_dict_i = load_fasta_as_dictionary(fasta_path_i)
    fasta_dict_j = load_fasta_as_dictionary(fasta_path_j)
    for id_i,seq_i in fasta_dict_i.items():
        for id_j,seq_j in fasta_dict_j.items():
            yield ((id_i,seq_i),(id_j,seq_j))

def load_idset_from_pairwise_file(path):
    idset = set()
    with gzip.open(path,"rt") as handle:
        for line in handle:
            ids = line.strip().split(",")
            idset.update(ids)
    return idset

def pairwise_generator_from_pairwise_list(pairwise_path,fasta_path):
    """ Computing pairwise SW scores from a list of paired sequence IDs.
    """
    idset = load_idset_from_pairwise_file(pairwise_path)
    fasta_dict = load_fasta_as_dictionary(fasta_path,idset)
    with gzip.open(pairwise_path, "rt") as handle:
        for line in handle:
            id_i,id_j = line.strip().split(",")
            yield ((id_i,fasta_dict[id_i]),(id_j,fasta_dict[id_j]))

def compute_sw_score_mp(pair,match=1,mismatch=1,gap_open=2,gap_extend=0):
    """ 
    When only computing alignment score, Bio.pairwise2 
    faster than Bio.Align.PairwiseAligner.

    This input format is configured for multiprocessing.
    """
    (id_i,seq_i),(id_j,seq_j) = pair
    sw = pairwise2.align.localms(seq_i,seq_j,match,-mismatch,-gap_open,-gap_extend,score_only=True)
    return (id_i,id_j,sw)

def compute_sw_score(seq_i,seq_j,match=1,mismatch=1,gap_open=2,gap_extend=0):
    """ 
    When only computing alignment score, Bio.pairwise2 
    faster than Bio.Align.PairwiseAligner.

    This input format is configured for multiprocessing.
    """
    return pairwise2.align.localms(seq_i,seq_j,match,-mismatch,-gap_open,-gap_extend,score_only=True)

def get_pairwise_comparisons(collisions_dict):
    comparisons = set()
    for query_id,collisions in collisions_dict.items():
        for candidate_id in collisions:
            if query_id == candidate_id: continue
            if (query_id,candidate_id) in comparisons or (candidate_id,query_id) in comparisons:
                continue
            comparisons.add((query_id,candidate_id))
    return list(comparisons)

# TODO: might need a more general way to specify forward/reverse strands
def get_complementary_id(seq_id,idset):
    if seq_id.endswith("_Reversed"): # reverse strand
        seq_id_ = seq_id.replace("_Reversed","")
    else: # forward strand
        seq_id_ = seq_id+"_Reversed"
    if seq_id_ in idset:
        return seq_id_
    else:
        return

def assert_membership(seq_id,group_set):
    if seq_id is not None:
        assert seq_id in group_set
    return

def write_splits_to_file(train_split,test_split,path):
    with gzip.open(path,"wt") as handle:
        for seq_id in sorted(train_split):
            handle.write(f"{seq_id},train\n")
        for seq_id in sorted(test_split):
            handle.write(f"{seq_id},test\n")
    return