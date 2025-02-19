import os
import gzip
import logging
import logging.config
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

def instantiate_logger():
    log_config = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default": {
                "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S"
            },
            "detailed": {
                "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s - %(module)s - %(funcName)s",
                "datefmt": "%Y-%m-%d %H:%M:%S"
            }
        },
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "level": "INFO", 
                "formatter": "default",
                "stream": "ext://sys.stdout" 
            },
            "error": {
                "class": "logging.StreamHandler",
                "level": "ERROR",
                "formatter": "detailed",
                "stream": "ext://sys.stderr"
            }
        },
        "loggers": {
            "": {
                "level": "DEBUG",
                "handlers": ["console", "error"]
            }
        }
    }

    logging.config.dictConfig(log_config)

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

def generate_reverse_complement_fasta(input_path,output_path,logger,compressed=True,suffix="_Reversed"):
    warn = False

    inhandle = gzip.open(input_path,"rt") if input_path.endswith(".gz") else open(input_path,"r")
    outhandle = gzip.open(output_path,"wt") if compressed else open(output_path,"w")

    for record in SeqIO.parse(inhandle,"fasta"):
        if record.id.endswith(suffix):
            warn = True
        outhandle.write(f">{record.id}\n{str(record.seq)}\n")
        outhandle.write(f">{record.id+suffix}\n{str(record.seq.reverse_complement())}\n")

    inhandle.close()
    outhandle.close()

    if warn:
        logger.warning(f"Encountered sequence headers already containing reverse complement suffix ('{suffix}')! Make sure to specify '--skip-revcomp' if reverse complements are already present.")
    return warn

def load_fasta_ids(path):
    ids  = []
    with gzip.open(path,"rt") as handle:
        for record in SeqIO.parse(handle,"fasta"):
            ids.append(record.description)
    return ids

def get_complementary_id(seq_id):
    if seq_id.endswith("_Reversed"): # reverse strand
        return seq_id.replace("_Reversed","")
    else: # forward strand
        return seq_id+"_Reversed"

def write_splits_to_tsv(train_split,test_split,path):
    with gzip.open(path,"wt") as handle:
        handle.write("id\tsplit\n")
        for sample_id in train_split:
            handle.write(f"{sample_id}\ttrain\n")
        for sample_id in test_split:
            handle.write(f"{sample_id}\ttest\n")
    return

def valdidate_subprocess_execution(process,command,logger):
    if process.returncode != 0:
        logger.error(f"Error executing the following command: {command}")
        if process.stdout:
            logger.info(f"BLASTn output: {process.stdout.decode('utf-8')}")
        if process.stderr:
            logger.error(f"BLASTn error: {process.stderr.decode('utf-8')}")
        raise subprocess.CalledProcessError()

    if process.stdout:
        logger.info(f"BLASTn output: {process.stdout.decode('utf-8')}")
    
    if process.stderr:
        logger.error(f"BLASTn error: {process.stderr.decode('utf-8')}")

    return

def blast_file_validation(df,expected_columns):
    if df.shape[1] != len(expected_columns):
        raise Exception
    if df.isnull().all(axis=0).any():
        raise Exception

def is_valid_fasta_file(path):
    return (path is not None) and os.path.exists(path)
        
def open_fasta_file(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")