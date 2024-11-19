import os
import gzip
import subprocess
from collections import defaultdict
import utils.helper_functions as helper

def valdidate_subprocess_execution(process,command):
    if process.returncode != 0:
        print(f"Error executing the following command: {command}")
        print(process.stderr.decode("utf-8"))
        raise subprocess.CalledProcessError()
    print(process.stdout.decode("utf-8"))
    return

def run(args):

    label = os.path.basename(args.fasta_path).replace(".fa.gz","")
    out_dir = os.path.dirname(args.output_path)
    blastdb_path = os.path.join(out_dir,f"{label}.blastdb")
    blastn_path = os.path.join(out_dir,f"{label}.blastn.out")

    with gzip.open(args.fasta_path,"rb") as fasta:
        command = " ".join([
            "makeblastdb",
            "-dbtype nucl",
            "-input_type 'fasta'",
            f"-title {label}",
            f"-out {blastdb_path}",
            " "
        ])
        if args.blastdb_args:
            command += args.blastdb_args
        process = subprocess.run(command,input=fasta.read(),shell=True,capture_output=True)
        valdidate_subprocess_execution(process,command)

    with gzip.open(args.fasta_path,"rb") as fasta:
        command = " ".join([
            "blastn",
            "-task blastn",
            f"-db {blastdb_path}",
            f"-out {blastn_path}",
            f"-word_size {args.word_size}",
            "-outfmt 6",
            " "
        ])
        if args.blastn_args:
            command += args.blastn_args
        process = subprocess.run(command,input=fasta.read(),shell=True,capture_output=True)
        valdidate_subprocess_execution(process,command)
        
    collisions_dict = defaultdict(list)
    with open(blastn_path,"r") as handle:
        for line in handle:
            entry = line.strip().split("\t")
            collisions_dict[entry[0]].append(entry[1])

    helper.write_collisions_dictionary(collisions_dict,args.output_path)
    print(f"BLASTn process finished and written to: {args.output_path}")
