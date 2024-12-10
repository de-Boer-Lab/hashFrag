import os
import sys
import gzip
import subprocess
import utils.helper_functions as helper

def run(args):

    label = os.path.basename(args.fasta_path).replace(".fa.gz","")
    blastdb_path = os.path.join(args.out_dir,f"{label}.blastdb")
    blastn_path = os.path.join(args.out_dir,f"{label}.blastn.out")

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
        helper.valdidate_subprocess_execution(process,command)

    print(f"BLAST DataBase construction finished and written to: {blastdb_path}")
    
    with gzip.open(args.fasta_path,"rb") as fasta:
        command = " ".join([
            "blastn",
            f"-db {blastdb_path}",
            f"-out {blastn_path}",
            f"-word_size {args.word_size}",
            f"-gapopen {args.gapopen}",
            f"-gapextend {args.gapextend}",
            f"-penalty {args.penalty}",
            f"-reward {args.reward}",
            f"-max_target_seqs {args.max_target_seqs}",
            f"-evalue {args.e_value}",
            f"-dust {args.dust}",
            f"-num_threads {args.threads}",
            "-outfmt '6 std score positive gaps'",
            "-strand plus",
            " "
        ])
        if args.blastn_args:
            command += args.blastn_args
        process = subprocess.run(command,input=fasta.read(),shell=True,capture_output=True)
        helper.valdidate_subprocess_execution(process,command)
    
    print(f"BLASTn process finished and written to: {blastn_path}")
