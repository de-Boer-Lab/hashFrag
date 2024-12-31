import os
import subprocess
import utils.helper_functions as helper
from glob import glob
import logging

from os.path import basename,join

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")
    if helper.is_valid_fasta_file(args.train_fasta_path) and helper.is_valid_fasta_file(args.test_fasta_path):
        logger.info("Train and test FASTA files detected. Computing pairwise BLAST comparisons across splits...")
        if args.blastdb_label is None:
            train_label  = basename(args.train_fasta_path).replace(".fa.gz","").replace(".fa","")
            test_label   = basename(args.test_fasta_path).replace(".fa.gz","").replace(".fa","")
            label        = train_label
            blastdb_path = join(args.output_dir,f"{train_label}.blastdb")
            blastn_path  = join(args.output_dir,f"{test_label}.blastn.out")
        else:
            label        = args.blastdb_label
            blastdb_path = join(args.output_dir,f"{args.blastdb_label}.blastdb")
            blastn_path  = join(args.output_dir,f"{args.blastdb_label}.blastn.out")

        reference_fasta_path = args.train_fasta_path
        query_fasta_path     = args.test_fasta_path

    elif helper.is_valid_fasta_file(args.fasta_path):
        logger.info("One FASTA files detected. Computing pairwise BLAST comparisons for all sequence-pairs...")
        if args.blastdb_label is None:
            label        = basename(args.fasta_path).replace(".fa.gz","").replace(".fa","")
            blastdb_path = join(args.output_dir,f"{label}.blastdb")
            blastn_path  = join(args.output_dir,f"{label}.blastn.out")
        else:
            label        = args.blastdb_label
            blastdb_path = join(args.output_dir,f"{args.blastdb_label}.blastdb")
            blastn_path  = join(args.output_dir,f"{args.blastdb_label}.blastn.out")

        reference_fasta_path = args.fasta_path
        query_fasta_path     = args.fasta_path

    else:

        raise FileNotFoundError("Error: no input FASTA file(s) detected.")

    if glob(blastdb_path+".*"):
        logger.info(f"Existing BLAST database found. Path: {blastdb_path}\n\tskipping `makeblastdb` call.")
    else:
        with helper.open_fasta_file(reference_fasta_path) as fasta:
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
            helper.valdidate_subprocess_execution(process,command,logger)
            logger.info(f"BLAST DataBase construction finished and written to: {blastdb_path}")
 
    if os.path.exists(blastn_path):
        logger.info(f"Existing BLAST results file found. Path: {blastn_path}\n\tskipping `blastn` call.")
    else:
        with helper.open_fasta_file(query_fasta_path) as fasta:
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
            helper.valdidate_subprocess_execution(process,command,logger)
    
        logger.info(f"BLASTn process finished and written to: {blastn_path}")
    
    logger.info(f"Module execution completed.\n")
