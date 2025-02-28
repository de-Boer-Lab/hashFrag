import os
import subprocess
import utils.helper_functions as helper
from glob import glob
import logging
from pathlib import Path

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")

    Path(args.output_dir).mkdir(parents=True,exist_ok=True)

    blast_module_files   = helper.blast_module_file_handler(args,logger)
    reference_label      = blast_module_files[0]
    reference_fasta_path = blast_module_files[1]
    query_label          = blast_module_files[2]
    query_fasta_path     = blast_module_files[3]
    blastdb_path         = blast_module_files[4]
    blastn_path          = blast_module_files[5]

    if glob(blastdb_path+"*") and not args.force:
        logger.info(f"Existing BLAST database found. Path: {blastdb_path}")
        logger.info(f"Skipping `makeblastdb` call...")
    else:
        if glob(blastdb_path+"*") and args.force:
            logger.warning(f"Overwriting previously computed BLAST database!")

        command = helper.construct_makeblastdb_command(
            fasta_path=reference_fasta_path,
            label=reference_label,
            blastdb_path=blastdb_path
        )
        process = subprocess.run(command,shell=True,capture_output=True)
        helper.validate_subprocess_execution(process,command,logger)
        logger.info(f"BLAST DataBase construction finished and written to: {blastdb_path}")
        
    if not args.exec_makeblastdb_only:

        if os.path.exists(blastn_path) and not args.force:
            logger.info(f"Existing BLAST results file found. Path: {blastn_path}")
            logger.info(f"Skipping `blastn` call.") 
        else:
            if os.path.exists(blastn_path) and args.force:
                logger.warning(f"Overwriting previously computed BLAST output file!") 

            command = helper.construct_blastn_command(
                query_path=query_fasta_path,
                blastdb_path=blastdb_path,
                blastn_path=blastn_path,
                word_size=args.word_size,
                gapopen=args.gapopen,
                gapextend=args.gapextend,
                penalty=args.penalty,
                reward=args.reward,
                max_target_seqs=args.max_target_seqs,
                xdrop_ungap=args.xdrop_ungap,
                xdrop_gap=args.xdrop_gap,
                xdrop_gap_final=args.xdrop_gap_final,
                evalue=args.evalue,
                dust=args.dust,
                threads=args.threads
            )
            process = subprocess.run(command,shell=True,capture_output=True)
            helper.validate_subprocess_execution(process,command,logger)
            logger.info(f"BLASTn process finished and written to: {blastn_path}")

    logger.info(f"Module execution completed.\n")
