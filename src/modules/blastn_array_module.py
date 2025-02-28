import os
import logging
import subprocess
import utils.helper_functions as helper
from glob import glob
from pathlib import Path

def run(args):

    logger = logging.getLogger(__name__.replace("modules.",""))
    logger.info("Calling module...")

    Path(args.output_dir).mkdir(parents=True,exist_ok=True)
    partition_dir = os.path.join(args.output_dir,"query_fastas")
    blast_dir     = os.path.join(args.output_dir,"blast_results")
    stdout_dir    = os.path.join(args.output_dir,"job_stdout_files")
    stderr_dir    = os.path.join(args.output_dir,"job_stderr_files")
    Path(partition_dir).mkdir(exist_ok=True,parents=False)
    Path(blast_dir).mkdir(exist_ok=True,parents=False)
    Path(stdout_dir).mkdir(exist_ok=True,parents=False)
    Path(stderr_dir).mkdir(exist_ok=True,parents=False)

    blast_module_files   = helper.blast_module_file_handler(args,logger)
    reference_label      = blast_module_files[0]
    reference_fasta_path = blast_module_files[1]
    query_label          = blast_module_files[2]
    query_fasta_path     = blast_module_files[3]
    blastdb_path         = blast_module_files[4]

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

    logger.info("Partitioning query FASTA file...")
    query_fasta_paths = helper.partition_fasta_by_seq_count(query_fasta_path,query_label,partition_dir,args.query_partition_size)
    queries_path = os.path.join(args.output_dir,"query_fasta_paths.txt")
    helper.write_filepaths_to_txt(query_fasta_paths,queries_path)

    """
    This initial instantiation aims to handle all of the arguments and setup specific
    to the specified job scheduler. It will setup a basic array job where each
    task calls blastn of a (partitioned) query FASTA file to the complete BLAST database.
    """
    logger.info("Constructing array job script...")
    job_script = helper.instantiate_job_script(
        queries_path=queries_path,
        n_queries=len(query_fasta_paths),
        blastdb_path=blastdb_path,
        blast_dir=blast_dir,
        job_scheduler=args.job_scheduler,
        job_account=args.job_account,
        job_name=args.job_name,
        job_memory=args.job_memory,
        num_cpus=args.num_cpus,
        job_time=args.job_time,
        stdout_dir=stdout_dir,
        stderr_dir=stderr_dir,
        environment_path=args.environment_path
    )

    blastn_command = helper.construct_blastn_command(
        query_path="$QUERY_PATH",
        blastdb_path="$BLASTDB_PATH",
        blastn_path="$BLASTN_PATH",
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
        threads=args.num_cpus
    )
    job_script.append(f"echo 'BLASTn command:\n{blastn_command}'\n")
    job_script.append(blastn_command)

    processing_command = helper.construct_process_blast_results_command(
        blastn_path="$BLASTN_PATH",
        processed_blastn_path="$PROCESSED_BLASTN_PATH"
    )
    job_script.append(f"echo 'hashFrag process_blast_results_module command:\n{processing_command}\n'")
    job_script.append(processing_command)

    script_path = os.path.join(args.output_dir,"hashFrag.blastn_array_module.array_jobs.sh")
    with open(script_path,"w") as handle:
        handle.write("\n".join(job_script)+"\n")

    logger.info(f"Job script written to: {script_path}")
    logger.info(f"Module execution completed.\n")
