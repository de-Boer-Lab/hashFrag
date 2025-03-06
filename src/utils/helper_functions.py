import os
import gzip
import logging
import logging.config
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

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
    if idset is None:
        for record in SeqIO.parse(path,"fasta"):
            fasta_dict[record.id] = str(record.seq)
    else:
        for record in SeqIO.parse(path,"fasta"):
            if record.id in idset:
                fasta_dict[record.id] = str(record.seq)   
    return fasta_dict

def generate_reverse_complement_fasta(input_path,output_path,logger,suffix="_Reversed"):
    warn = False

    with open(output_path,"w") as outhandle:
        for record in SeqIO.parse(input_path,"fasta"):
            if record.id.endswith(suffix):
                warn = True
            outhandle.write(f">{record.id}\n{str(record.seq)}\n")
            outhandle.write(f">{record.id+suffix}\n{str(record.seq.reverse_complement())}\n")

    if warn:
        logger.warning(f"Encountered sequence headers already containing reverse complement suffix ('{suffix}')! Make sure to specify '--skip-revcomp' if reverse complements are already present.")
    return warn

def load_fasta_ids(path):
    ids  = []
    for record in SeqIO.parse(path,"fasta"):
        ids.append(record.id)
    return ids

def get_complementary_id(seq_id):
    if seq_id.endswith("_Reversed"): # reverse strand
        return seq_id.replace("_Reversed","")
    else: # forward strand
        return seq_id+"_Reversed"

def write_splits_to_tsv(train_split,test_split,path):
    with open(path,"w") as handle:
        handle.write("id\tsplit\n")
        for sample_id in train_split:
            handle.write(f"{sample_id}\ttrain\n")
        for sample_id in test_split:
            handle.write(f"{sample_id}\ttest\n")
    return


def validate_subprocess_execution(process,command,logger):
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

def partition_fasta_by_seq_count(fasta_path,label,partition_dir,partition_size):

    p = 1
    pathlist = []
    records = []
    for record in SeqIO.parse(fasta_path,"fasta"):
        records.append(record)
        if len(records) >= partition_size:
            path = os.path.join(partition_dir,f"{label}.part_{p:0>8}.fa")
            SeqIO.write(records,path,"fasta")
            pathlist.append(path)
            records.clear()
            p += 1

    if records:
        path = os.path.join(partition_dir,f"{label}.part_{p:0>8}.fa")
        SeqIO.write(records,path,"fasta")
        pathlist.append(path)

    return pathlist

def write_filepaths_to_txt(pathlist,output_file):
    with open(output_file,"w") as handle:
        handle.write("\n".join(pathlist)+"\n")
    return

def extract_fasta_file_label(fasta_path):
    label = os.path.basename(fasta_path)
    label = label.replace(".fa","").replace(".fasta","").replace(".fna","")
    return label

def construct_makeblastdb_command(fasta_path,label,blastdb_path):
    command = " ".join([
        "makeblastdb",
        f"-in {fasta_path}",
        "-dbtype nucl",
        "-input_type 'fasta'",
        f"-title {label}",
        f"-out {blastdb_path}"
    ])
    return command

def construct_blastn_command(query_path,blastdb_path,blastn_path,
                             word_size,gapopen,gapextend,penalty,reward,
                             max_target_seqs,xdrop_ungap,xdrop_gap,xdrop_gap_final,
                             evalue,dust,threads):
    command = " ".join([
        "blastn",
        f"-query {query_path}",
        f"-db {blastdb_path}",
        f"-out {blastn_path}",
        f"-word_size {word_size}",
        f"-gapopen {gapopen}",
        f"-gapextend {gapextend}",
        f"-penalty {penalty}",
        f"-reward {reward}",
        f"-max_target_seqs {max_target_seqs}",
        f"-xdrop_ungap {xdrop_ungap}",
        f"-xdrop_gap {xdrop_gap}",
        f"-xdrop_gap_final {xdrop_gap_final}",
        f"-evalue {evalue}",
        f"-dust {dust}",
        f"-num_threads {threads}",
        "-outfmt '6 std score positive gaps'",
        "-strand plus"
    ])
    return command

def construct_process_blast_results_command(blastn_path,processed_blastn_path):
    command = " ".join([
        "hashFrag process_blast_results_module",
        f"--blastn-path {blastn_path}",
        f"--processed-blastn-path {processed_blastn_path}"
    ])
    return command

def instantiate_job_script(queries_path,n_queries,blastdb_path,blast_dir,
                           job_scheduler,job_account,job_name,job_memory,num_cpus,job_time,
                           stdout_dir,stderr_dir,environment_path):
    if job_scheduler == "slurm":
        job_script = [
            "#!/bin/bash",
            f"#SBATCH --account={job_account}",
            f"#SBATCH --job-name={job_name}",
            "#SBATCH --ntasks=1",
            "#SBATCH --nodes=1",
            f"#SBATCH --mem={job_memory}",
            f"#SBATCH --cpus-per-task={num_cpus}",
            f"#SBATCH --time={job_time}",
            f"#SBATCH --array=1-{n_queries}",
            f"#SBATCH --output={os.path.join(stdout_dir,f'{job_name}.jobid_%A_%a.out')}",
            f"#SBATCH --error={os.path.join(stderr_dir,f'{job_name}.jobid_%A_%a.err')}\n",
            f"source {environment_path}\n",
            f"BLASTDB_PATH={blastdb_path}",
            f"QUERIES_PATH={queries_path}",
            f"BLAST_DIR={blast_dir}",
            "QUERY_PATH=$( sed -n ${SLURM_ARRAY_TASK_ID}p $QUERIES_PATH )",
            "LABEL=$( basename -s '.fa' $QUERY_PATH )",
            "BLASTN_PATH=$BLAST_DIR/${LABEL}.blastn.out",
            "PROCESSED_BLASTN_PATH=$BLAST_DIR/${LABEL}.blastn.processed.tsv\n",
            'echo "hashFrag blastn_array_module"',
            'echo "BLAST database path: $BLASTDB_PATH"',
            'echo "Query path: $QUERY_PATH"',
            'echo "BLAST results path: $BLASTN_PATH"',
            'echo "Processed BLAST results path: $PROCESSED_BLASTN_PATH"\n'
        ]
    elif job_scheduler == "sge":
        job_script = [
            "#!/bin/bash",
            f"#$ -N {job_name}",
            f"#$ -l s_vmem={job_memory}",
            f"#$ -pe def_slot {num_cpus}",
            f"#$ -t 1:{n_queries}",
            f"#$ -o {stdout_dir}",
            f"#$ -e {stderr_dir}\n",
            f"source {environment_path}\n",
            f"BLASTDB_PATH={blastdb_path}",
            f"QUERIES_PATH={queries_path}",
            f"BLAST_DIR={blast_dir}",
            "QUERY_PATH=$( sed -n ${SGE_TASK_ID}p $QUERIES_PATH )",
            "LABEL=$( basename -s '.fa' $QUERY_PATH )",
            "BLASTN_PATH=$BLAST_DIR/${LABEL}.blastn.out",
            "PROCESSED_BLASTN_PATH=$BLAST_DIR/${LABEL}.blastn.processed.tsv\n",
            'echo "hashFrag blastn_array_module"',
            'echo "BLAST database path: $BLASTDB_PATH"',
            'echo "Query path: $QUERY_PATH"',
            'echo "BLAST results path: $BLASTN_PATH"',
            'echo "Processed BLAST results path: $PROCESSED_BLASTN_PATH"\n',
        ]
    else:
        raise Exception(f"Unrecognized `job_scheduler` specified ({job_scheduler} is invalid or currently not supported)!")

    return job_script

def blast_module_file_handler(args,logger):

    output_dir = Path(args.output_dir).resolve()

    # Existing train AND test splits provided as FASTA files 
    if is_valid_fasta_file(args.train_fasta_path) and is_valid_fasta_file(args.test_fasta_path):
        logger.info("Train and test FASTA files detected. Computing pairwise BLAST comparisons across splits...")

        train_fasta_path = Path(args.train_fasta_path).resolve()
        test_fasta_path = Path(args.test_fasta_path).resolve()

        if args.train_fasta_path.endswith(".gz"):
            logger.error(f"Error with the following input: {train_fasta_path}")
            raise Exception(f"FASTA file must be uncompressed!")

        if args.test_fasta_path.endswith(".gz"):
            logger.error(f"Error with the following input: {test_fasta_path}")
            raise Exception(f"FASTA file must be uncompressed!")
    
        if args.blastdb_label is None:
            ref_label    = extract_fasta_file_label(train_fasta_path)
            query_label  = extract_fasta_file_label(test_fasta_path)
            blastdb_path = os.path.join(output_dir,f"{ref_label}.blastdb")
            blastn_path  = os.path.join(output_dir,f"{query_label}.blastn.out")
        else:
            # Useful when running multiple BLASTs against the same data base
            ref_label    = args.blastdb_label
            query_label  = extract_fasta_file_label(test_fasta_path)
            blastdb_path = os.path.join(output_dir,f"{ref_label}.blastdb")
            blastn_path  = os.path.join(output_dir,f"{query_label}.blastn.out")

        if args.skip_revcomp:
            ref_fasta_path = train_fasta_path # BLAST database
        else:
            ref_fasta_path = os.path.join(output_dir,f"{ref_label}.revcomps.fa")
            if not args.force and os.path.exists(ref_fasta_path):
                logger.info("File with generated reverse complement ref sequences already exists (skipping)...")
            else:
                logger.info("Generating reverse complement ref sequences...")
                generate_reverse_complement_fasta(train_fasta_path,ref_fasta_path,logger)

        query_fasta_path = test_fasta_path # Query each of these sequences against the BLAST database

    # All sequences in dataset in a single file
    elif is_valid_fasta_file(args.fasta_path):
        logger.info("One FASTA files detected. Computing pairwise BLAST comparisons for all sequence-pairs...")

        fasta_path = Path(args.fasta_path).resolve()

        if args.fasta_path.endswith(".gz"):
            logger.error(f"Error with the following input: {fasta_path}")
            raise Exception(f"FASTA file must be uncompressed!")

        if args.blastdb_label is None:
            ref_label    = extract_fasta_file_label(fasta_path)
            query_label  = ref_label
            blastdb_path = os.path.join(output_dir,f"{ref_label}.blastdb")
            blastn_path  = os.path.join(output_dir,f"{query_label}.blastn.out")
        else:
            # Useful when running multiple BLASTs against the same data base
            ref_label    = args.blastdb_label
            query_label  = extract_fasta_file_label(fasta_path)
            blastdb_path = os.path.join(output_dir,f"{args.blastdb_label}.blastdb")
            blastn_path  = os.path.join(output_dir,f"{args.blastdb_label}.blastn.out")

        if args.skip_revcomp:
            ref_fasta_path = fasta_path # BLAST database
            query_fasta_path = ref_fasta_path
        else:
            ref_fasta_path = os.path.join(output_dir,f"{ref_label}.revcomps.fa")
            if os.path.exists(ref_fasta_path) and not args.force:
                logger.info("File with generated reverse complement ref sequences already exists (skipping)...")
            else:
                logger.info("Generating reverse complement ref sequences...")
                generate_reverse_complement_fasta(fasta_path,ref_fasta_path,logger)
            query_fasta_path = fasta_path

    else:
 
        raise FileNotFoundError("Error: no input FASTA file(s) detected.")

    return (ref_label,ref_fasta_path,query_label,query_fasta_path,blastdb_path,blastn_path)


