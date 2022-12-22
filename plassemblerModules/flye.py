import os
import sys
import subprocess as sp
import pandas as pd
import plassemblerModules



def run_flye(out_dir, threads, raw_flag, logger):
    """Runs flye on trimmed long reads

    :param out_dir: output directory
    :param threads: threads
    :param raw_flag: boolean - true if --nano-raw used for flue
    :param logger: logger
    :return:
    """
    trim_long = os.path.join(out_dir, "filtered_long_reads.fastq.gz")
    nanopore_flye_model = "--nano-hq"
    if raw_flag == True:
        nanopore_flye_model = "--nano-raw"
    try:
        flye = sp.Popen(["flye", nanopore_flye_model, trim_long, "--out-dir", out_dir, "--threads", threads], stdout=sp.PIPE, stderr=sp.PIPE) 
        plassemblerModules.write_to_log(flye.stdout, logger)
    except:
        sys.exit("Error with Flye\n")  

def contig_count(out_dir):
    """ Counts the number of contigs assembled by flye
    :param out_dir: output directory
    :param logger: logger
    :return:
    """
    info_file =  os.path.join(out_dir, "assembly_info.txt")
    col_list = ["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"] 
    info_df = pd.read_csv(info_file, delimiter= '\t', index_col=False , names=col_list, skiprows=1) 
    contig_count = len(info_df['seq_name'])
    print("Flye assembled " + str(contig_count) + " contigs.")
    return contig_count

