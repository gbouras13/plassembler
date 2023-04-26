import os
import sys
import subprocess as sp
import pandas as pd
import log
from plass_class import Plass 



def run_flye(out_dir, threads, raw_flag, logger):
    """Runs flye on trimmed long reads

    :param out_dir: output directory
    :param threads: threads
    :param raw_flag: boolean - true if --nano-raw used for flue
    :param logger: logger
    :return:
    """
    trim_long = os.path.join(out_dir, "final_filtered_long_reads.fastq.gz")
    nanopore_flye_model = "--nano-hq"
    if raw_flag == True:
        nanopore_flye_model = "--nano-raw"
    try:
        flye = sp.Popen(["flye", nanopore_flye_model, trim_long, "--out-dir", out_dir, "--threads", threads], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(flye.stdout, logger)
    except:
        sys.exit("Error with Flye\n")  



