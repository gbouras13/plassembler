import os
import sys
import subprocess as sp
import pandas as pd
import log
from plass_class import Plass 



def run_flye(out_dir, threads, raw_flag, pacbio_model, logger):
    """Runs flye on trimmed long reads

    :param out_dir: output directory
    :param threads: threads
    :param raw_flag: boolean - true if --nano-raw used for flue
    :param logger: logger
    :return:
    """
    trim_long = os.path.join(out_dir, "final_filtered_long_reads.fastq.gz")
    flye_model = "--nano-hq"
    if raw_flag == True:
        flye_model = "--nano-raw"
    if pacbio_model != 'nothing':
        flye_model = pacbio_model
    try:
        flye = sp.Popen(["flye", flye_model, trim_long, "--out-dir", out_dir, "--threads", threads], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(flye.stdout, logger)
    except:
        sys.exit("Error with Flye\n")  



