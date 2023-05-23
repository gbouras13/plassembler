import os
import sys
import subprocess as sp
import pandas as pd
import log
import logging



def run_flye(out_dir, threads, raw_flag, pacbio_model, logger):
    """Runs flye on trimmed long reads

    :param out_dir: output directory
    :param raw_flag: boolean - true if --nano-raw used for flue
    :param logger: logger
    :return:
    """
    trim_long = os.path.join(out_dir, "chopper_long_reads.fastq.gz")
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

def run_raven(out_dir, threads, logger):
    """Runs raven on trimmed long reads

    :param out_dir: output directory
    :param raw_flag: boolean - true if --nano-raw used for flye
    :param logger: logger
    :return:
    """
    trim_long = os.path.join(out_dir, "chopper_long_reads.fastq.gz")
    # gfa
    gfa = os.path.join(out_dir, "assembly_graph.gfa")

    f = open(os.path.join(out_dir, "assembly.fasta"), "w")

    # default to 8 threads
    command = "raven -t " + str(threads) + " " + trim_long +  " --graphical-fragment-assembly " + gfa

    try:
        # Create a subprocess using Popen
        process = sp.Popen(command, shell=True, stdout=f, stderr=sp.PIPE)
        # Read the output and error from the subprocess
        output, raven_log = process.communicate()
        logger.log(logging.INFO, raven_log)
    
    except:
        sys.exit("Error with Raven\n")  



