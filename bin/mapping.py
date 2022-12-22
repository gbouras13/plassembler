import os
import sys
import subprocess as sp
import logging

#################################
# original mapping
#################################


def index_fasta(fasta,  logger):
    """ indexes fasta for bwa
	:param fasta: chromosome fasta
    :param logger: logger
    :return: 
    """
    try:
        bwa_index = sp.Popen(["bwa", "index", fasta ], stdout=sp.PIPE, stderr=sp.PIPE) 
        logging.write_to_log(bwa_index.stdout, logger)
    except:
        sys.exit("Error with bwa index\n")  

def minimap_long_reads(chromosome_flag, out_dir, threads, logger):
    """ maps long reads using minimap2
    :param chromosome_flag: True if mapping to chromosome
	:param out_dir: output directory path
    :param threads: threads
    :param logger: logger
    :return: 
    """
    input_long_reads = os.path.join(out_dir, "filtered_long_reads.fastq.gz")
    # chromosome is a flag for mapping chromosome or not
    if chromosome_flag == True:
        fasta = os.path.join(out_dir, "chromosome.fasta")
        sam = os.path.join(out_dir, "long_read_chromosome.sam")
    else:
        fasta = os.path.join(out_dir, "non_chromosome.fasta")
        sam = os.path.join(out_dir, "long_read_non_chromosome.sam")
    f = open(sam, "w")
    try:
        minimap = sp.Popen(["minimap2", "-ax", "map-ont", "-t", threads, fasta, input_long_reads ], stdout=f, stderr=sp.PIPE) 
        logging.write_to_log(minimap.stderr, logger)
    except:
        sys.exit("Error with minimap2\n")  

# bwa for short reads 

def bwa_map_short_reads(out_dir, chromosome_flag, threads, logger):
    """ maps short reads using bwa
    :param chromosome_flag: True if mapping to chromosome
	:param out_dir: output directory path
    :param threads: threads
    :param logger: logger
    :return: 
    """
    trim_one = os.path.join(out_dir, "trimmed_R1.fastq")
    trim_two = os.path.join(out_dir, "trimmed_R2.fastq")
    if chromosome_flag == False:
        fasta = os.path.join(out_dir, "non_chromosome.fasta")
        sam = os.path.join(out_dir, "short_read_non_chromosome.sam")
    else:
        fasta = os.path.join(out_dir, "chromosome.fasta")
        sam = os.path.join(out_dir, "short_read_chromosome.sam") 
    f = open(sam, "w")
    try:
        bwa_map = sp.Popen(["bwa", "mem", "-t", threads, fasta, trim_one, trim_two ], stdout=f, stderr=sp.PIPE) 
        logging.write_to_log(bwa_map.stderr, logger)
    except:
        sys.exit("Error with bwa mem\n")  