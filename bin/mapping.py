import os
import sys
import subprocess as sp
import log

#################################
# original mapping
#################################


def minimap_long_reads(out_dir, threads, logger):
    """ maps long reads using minimap2
    :param chromosome_flag: True if mapping to chromosome
	:param out_dir: output directory path
    :param threads: threads
    :param logger: logger
    :return: 
    """
    input_long_reads = os.path.join(out_dir, "final_filtered_long_reads.fastq.gz")

    fasta = os.path.join(out_dir, "flye_renamed.fasta")
    sam = os.path.join(out_dir, "long_read.sam")

    f = open(sam, "w")
    try:
        minimap = sp.Popen(["minimap2", "-ax", "map-ont", "-t", threads, fasta, input_long_reads ], stdout=f, stderr=sp.PIPE) 
        log.write_to_log(minimap.stderr, logger)
    except:
        sys.exit("Error with minimap2\n")  

# bwa for short reads 

def minimap_short_reads(out_dir,  threads, logger):
    """ maps short reads using minimap
	:param out_dir: output directory path
    :param threads: threads
    :param logger: logger
    :return: 
    """
    trim_one = os.path.join(out_dir, "trimmed_R1.fastq")
    trim_two = os.path.join(out_dir, "trimmed_R2.fastq")
    fasta = os.path.join(out_dir, "flye_renamed.fasta")
    sam = os.path.join(out_dir, "short_read.sam")

    f = open(sam, "w")
    try:
        minimap2_map = sp.Popen(["minimap2", "-ax", "sr", "-t", threads, fasta, trim_one, trim_two ], stdout=f, stderr=sp.PIPE) 
        log.write_to_log(minimap2_map.stderr, logger)
    except:
        sys.exit("Error with bwa mem\n")  