
import os
import sys
import subprocess as sp
import log


def deduplicate_fastqs(out_dir, threads, logger):
    """ deduplicates fastq
    :param out_dir:  Output Directory
    :param logger: logger
    :param threads: threads
    :return: 
    """
    # sets file paths
    concat_long = os.path.join(out_dir, "long_read_concat.fastq")
    dedup_long = os.path.join(out_dir, "long_read_dedup.fastq")
    concat_short_one = os.path.join(out_dir, "short_read_concat_R1.fastq.gz")
    dedup_short_one = os.path.join(out_dir, "short_read_dedup_R1.fastq.gz")
    concat_short_two = os.path.join(out_dir, "short_read_concat_R2.fastq.gz")
    dedup_short_two = os.path.join(out_dir, "short_read_dedup_R2.fastq.gz")
    try:
        dedup_long= sp.Popen(["seqkit", "rmdup", concat_long, "-j", threads, "-n", "-o",dedup_long], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(dedup_long.stdout, logger)
        dedup_s1= sp.Popen(["seqkit", "rmdup", concat_short_one, "-j", threads, "-n", "-o",dedup_short_one], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(dedup_s1.stdout, logger)
        dedup_s2= sp.Popen(["seqkit", "rmdup", concat_short_two, "-j", threads, "-n", "-o",dedup_short_two], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(dedup_s2.stdout, logger)
    except:
        sys.exit("Error with seqkit\n")  

def deduplicate_long_fastqs(out_dir, threads, logger):
    """ deduplicates fastq
    :param out_dir:  Output Directory
    :param logger: logger
    :param threads: threads
    :return: 
    """
    # sets file paths
    concat_long = os.path.join(out_dir, "long_read_concat.fastq")
    dedup_long = os.path.join(out_dir, "long_read_dedup.fastq")
    try:
        dedup_long= sp.Popen(["seqkit", "rmdup", concat_long, "-j", threads, "-n", "-o",dedup_long], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(dedup_long.stdout, logger)
    except:
        sys.exit("Error with seqkit\n")  