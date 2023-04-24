
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
    concat_short_one = os.path.join(out_dir, "short_read_concat_chrom_non_chrom_R1.fastq")
    dedup_short_one = os.path.join(out_dir, "multimap_plasmid_chromosome_R1.fastq")
    rubbish_short_one = os.path.join(out_dir, "rubbish_R1.fastq")
    concat_short_two = os.path.join(out_dir, "short_read_concat_chrom_non_chrom_R2.fastq")
    dedup_short_two = os.path.join(out_dir, "multimap_plasmid_chromosome_R2.fastq")
    rubbish_short_two = os.path.join(out_dir, "rubbish_R2.fastq")
    try:
        dedup_s1= sp.Popen(["seqkit", "rmdup", concat_short_one, "-j", threads, "-n", "-d",dedup_short_one, "-o", rubbish_short_one], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(dedup_s1.stdout, logger)
        dedup_s2= sp.Popen(["seqkit", "rmdup", concat_short_two, "-j", threads, "-n", "-d",dedup_short_two, "-o", rubbish_short_two], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(dedup_s2.stdout, logger)
    except:
        sys.exit("Error with seqkit\n")  

