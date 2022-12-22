import os
import sys
import subprocess as sp
import log

# sam to bam
def sam_to_bam(out_dir, combo, threads, logger):
    """ converts sam to bam using pysam
    :param combo: combination of chromosome and read: chromosome_long, non_chromosome_long, non_chromosome_short or chromosome_short
	:param out_dir: output directory path
    :param threads: threads
    :param logger: logger
    :return: 
    """
    #  combo for what combination to index
    if combo == "chromosome_long":
        sam = os.path.join(out_dir, "long_read_chromosome.sam")
        bam = os.path.join(out_dir, "long_read_chromosome.bam")
    elif combo == "non_chromosome_long":
        sam = os.path.join(out_dir, "long_read_non_chromosome.sam")
        bam = os.path.join(out_dir, "long_read_non_chromosome.bam")
    elif combo == "non_chromosome_short":
        sam = os.path.join(out_dir, "short_read_non_chromosome.sam")
        bam = os.path.join(out_dir, "short_read_non_chromosome.bam")
    else: # chromosome_short
        sam = os.path.join(out_dir, "short_read_chromosome.sam")
        bam = os.path.join(out_dir, "short_read_chromosome.bam")
    outFile = open(bam, "w")
    try:
        sam_to_bam = sp.Popen(["samtools", "view", "-h", "-@", threads, "-b", sam], stdout=outFile, stderr=sp.PIPE) 
        log.write_to_log(sam_to_bam.stderr, logger)
    except:
        sys.exit("Error with samtools view.\n")  


#### extract mapped or unmapped bams ######

def bam_to_mapped_or_unmapped(out_dir, combo, threads, logger):
    """ extracts mapped or unmapped reads to bam depending on the input.
    :param combo: combination of chromosome and read: chromosome_long, non_chromosome_long, non_chromosome_short or chromosome_short
	:param out_dir: output directory path
    :param threads: threads
    :param logger: logger
    :return: 
    """
    if combo == "chromosome_long":
        bam = os.path.join(out_dir, "long_read_chromosome.bam")
        map_bam = os.path.join(out_dir, "long_read_chromosome_unmapped.bam")
        map_flag = "-f"
    elif combo == "chromosome_short":
        bam = os.path.join(out_dir, "short_read_chromosome.bam")
        map_bam = os.path.join(out_dir, "short_read_chromosome_unmapped.bam")
        map_flag = "-f"
    elif combo == "non_chromosome_short":
        bam = os.path.join(out_dir, "short_read_non_chromosome.bam")
        map_bam = os.path.join(out_dir, "short_read_non_chromosome_mapped.bam")
        map_flag = "-F"
    elif combo == "non_chromosome_long":
        bam = os.path.join(out_dir, "long_read_non_chromosome.bam")
        map_bam = os.path.join(out_dir, "long_read_non_chromosome_mapped.bam")
        map_flag = "-F"
    outFile = open(map_bam, "w")
    try:
        bam_to_map = sp.Popen(["samtools", "view", "-h", "-@", threads, map_flag, "4", bam], stdout=outFile, stderr=sp.PIPE) 
        log.write_to_log(bam_to_map.stderr, logger)
    except:
        sys.exit("Error with samtools view. \n")  


#### extract fastqs ######
# long
def extract_long_fastq(out_dir, chrom, logger):
    """ extracts long mapped or unmapped reads from bam to fastq 
    :param chrom: either "chromosome" or "non_chromosome"
	:param out_dir: output directory path
    :param threads: threads
    :param logger: logger
    :return: 
    """
    # want unmapped for chromosome
    if chrom == "chromosome":
        bam = os.path.join(out_dir, "long_read_chromosome_unmapped.bam")
        long_fastq = os.path.join(out_dir, "long_read_chromosome_unmapped.fastq")
    else: # mapped for non-chromosome
        bam = os.path.join(out_dir, "long_read_non_chromosome_mapped.bam")
        long_fastq = os.path.join(out_dir, "long_read_non_chromosome_mapped.fastq")
    out_fastq = open(long_fastq, "w")
    try:
        extract_long_fastq = sp.Popen(["samtools", "bam2fq", bam], stdout=out_fastq, stderr=sp.PIPE) 
        log.write_to_log(extract_long_fastq.stderr, logger)
    except:
        sys.exit("Error with samtools bam2fastq\n")  

# short

def extract_short_fastq(out_dir, chrom, threads, logger):
    """ extracts short mapped or unmapped reads from bam to fastq 
    :param chrom: either "chromosome" or "non_chromosome"
	:param out_dir: output directory path
    :param threads: threads
    :param logger: logger
    :return: 
    """
    if chrom == "chromosome":
        bam = os.path.join(out_dir, "short_read_chromosome_unmapped.bam")
        fastq_one = os.path.join(out_dir, "unmapped_chromosome_R1.fastq.gz")
        fastq_two = os.path.join(out_dir, "unmapped_chromosome_R2.fastq.gz")
        map_flag = "-f"
    else: # non_chromosome
        bam = os.path.join(out_dir, "short_read_non_chromosome_mapped.bam")
        fastq_one = os.path.join(out_dir, "mapped_non_chromosome_R1.fastq.gz")
        fastq_two = os.path.join(out_dir, "mapped_non_chromosome_R2.fastq.gz")
        map_flag = "-F"
    try:
        extract_short_fastq= sp.Popen(["samtools", "fastq", "-@", threads, map_flag, "4", bam, "-1", fastq_one, "-2", fastq_two, "-0", "/dev/null", "-s", "/dev/null", "-n"], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(extract_short_fastq.stdout, logger)
    except:
        sys.exit("Error with samtools fastq.\n")  
