import os
import sys
import subprocess as sp
import plassemblerModules

def concatenate_all_fastqs(out_dir,logger):
    """ moves and copies files
    :param out_dir:  Output Directory
    :param logger: logger
    :return: 
    """
    # list all the inputs for concatenation
    long_fastq_chrom = os.path.join(out_dir, "long_read_chromosome_unmapped.fastq")
    long_fastq_non_chrom = os.path.join(out_dir, "long_read_non_chromosome_mapped.fastq")
    chrom_fastq_one_short = os.path.join(out_dir, "unmapped_chromosome_R1.fastq.gz")
    chrom_fastq_two_short = os.path.join(out_dir, "unmapped_chromosome_R2.fastq.gz")
    non_chrom_fastq_one_short = os.path.join(out_dir, "mapped_non_chromosome_R1.fastq.gz")
    non_chrom_fastq_two_short = os.path.join(out_dir, "mapped_non_chromosome_R2.fastq.gz")
    # final outputs
    long_file = open(os.path.join(out_dir, "long_read_concat.fastq"), "w")
    short_one_file = open(os.path.join(out_dir, "short_read_concat_R1.fastq.gz"), "w")
    short_two_file = open(os.path.join(out_dir, "short_read_concat_R2.fastq.gz"), "w")
    try:
        concatenate_single(long_fastq_chrom, long_fastq_non_chrom, long_file, logger)
        concatenate_single(chrom_fastq_one_short, non_chrom_fastq_one_short, short_one_file, logger)
        concatenate_single(chrom_fastq_two_short, non_chrom_fastq_two_short, short_two_file, logger)
    except:
        sys.exit("Error with concatenate_fastqs\n")  
    
        
def concatenate_single(fastq_in1, fastq_in2, fastq_out, logger):
    """ concatenates 2 fastq files
    :param fastq_in1:  fastq_in1 input fastq 1
    :param fastq_in2: fastq_in1 input fastq 2
    :param fastq_out: fastq_out output fastq 2
    :param logger: logger 
    :return: 
    """
    concat_fastq = sp.Popen(["cat", fastq_in1, fastq_in2 ], stdout=fastq_out, stderr=sp.PIPE)
    plassemblerModules.write_to_log(concat_fastq.stderr, logger) 

