import os
import sys
import subprocess as sp
import log

def concatenate_all_fastqs(out_dir,logger):
    """ moves and copies files
    :param out_dir:  Output Directory
    :param logger: logger
    :return: 
    """
    # list all the inputs for concatenation
    unmapped_fastq_one_short = os.path.join(out_dir, "unmapped_R1.fastq")
    unmapped_fastq_two_short = os.path.join(out_dir, "unmapped_R2.fastq")
    non_chrom_fastq_one_short = os.path.join(out_dir, "mapped_non_chromosome_R1.fastq")
    non_chrom_fastq_two_short = os.path.join(out_dir, "mapped_non_chromosome_R2.fastq")

    # final outputs
    short_one_file = open(os.path.join(out_dir, "short_read_concat_R1.fastq"), "w")
    short_two_file = open(os.path.join(out_dir, "short_read_concat_R2.fastq"), "w")

    try:
        concatenate_single(unmapped_fastq_one_short, non_chrom_fastq_one_short, short_one_file, logger)
        concatenate_single(unmapped_fastq_two_short, non_chrom_fastq_two_short, short_two_file, logger)
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
    log.write_to_log(concat_fastq.stderr, logger) 


    