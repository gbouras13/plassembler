import os
import sys
import subprocess as sp
import gzip
from Bio import SeqIO


def concatenate_short_fastqs(out_dir):
    """moves and copies files
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
        concatenate_single(
            unmapped_fastq_one_short, non_chrom_fastq_one_short, short_one_file
        )
        concatenate_single(
            unmapped_fastq_two_short, non_chrom_fastq_two_short, short_two_file
        )
    except:
        sys.exit("Error with concatenate_fastqs\n")


def concatenate_single(fastq_in1, fastq_in2, fastq_out):
    """concatenates 2 fastq files
    :param fastq_in1:  fastq_in1 input fastq 1
    :param fastq_in2: fastq_in1 input fastq 2
    :param fastq_out: fastq_out output fastq 2
    :param logger: logger
    :return:
    """
    
    records = []

    # Read and append records from the first FASTQ file
    if fastq_in1.endswith(".gz"):
        with gzip.open(fastq_in1, "rt") as handle:
            records.extend(SeqIO.parse(handle, "fastq"))
    else:
        with open(fastq_in1, "r") as handle:
            records.extend(SeqIO.parse(handle, "fastq"))

    # Read and append records from the second FASTQ file
    if fastq_in2.endswith(".gz"):
        with gzip.open(fastq_in2, "rt") as handle:
            records.extend(SeqIO.parse(handle, "fastq"))
    else:
        with open(fastq_in2, "r") as handle:
            records.extend(SeqIO.parse(handle, "fastq"))

    # Write the concatenated records to the output FASTQ file
    with open(fastq_out, "w") as handle:
        SeqIO.write(records, handle, "fastq")

