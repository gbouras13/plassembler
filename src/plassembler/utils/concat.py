import gzip
from pathlib import Path

from Bio import SeqIO
from loguru import logger


def concatenate_short_fastqs(out_dir):
    """moves and copies files
    :param out_dir:  Output Directory
    :param logger: logger
    :return:
    """
    # list all the inputs for concatenation
    unmapped_fastq_one_short: Path = Path(out_dir) / "unmapped_R1.fastq"
    unmapped_fastq_two_short: Path = Path(out_dir) / "unmapped_R2.fastq"
    non_chrom_fastq_one_short: Path = Path(out_dir) / "mapped_non_chromosome_R1.fastq"
    non_chrom_fastq_two_short: Path = Path(out_dir) / "mapped_non_chromosome_R2.fastq"

    # final outputs
    short_one_file: Path = Path(out_dir) / "short_read_concat_R1.fastq"
    short_two_file: Path = Path(out_dir) / "short_read_concat_R2.fastq"

    try:
        concatenate_single_fastq(
            unmapped_fastq_one_short, non_chrom_fastq_one_short, short_one_file
        )
        concatenate_single_fastq(
            unmapped_fastq_two_short, non_chrom_fastq_two_short, short_two_file
        )
    except Exception:
        logger.error("Error with concatenate_fastqs\n")


def concatenate_single_fastq(fastq_in1: Path, fastq_in2: Path, fastq_out: Path):
    """concatenates 2 fastq files
    :param fastq_in1:  fastq_in1 input fastq 1
    :param fastq_in2: fastq_in1 input fastq 2
    :param fastq_out: fastq_out output fastq 2
    :param logger: logger
    :return:
    """

    records = []

    # Read and append records from the first FASTQ file
    if fastq_in1.suffix == ".gz":
        with gzip.open(fastq_in1, "rt") as handle:
            records.extend(SeqIO.parse(handle, "fastq"))
    else:
        with open(fastq_in1, "r") as handle:
            records.extend(SeqIO.parse(handle, "fastq"))

    # Read and append records from the second FASTQ file
    if fastq_in2.suffix == ".gz":
        with gzip.open(fastq_in2, "rt") as handle:
            records.extend(SeqIO.parse(handle, "fastq"))
    else:
        with open(fastq_in2, "r") as handle:
            records.extend(SeqIO.parse(handle, "fastq"))

    # Write the concatenated records to the output FASTQ file
    with open(fastq_out, "w") as handle:
        SeqIO.write(records, handle, "fastq")


def concatenate_single_fasta(file1: Path, file2: Path, output_file: Path):
    sequences = []

    # Read sequences from the first file
    with open(file1, "r") as f1:
        sequences.extend(SeqIO.parse(f1, "fasta"))

    # Read sequences from the second file
    with open(file2, "r") as f2:
        sequences.extend(SeqIO.parse(f2, "fasta"))

    # Write concatenated sequences to the output file
    with open(output_file, "w") as output:
        SeqIO.write(sequences, output, "fasta")
