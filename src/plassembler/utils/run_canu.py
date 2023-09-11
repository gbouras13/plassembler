import math
from itertools import product
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger

from plassembler.utils.external_tools import ExternalTool


def run_canu(
    threads,
    logdir,
    longreads,
    canu_output_dir,
    canu_nano_or_pacbio,
    total_flye_plasmid_length,
):
    """runs canu
    :param long: long read fastq
    :param canu_output_dir: canu Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """
    # canu -p C308_canu -d C308_canu genomeSize=0.01m maxInputCoverage=200 maxThreads=8 -nanopore plasmid_long.fastq
    # for the assembly param need to divide by a million
    total_flye_plasmid_length = round(total_flye_plasmid_length / 1000000, 5)
    try:
        canu = ExternalTool(
            tool="canu",
            input="",
            output="",
            params=f" -p canu -d {canu_output_dir} genomeSize={total_flye_plasmid_length}m maxInputCoverage=250 stopOnLowCoverage=1 maxThreads={threads} -{canu_nano_or_pacbio} {longreads}",
            logdir=logdir,
            outfile="",
        )

        ExternalTool.run_tool(canu, to_stdout=False)
    except Exception:
        logger.info(
            f"canu failed. This likely means that you have non-chromosomal reads to assemble anything. It is likely that you have no plasmids in this sample, but check the canu output in {logdir}."
        )


"""
5mer entropy
"""


def shannon_entropy_5mers(dna_sequence):
    """gets 5mer shannon entropy
    :param dna_sequence - str of DNA sequence
    :return: float entropy
    """

    # nucleotides and kmers
    nucleotides = "ACGT"
    k = 4

    # Generate all k-mer combinations
    five_mers_all = ["".join(p) for p in product(nucleotides, repeat=k)]

    kmer_counts = {}
    # total number of k-mers
    total_kmers = 0

    for kmer in five_mers_all:
        count = dna_sequence.count(kmer)
        kmer_counts[kmer] = count
        total_kmers += count

    # Calculate the probability of each k-mer
    probabilities = [
        int(count) / int(total_kmers) for count in kmer_counts.values() if count > 0
    ]

    # Calculate the Shannon entropy using the formula
    entropy = -sum(p * math.log2(p) for p in probabilities)

    # max entropy is log2(4^k). So in this case it is log2(4^4) = log2(256) = 8

    return entropy


"""
determine entropy
"""


def filter_entropy(canu_fasta, outdir):
    """runs canu
    :param long: long read fastq
    :param canu_output_dir: canu Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    filtered_records = []

    for record in SeqIO.parse(canu_fasta, "fasta"):
        entropy = shannon_entropy_5mers(str(record.seq))
        if (
            entropy > 5
        ):  # reject low entropy contigs - any real sequence should have entropy more than this (close to 8)
            filtered_records.append(record)

    # Write the records to a FASTA file
    output_filename: Path = Path(outdir) / "canu_filtered_contigs.fasta"
    with open(output_filename, "w") as output_handle:
        SeqIO.write(filtered_records, output_handle, "fasta")

    return output_filename


"""
trim contigs
also in trycycler fwiw
https://github.com/rrwick/Trycycler/blob/main/scripts/canu_trim.py
"""


def trim_contigs(canu_filtered_fasta, outdir):
    """trims contigs based on canu trim designation
    :param long: long read fastq
    :param canu_output_dir: canu Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    trimmed_records = []

    i = 1
    for record in SeqIO.parse(canu_filtered_fasta, "fasta"):
        # Extract the header
        header = record.description

        # Check if 'trim=' is present in the header
        if "trim=" in header and "suggestCircular=yes" in header:
            # Extract the coordinates from the header
            trim_info = header.split("trim=")[1]
            trim_coords = trim_info.split("-")

            # Parse the start and end coordinates
            start_coord = int(trim_coords[0])
            end_coord = int(trim_coords[1])

        circular_info = header.split("suggestCircular=")[1]
        circular_info = circular_info.split(" ")[0]

        if circular_info == "yes":
            desc = "circular=True"
        else:
            desc = ""

        # Extract the sequence
        subsequence = record.seq[start_coord:end_coord]
        record_trim = SeqRecord(seq=subsequence, id=str(i), description=desc)
        trimmed_records.append(record_trim)
        i += 1

    # Write the records to a FASTA file
    output_filename: Path = Path(outdir) / "canu_filtered_trimmed_contigs.fasta"
    with open(output_filename, "w") as output_handle:
        SeqIO.write(trimmed_records, output_handle, "fasta")

    return output_filename
