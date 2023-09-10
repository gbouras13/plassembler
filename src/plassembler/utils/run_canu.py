import math
from collections import Counter
from itertools import product
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
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


def make_blastdb(canu_output_dir, plasmid_fasta, logdir):
    db: Path = Path(canu_output_dir) / "db"

    makeblastdb = ExternalTool(
        tool="makeblastdb",
        input=f"-in {plasmid_fasta}",
        output=f"-out {db}",
        params="-dbtype nucl ",
        logdir=logdir,
        outfile="",
    )

    ExternalTool.run_tool(makeblastdb, to_stdout=False)


"""
some of this adapted from dnaapler 

https://github.com/gbouras13/dnaapler

"""


def run_blast(canu_output_dir, plasmid_fasta, threads, logdir):
    blast_output: Path = Path(canu_output_dir) / "blast_output.txt"
    db: Path = Path(canu_output_dir) / "db"

    blast = ExternalTool(
        tool="blastn",
        input=f"-query {plasmid_fasta}",
        output=f"-out {blast_output}",
        params=f'-db {db} -evalue  1e-05 -num_threads {threads} -outfmt " 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq "',
        logdir=logdir,
        outfile="",
    )
    ExternalTool.run_tool(blast, to_stdout=False)


def process_blast_output(canu_output_dir, combined_plasmid_file, outdir):
    """Processes

    :param input: input file
    :param blast_file: blast output file
    :param out_file: output file
    :return: output_filename: Path of output dedeup Fasta
    """

    blast_output_file: Path = Path(canu_output_dir) / "blast_output.txt"

    # define colnames
    col_list = [
        "qseqid",
        "qlen",
        "sseqid",
        "slen",
        "length",
        "qstart",
        "qend",
        "sstart",
        "send",
        "pident",
        "nident",
        "gaps",
        "mismatch",
        "evalue",
        "bitscore",
        "qseq",
        "sseq",
    ]

    # read in the dataframe from BLAST
    try:
        blast_df = pd.read_csv(
            blast_output_file, delimiter="\t", index_col=False, names=col_list
        )
    except Exception:
        logger.error("There was an issue with parsing the BLAST output file.")

    # if the BLAST input is empty
    if isinstance(blast_df, pd.DataFrame) and blast_df.empty:
        logger.error("There were 0 BLAST hits. This must be a BLAST error.")

    # keep only the columns where the same contig is blasted against each other

    # Filter rows where qseqid is equal to sseqid - the ones fo BLASTing against themselves
    blast_df = blast_df[blast_df["qseqid"] == blast_df["sseqid"]]

    # read in the canu FASTA and save as a dictionary

    fasta_dict = {}
    i = 1
    for record in SeqIO.parse(combined_plasmid_file, "fasta"):
        fasta_dict[record.id] = {
            "count": i,
            "sequence": str(record.seq),
            "dupe": False,  # stores the dupe status
            "start": 1,  # stores the start
            "end": len(record.seq),  # stores the end
        }
        i += 1

    for contig in fasta_dict.keys():
        tmp_df = blast_df[blast_df["qseqid"] == contig]
        # Sort by 'length' column in descending order
        tmp_df_sorted = tmp_df.sort_values(by="length", ascending=False)
        # get rid of the 100% match row
        tmp_df_sorted = tmp_df_sorted[tmp_df_sorted["qlen"] != tmp_df_sorted["length"]]
        # more than 99% identical
        tmp_df_sorted = tmp_df_sorted[tmp_df_sorted["pident"] > 99]
        # starts need to be < 100
        tmp_df_sorted = tmp_df_sorted[tmp_df_sorted["qstart"] < 100]
        num_rows = tmp_df_sorted.shape[0]
        if num_rows == 0:  # where there is no dupe at all
            fasta_dict[contig]["dupe"] = False
            # exit
        else:
            # Get the first row
            first_row = tmp_df_sorted.iloc[0]
            # ensure the match is good
            if (
                first_row["length"] < 500
            ):  # less than 500bp repeat in the top hit - probably Insertion Seq not a real plasmid dupe - otherwise probably legit
                fasta_dict[contig]["dupe"] = False
            else:
                try:
                    # the repeat will be in the longest hit with qstart < 100  (usually 1 or very close to it)
                    # heuristic i need to check i guess
                    # just take until the next repeat element
                    # this has been filtered for prior
                    best_row = tmp_df_sorted.iloc[0]
                    fasta_dict[contig]["dupe"] = True
                    fasta_dict[contig]["start"] = best_row["qstart"]
                    fasta_dict[contig]["end"] = best_row["sstart"]

                    # if the query end is larger than the sstart - there is an overlap
                    # take 1 as the start and then the sstart as the end
                    # otherwise check for  concatenation (within 1000bp))
                    # otherwise exit just the whole plasmid
                    # if best_row["qend"] > best_row["sstart"]:
                    #     fasta_dict[contig]["dupe"] = True
                    #     fasta_dict[contig]["start"] = best_row["qstart"]
                    #     fasta_dict[contig]["end"] = best_row["sstart"]
                    # #elif (best_row["qend"] + 1000) > best_row[
                    #     #"sstart"
                    # #]:  # the longest match is likely to be a duplication
                    # else:
                    #     fasta_dict[contig]["dupe"] = True
                    #     fasta_dict[contig]["start"] = best_row["qstart"]
                    #     fasta_dict[contig]["end"] = best_row["sstart"]
                    # else:
                    #     fasta_dict[contig]["dupe"] = False
                except Exception:
                    logger.error("Flye not found. Please reinstall Plassembler.")

    # Create a list of SeqRecord objects
    records = []
    for entry_id, entry_data in fasta_dict.items():
        subsequence = entry_data["sequence"][
            entry_data["start"] - 1 : entry_data["end"]
        ]
        count = str(entry_data["count"])
        l = entry_data["end"]

        record = SeqRecord(
            seq=Seq(subsequence), id=count, description=f"{count} len={l}"
        )
        records.append(record)

    # Write the records to a FASTA file
    output_filename: Path = Path(outdir) / "combined_plasmids_dedup.fasta"
    with open(output_filename, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    return output_filename


# then BLAST output
# # then figure out for overlaps
# parse all blast hits as a dictionary
# need more than 1 hit (itself)
# if the blast hit is more than 50% and less than 90 % of contig length, take as duplicate
# elif the blast hit is more than 2000bp (lower could be e.g. IS element) and far away (more than 50% length away)
# then assume partial duplication too
# get all dupe regions this way
