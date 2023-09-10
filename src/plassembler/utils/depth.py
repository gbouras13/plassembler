#########################################
# some code taken from & modified
# https://github.com/rrwick/Small-plasmid-Nanopore/blob/main/scripts/get_depths.py
#########################################
import os
import statistics
import subprocess as sp
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from loguru import logger

from plassembler.utils.concat import concatenate_single_fasta


def concatenate_chrom_plasmids(outdir):
    """concatenates chromosome and plasmids
    :param outdir:  Output Directory
    :param logger: logger
    :return:
    """
    chrom_fasta: Path = Path(outdir) / "chromosome.fasta"
    plas_fasta: Path = Path(outdir) / "unicycler_output/assembly.fasta"
    concat_fasta: Path = Path(outdir) / "combined.fasta"

    try:
        concatenate_single_fasta(chrom_fasta, plas_fasta, concat_fasta)
    except Exception:
        logger.error("Error with concatenate_fastas\n")


# get lengths of contigs
def get_contig_lengths(fasta):
    """gets contig lengths of combined chrom and plasmids fastas
    :param fasta:  input fasta
    :return: contig_lengths: dictionary of headers and lengths
    """
    contig_lengths = {}
    for dna_record in SeqIO.parse(fasta, "fasta"):
        plas_len = len(dna_record.seq)
        dna_header = dna_record.id
        contig_lengths[dna_header] = plas_len
    return contig_lengths


# get circular status of contigs
def get_contig_circularity(fasta):
    """gets circularity of contigs
    :param outdir:  Output Directory
    :return: circular_status: dictionary of contig header and circular status
    """
    circular_status = {}
    # add circularity
    for dna_record in SeqIO.parse(fasta, "fasta"):
        dna_header = dna_record.id
        # check if circular is in unicycler output description
        if "circular" in dna_record.description:
            circular_status[dna_header] = "circular"
        # circular chromsome
        elif "chromosome" in dna_record.id:
            circular_status[dna_header] = "circular"
        else:
            circular_status[dna_header] = "not_circular"
    return circular_status


def get_depths_from_bam(bam_file: Path, contig_lengths: pd.DataFrame):
    """maps runs samtools depth on bam
    :param bam_file: Path
    :param: contig_lengths: dictionary of headers and contig lengths
    :return: depths: dictionary of contigs and depths
    """
    depths = {}
    for repName, repLength in contig_lengths.items():
        depths[repName] = [0] * repLength
    depthCommand = ["samtools", "depth", bam_file]
    with open(os.devnull, "wb") as devNull:
        depthOutput = sp.check_output(depthCommand, stderr=devNull).decode()
    for line in depthOutput.splitlines():  # parse output
        parts = line.strip().split("\t")
        repName = parts[0]
        depths[repName][int(parts[1]) - 1] = int(parts[2])
    return depths


def collate_depths(depths, shortFlag, contig_lengths):
    """calculates summary statistics for all depths
    :param depths:  dictionary of contigs and depths from get_depths_from_bam
    :param: shortFlag: string either "short" or "long"
    :param: contig_lengths: dictionary of headers and contig lengths
    :return: summary_df: pandas df of depth summary statistics
    """
    # define the columns of dataframe
    contig_names = []
    contig_length = []
    mean_depth_col = []
    sd_depth_col = []
    q25_depth = []
    q75_depth = []
    # iterate over the conitgs
    for replicon_name, base_depths in depths.items():
        replicon_length = contig_lengths[replicon_name]
        unmaskedDepths = []
        for depth in enumerate(base_depths):
            unmaskedDepths.append(depth)
        try:
            mean_depth = round(statistics.mean(base_depths), 2)
            depth_stdev = round(statistics.stdev(base_depths), 2)
            q25, q75 = np.percentile(base_depths, [25, 75])
            q25, q75 = int(q25), int(q75)
            # save the chromosome depth
            if replicon_name == "chromosome":
                chromosome_depth = mean_depth
        except statistics.StatisticsError:  # if can't calculate
            mean_depth, depth_stdev, q25, q75 = "NA", "NA", "NA", "NA"
        # append to list
        contig_names.append(replicon_name)
        contig_length.append(replicon_length)
        mean_depth_col.append(mean_depth)
        sd_depth_col.append(depth_stdev)
        q25_depth.append(q25)
        q75_depth.append(q75)
    # make summary df
    if shortFlag == "short":
        summary_df = pd.DataFrame(
            {
                "contig": contig_names,
                "length": contig_length,
                "mean_depth_short": mean_depth_col,
                "sd_depth_short": sd_depth_col,
                "q25_depth_short": q25_depth,
                "q75_depth_short": q75_depth,
            }
        )
        summary_df["plasmid_copy_number_short"] = round(
            summary_df["mean_depth_short"] / chromosome_depth, 2
        )
    else:  # long
        summary_df = pd.DataFrame(
            {
                "contig": contig_names,
                "length": contig_length,
                "mean_depth_long": mean_depth_col,
                "sd_depth_long": sd_depth_col,
                "q25_depth_long": q25_depth,
                "q75_depth_long": q75_depth,
            }
        )
        summary_df["plasmid_copy_number_long"] = round(
            summary_df["mean_depth_long"] / chromosome_depth, 2
        )
    # return df
    return summary_df


def combine_depth_dfs(df_short, df_long, circular_status):
    """combines long and short depths
    :param outdir:  output directory
    :param df_short: short depth summary df
    :param df_long: long depth summary df
    :param: prefix: prefix - default plassembler
    :param: circular_status: dictionary of contig header and circular status
    """
    # drop double up len
    df_long = df_long.drop(["length"], axis=1)
    combined_df = pd.merge(df_short, df_long, on="contig", how="outer")
    # add in circularity info
    combined_df["circularity"] = combined_df["contig"].map(circular_status)
    return combined_df


def depth_df_single(df, circular_status):
    """final output for kmer mode
    :param outdir:  output directory
    :param df: short or long depth summary df
    :param: prefix: prefix - default plassembler
    :param: circular_status: dictionary of contig header and circular status
    """

    # add in circularity info
    df["circularity"] = df["contig"].map(circular_status)
    return df
