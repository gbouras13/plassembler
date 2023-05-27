import os
import sys
from src.external_tools import ExternalTool
from pathlib import Path


#################################
# original mapping
#################################


def minimap_long_reads(outdir, threads, pacbio_model, logdir):
    """maps long reads using minimap2
    :param threads: threads
    :param pacbio_model: pacbio_model
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    input_long_reads: Path =  outdir/ f"chopper_long_reads.fastq.gz"
    fasta: Path =  outdir/ f"flye_renamed.fasta"
    sam: Path = outdir/ f"long_read.sam"

    # ONT
    minimap2_model = "map-ont"

    # Pacbio
    if pacbio_model == "--pacbio-raw" or "--pacbio-corr":
        minimap2_model = "map-pb"

    if pacbio_model == "--pacbio-hifi":
        minimap2_model = "map-hifi"


    minimap2 = ExternalTool(
        tool="minimap2",
        input=f"",
        output=f"",
        params=f" -ax {minimap2_model} -t {threads} {fasta} {input_long_reads}",
        logdir=logdir,
        outfile = sam
    )

    # need to write to stdout
    ExternalTool.run_tool(minimap2, to_stdout = True)

# short reads


def minimap_short_reads(outdir, threads, logdir):
    """maps short reads using minimap2
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    trim_one: Path =  outdir/ f"trimmed_R1.fastq"
    trim_two: Path =  outdir/ f"trimmed_R2.fastq"
    fasta: Path =  outdir/ f"flye_renamed.fasta"
    sam: Path = outdir/ f"short_read.sam"

    minimap2 = ExternalTool(
        tool="minimap2",
        input=f"",
        output=f"",
        params=f" -ax sr -t {threads} {fasta} {trim_one} {trim_two}",
        logdir=logdir,
        outfile = sam
    )

    # need to write to stdout
    ExternalTool.run_tool(minimap2, to_stdout = True)
