from src.external_tools import ExternalTool
from pathlib import Path


#################################
# original mapping
#################################


def minimap_long_reads(input_long_reads,fasta,sam, threads, pacbio_model, logdir):
    """maps long reads using minimap2
    :param threads: threads
    :param pacbio_model: pacbio_model
    :param threads: threads
    :param logdir: logdir
    :return:
    """

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


def minimap_short_reads(r1, r2, fasta, sam, threads, logdir):
    """maps short reads using minimap2
    :param outdir: output directory path
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    minimap2 = ExternalTool(
        tool="minimap2",
        input=f"",
        output=f"",
        params=f" -ax sr -t {threads} {fasta} {r1} {r2}",
        logdir=logdir,
        outfile = sam
    )

    # need to write to stdout
    ExternalTool.run_tool(minimap2, to_stdout = True)
