from pathlib import Path

from plassembler.utils.external_tools import ExternalTool


def run_flye(outdir, threads, raw_flag, pacbio_model, logdir):
    """Runs flye on trimmed long reads

    :param outdir: output directory
    :param raw_flag: boolean - true if --nano-raw used for flue
    :param logger: logger
    :return:
    """
    trim_long: Path = Path(outdir) / "chopper_long_reads.fastq.gz"
    flye_model = "--nano-hq"
    if raw_flag is True:
        flye_model = "--nano-raw"
    if pacbio_model != "nothing":
        flye_model = pacbio_model

    flye = ExternalTool(
        tool="flye",
        input=" ",
        output=" ",
        params=f"{flye_model} {trim_long} --out-dir {outdir} --threads {threads} ",
        logdir=logdir,
        outfile="",
    )

    ExternalTool.run_tool(flye)


def run_raven(outdir, threads, logdir):
    """Runs raven on trimmed long reads

    :param outdir: output directory
    :param raw_flag: boolean - true if --nano-raw used for flye
    :param logger: logger
    :return:
    """
    trim_long: Path = Path(outdir) / "chopper_long_reads.fastq.gz"

    # gfa
    gfa: Path = Path(outdir) / "assembly_graph.gfa"

    # define outfile
    outfile: Path = Path(outdir) / "assembly.fasta"

    raven = ExternalTool(
        tool="raven",
        input="",
        output="",
        params=f" -t {str(threads)} {trim_long} --graphical-fragment-assembly {gfa}",
        logdir=logdir,
        outfile=outfile,
    )

    # need to write to stdout
    ExternalTool.run_tool(raven, to_stdout=True)
