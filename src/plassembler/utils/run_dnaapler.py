import os
from pathlib import Path

from plassembler.utils.external_tools import ExternalTool


def run_dnaapler(threads, plasmid_fasta, logdir, outdir):
    """runs dnaapler bulk
    :param long: long read fastq
    :param plasmid_fasta
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    dnaapler_outdir: Path = Path(outdir) / "dnaapler"

    dnaapler = ExternalTool(
        tool="dnaapler all",
        input="",
        output="",
        params=f" -i {plasmid_fasta} -o {dnaapler_outdir} -t {threads}",
        logdir=logdir,
        outfile="",
    )

    ExternalTool.run_tool(dnaapler, to_stdout=False)

    plasmids_for_sketching: Path = (
        Path(outdir) / "dnaapler" / "dnaapler_reoriented.fasta"
    )

    # if dnaapler failed to reorient anything - then return the original plasmid fasta
    if os.path.exists(plasmids_for_sketching) is False:
        plasmids_for_sketching = plasmid_fasta

    return plasmids_for_sketching
