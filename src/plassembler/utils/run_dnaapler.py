from pathlib import Path

from loguru import logger

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

    try:
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
            Path(outdir) / "dnaapler" / "dnaapler_all_reoriented.fasta"
        )
        return plasmids_for_sketching
    except Exception:
        logger.warning("Dnaapler failed to reorient any plasmids.")
        plasmids_for_sketching = plasmid_fasta
        return plasmids_for_sketching
