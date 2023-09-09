from pathlib import Path

from loguru import logger

from plassembler.utils.external_tools import ExternalTool


def run_dnaapler(threads, logdir, outdir):
    """runs dnaapler bulk
    :param long: long read fastq
    :param canu_output_dir: canu Output Directory
    :param threads: threads
    :param logdir: logdir
    :return:
    """

    canu_plasmid_fasta: Path = Path(outdir) / "plasmids_canu.fasta"
    dnaapler_outdir: Path = Path(outdir) / "dnaapler"

    try:
        dnaapler = ExternalTool(
            tool="dnaapler all",
            input="",
            output="",
            params=f" -i {canu_plasmid_fasta} -o {dnaapler_outdir} -t {threads}",
            logdir=logdir,
            outfile="",
        )

        ExternalTool.run_tool(dnaapler, to_stdout=False)
    except Exception:
        logger.error(f"Dnaapler failed.")
