import os
import shutil
from pathlib import Path

from Bio import SeqIO

from plassembler.utils.external_tools import ExternalTool


def mash_sketch(out_dir, fasta_file, logdir):
    """
    Runs mash to output fastas
    :param out_dir: output directory
    :param logger: logger
    :return:
    """

    plasmid_fasta: Path = Path(f"{out_dir}/plasmids.fasta")
    shutil.copy2(fasta_file, plasmid_fasta)

    # mash command
    mash = ExternalTool(
        tool="mash",
        input="",
        output="",
        params=f" sketch {plasmid_fasta} -i ",
        logdir=logdir,
        outfile="",
    )

    # need to write to stdout
    ExternalTool.run_tool(mash, to_stdout=False)


def run_mash(out_dir, plassembler_db_dir, logdir):
    """
    Runs mash to output fastas
    :param out_dir: output directory
    :param plassembler_db_dir: plassembler db directory
    :param logger: logger
    :return:
    """

    plsdb_sketch: Path = Path(f"{plassembler_db_dir}/plsdb_2023_11_03_v2.msh")
    plasmid_sketch: Path = Path(f"{out_dir}/plasmids.fasta.msh")
    mash_tsv: Path = Path(f"{out_dir}/mash.tsv")

    mash = ExternalTool(
        tool="mash",
        input="",
        output="",
        params=f" dist  {plasmid_sketch} {plsdb_sketch} -v 0.1 -d 0.1 -i ",
        logdir=logdir,
        outfile=mash_tsv,
    )

    # need to write to stdout
    ExternalTool.run_tool(mash, to_stdout=True)


def get_contig_count(plasmid_fasta):
    """
    Process mash output
    :param out_dir: output directory
    :return: i: int contig_count
    """
    i = 0
    for dna_record in SeqIO.parse(plasmid_fasta, "fasta"):
        i += 1
    return i


# check if a file has more than 1 line (not empty)
def is_file_empty(file):
    """
    Determines if file is empty
    :param file: file path
    :return: empty Boolean
    """
    empty = False
    if os.stat(file).st_size == 0:
        empty = True
    return empty
