from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def create_fake_flye_chromosome_assembly(assembly_fasta_file: Path) -> None:
    """
    creates 3MB chromosome of A's as the fake chromosome
    """

    # Create a sequence of 5 million 'A's
    sequence = Seq("A" * 3000000)

    # Create a SeqRecord for the sequence
    record = SeqRecord(sequence, id="chromosome", description="Chromosome Sequence")

    # Write the SeqRecord to a FASTA file
    with open(assembly_fasta_file, "w") as output_file:
        SeqIO.write(record, output_file, "fasta")


def create_fake_flye_chromosome_info(assembly_info_file: Path) -> None:
    """
    creates 3MB chromosome flye info
    """

    data = {
        "seq_name": ["chromosome"],
        "length": [3000000],
        "cov.": [60],
        "circ.": ["Y"],
        "repeat": ["N"],
        "mult.": [1],
        "alt_group": ["*"],
        "graph_path": [2],
    }

    flye_info_df = pd.DataFrame(data)
    flye_info_df.to_csv(assembly_info_file, index=False)
