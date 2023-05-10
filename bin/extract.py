import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


def extract_chromosome_fasta(out_dir, contig_name):
    """Extracts Chromosome fasta

    :param out_dir: output directory
    :param contig_name: name of the chromosome contig
    :return: 
    """
    with open(os.path.join(out_dir, "chromosome.fasta"), 'w') as chrom_fa:
        for dna_record in SeqIO.parse(os.path.join(out_dir, "assembly.fasta"), 'fasta'): 
            # has to be contig of the chromosome
            if dna_record.id == contig_name:
                dna_header = "chromosome"
                dna_description = ""
                dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                SeqIO.write(dna_record, chrom_fa, 'fasta')


def extract_plasmid_fastas(out_dir, contig_name):
    """Extracts plasmid fastas

    :param out_dir: output directory
    :param contig_name: name of the chromosome contig
    :return: 
    """
    with open(os.path.join(out_dir, "non_chromosome.fasta"), 'w') as non_chrom_fa:
        for dna_record in SeqIO.parse(os.path.join(out_dir, "assembly.fasta"), 'fasta'): 
            # has to be contig of the chromosome
            if dna_record.id != contig_name:
                dna_header = dna_record.id
                dna_description = ""
                dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
                SeqIO.write(dna_record, non_chrom_fa, 'fasta')

