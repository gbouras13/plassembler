import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd



def extract_chromosome(out_dir, chromosome_len, no_plasmids_flag):
    """Extracts Chromosome

    :param out_dir: output directory
    :param chromosome_len: lower bound on length of chromosome
    :param no_plasmids_flag: boolean - false means yes, there is a plasmid
    :return: correct_chromosome = a flag saying whether or not there was an identified chromosome
    """
    info_file =  os.path.join(out_dir, "assembly_info.txt")
    col_list = ["seq_name", "length", "cov", "circ", "repeat", "mult", "alt_group", "graph_path"] 
    info_df = pd.read_csv(info_file, delimiter= '\t', index_col=False , names=col_list, skiprows=1) 
    max_length = max(info_df['length'])
    # get putative chromosome contig
    chrom_contig = info_df[info_df['length'] == max_length].iloc[0]['seq_name']
    # chrom_circ = info_df[info_df['length'] == max_length].iloc[0]['circ']

    # of chromosome is at least 90% of inputted chromosome length, flag 
    # to say that the chromosome has been correctly identified 
    correct_chromosome = True
    if max_length < int(chromosome_len)*0.9:  
        correct_chromosome = False
    else:
        extract_chromosome_fasta(out_dir, chrom_contig)
        # only extract the plasmids if more than 1 contig was assembled
        if no_plasmids_flag == False:
            extract_plasmid_fastas(out_dir, chrom_contig)
    return correct_chromosome

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

