from Bio import SeqIO
import os
import subprocess as sp
import log
import sys
import shutil



def mash_sketch(out_dir, fasta_file,  logger):
    """
    Runs mash to output fastas
    :param out_dir: output directory
    :param logger: logger
    :return:
    """
    plasmid_fasta = os.path.join(out_dir, "plasmids.fasta")
    shutil.copy2(fasta_file, plasmid_fasta)
    mash_sketch = sp.Popen(["mash", "sketch",  plasmid_fasta, "-i" ], stdout=sp.PIPE, stderr=sp.PIPE) 
    log.write_to_log(mash_sketch.stdout, logger)
    try:
        mash_sketch = sp.Popen(["mash", "sketch",  plasmid_fasta, "-i" ], stdout=sp.PIPE, stderr=sp.PIPE) 
        log.write_to_log(mash_sketch.stdout, logger)
    except:
        sys.exit("Error with mash sketch.\n")  


def run_mash(out_dir,  plassembler_db_dir, logger):
    """
    Runs mash to output fastas
    :param out_dir: output directory
    :param plassembler_db_dir: plassembler db directory
    :param logger: logger
    :return:
    """

    plsdb_sketch = os.path.join(plassembler_db_dir, "plsdb.msh")
    plasmid_sketch = os.path.join(out_dir, "plasmids.fasta.msh")
    mash_tsv = os.path.join(out_dir,"mash.tsv")
    outFile = open(mash_tsv, "w")

    try:
        mash_sketch = sp.Popen(["mash", "dist",  plasmid_sketch, plsdb_sketch, "-v", "0.1", "-d", "0.1", "-i" ], stdout=outFile, stderr=sp.PIPE) 
        log.write_to_log(mash_sketch.stderr, logger)
    except:
        sys.exit("Error with mash dist.\n")  

def get_contig_count( plasmid_fasta):
    """
    Process mash output
    :param out_dir: output directory
    :return: i: int contig_count
    """
    i = 0
    for dna_record in SeqIO.parse(plasmid_fasta, 'fasta'): 
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



