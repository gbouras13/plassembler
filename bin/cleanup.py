import os
import subprocess as sp
from Bio import SeqIO
import pandas as pd
import sys
import glob
from Bio.SeqRecord import SeqRecord
import shutil


####################################################
# cleanup
##########################################################

def remove_intermediate_files(out_dir, keep_chromosome, assembled_mode, long_only):
    """ removes intermediate files
    :param out_dir:  Output Directory
    :return: 
    """

    # find all files with the suffix "fastq"
    # find all files with the specified suffixes
    suffixes = ['fastq', 'bam', 'sa', 'sam', 'json', 'bed', 'msh' ]
    files = []
    for suffix in suffixes:
        files.extend(glob.glob(os.path.join(out_dir, "*." + suffix)))

    # loop through the files and remove them
    for file in files:
        remove_file(file)
    
    if assembled_mode == False:
        shutil.rmtree(os.path.join(out_dir,"00-assembly"))
        shutil.rmtree(os.path.join(out_dir,"10-consensus"))
        shutil.rmtree(os.path.join(out_dir,"20-repeat"))
        shutil.rmtree(os.path.join(out_dir,"30-contigger"))
        shutil.rmtree(os.path.join(out_dir, "40-polishing"))

    if long_only == True:
        shutil.rmtree(os.path.join(out_dir,"unicycler_output"))

    # delete intermediate mash file
    remove_file(os.path.join(out_dir,"mash.tsv") )

    # delete intermediate fasta assemble files
    remove_file(os.path.join(out_dir,"combined.fasta"))
    remove_file(os.path.join(out_dir,"flye_renamed.fasta"))
    remove_file(os.path.join(out_dir,"plasmids.fasta"))

    # delete fastq intermediate files
    remove_file(os.path.join(out_dir,"final_filtered_long_reads.fastq.gz"))
    remove_file(os.path.join(out_dir,"chopper_long_reads.fastq.gz"))
    remove_file(os.path.join(out_dir, "multimap_plasmid_chromosome_long.fastq"))

    # multimer
    remove_file(os.path.join(out_dir,"mapping.paf"))

    # chromosome
    if keep_chromosome == False:
        remove_file(os.path.join(out_dir,"chromosome.fasta"))


def move_and_copy_files(out_dir, prefix, unicycler_success_flag, keep_fastqs, assembled_mode):
    """ moves and copies files
    :param out_dir:  Output Directory
    :param prefix: prefix
    :param unicycler_success_flag: whether or not unicycler worked
    :return: 
    """
    # make flye dir 
    flye_dir = os.path.join(out_dir,"flye_output")
    if not os.path.exists(flye_dir):
        os.mkdir(flye_dir)

    # move flye files
    if assembled_mode == False:
        shutil.move(os.path.join(out_dir,"assembly.fasta"), flye_dir) 
        shutil.move(os.path.join(out_dir,"assembly_info.txt"), flye_dir)
        shutil.move(os.path.join(out_dir,"flye.log"), flye_dir)
        shutil.move(os.path.join(out_dir,"assembly_graph.gfa"), flye_dir)
        shutil.move(os.path.join(out_dir,"assembly_graph.gv"), flye_dir)

    if unicycler_success_flag == True:
         # move unicycler graph output to main directory
        shutil.copy2( os.path.join(out_dir,"unicycler_output", "assembly.gfa"), os.path.join(out_dir, prefix + "_plasmids.gfa"))
    else:
        # to touch empty versions of the output files if no plasmids 
        touch_output_fail_files(out_dir, prefix)

    # put kept fastqs into separate directory
    # make fastqs dir 
    if keep_fastqs == True:
        fastqs_dir = os.path.join(out_dir,"plasmid_fastqs")
        if not os.path.exists(fastqs_dir):
            os.mkdir(fastqs_dir)

        # move flye files
        shutil.move(os.path.join(out_dir,"short_read_concat_R1.fastq"), os.path.join(fastqs_dir,"plasmids_R1.fastq")) 
        shutil.move(os.path.join(out_dir,"short_read_concat_R2.fastq"), os.path.join(fastqs_dir,"plasmids_R2.fastq")) 
        shutil.move(os.path.join(out_dir,"plasmid_long.fastq"), os.path.join(fastqs_dir,"plasmids_long.fastq")) 
        shutil.move(os.path.join(out_dir, "multimap_plasmid_chromosome_long.fastq"), os.path.join(out_dir, "multimap_long.fastq")) 


# function to touch create a file 
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch_file(path):
    with open(path, 'a'):
        os.utime(path, None)

# to create empty plasmids fasta and gfa files
def touch_output_fail_files(out_dir, prefix):
    touch_file(os.path.join(out_dir, prefix + "_plasmids.fasta"))
    touch_file(os.path.join(out_dir, prefix + "_plasmids.gfa"))
    touch_file(os.path.join(out_dir, prefix + "_summary.tsv"))


def remove_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)